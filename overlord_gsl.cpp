#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <gsl/gsl_spline.h>
#include <slepc.h>
#include <mpi.h>
#include "mkl.h"
#include "utils.h"
#include "sasha.h"
#include "bieleclist.h"
#include "class_grid.h"
#include "orb.h"
int main(int argc,char **argv){
  FILE *file;
  int i,j,k,norb,nirrep,l,tag=99;
  int id=0,numprocs,root=0,nbas,lamax,lmax,io,*lo=NULL,*mo=NULL,lmorb;
  int *la,*ma,*le,*me,*no,lmerb=0;
  double z1,z2,z3,z4,start_time,end_time,Vnuc;
  cmplx *Xi,*Wi,*TXX,*TIXX,*c_nl,*e_nl,*cilm=NULL,*K,*J,*NE,*CKj,*CKi,*CJj,*CJi;
  char command[150],fichero[50],fileel[50]="angular_";
  char symm[5]="ylm",**irrep,str[5]="y",**irreps,strs[5]="c";
  orbitals *od;
  hamilton *hd;

  int lwork,info=0,nb,n,nf=5,mf=5,*nmkl,*ini;
  double *AT,norm;
  char solution[50]="bound.m";
  
  Mat AA,BB;
  Vec x,y,*Q;
  PetscViewer view;
  EPS eps;
  PetscReal error;
  PetscScalar cval,kr,ki,vAA;
  PetscInt d_nz,o_nz;
  int nvec,nloc,its,nconv,start,end,il,jl,i1,j1,type;
  
  SlepcInitialize(&argc,&argv,NULL,NULL);  
  MPI_Comm_size(PETSC_COMM_WORLD,&numprocs);
  MPI_Comm_rank(PETSC_COMM_WORLD,&id);

  file=fopen("quadrature.dat","r");
  fscanf(file,"%d",&lamax);
  fscanf(file,"%d",&nbas);
  Xi=(cmplx*)calloc(nbas,sizeof(cmplx));
  Wi=(cmplx*)calloc(nbas,sizeof(cmplx));
  for(i=0;i<nbas;i++){
    fscanf(file,"%le%le%le%le",&z1,&z2,&z3,&z4);
    Xi[i]=z1+z2*I;
    Wi[i]=z3+z4*I;
  }
  fclose(file);

  lmax=floor(lamax/2);
  TXX=(cmplx*)calloc(nbas*nbas*lamax,sizeof(cmplx));
  TIXX=(cmplx*)calloc(nbas*nbas*lamax,sizeof(cmplx));
  file=fopen("kinetic.dat","r");
  for(i=0;i<nbas*nbas*lamax;i++){
    fscanf(file,"%le%le",&z1,&z2);
    TXX[i]=z1+z2*I;
  }
  for(i=0;i<nbas*nbas*lamax;i++){
    fscanf(file,"%le%le",&z1,&z2);
    TIXX[i]=z1+z2*I;
  }
  fclose(file);  
  GRID CORR(nbas,lamax,Xi,Wi,TXX,TIXX);
  free(Xi);
  free(Wi);
  free(TXX);
  free(TIXX);

  file=fopen("orb.inp","r");
  fscanf(file,"%s",fichero);
  fscanf(file,"%s%d",symm,&nirrep);
  irrep=(char**)calloc(nirrep,sizeof(char*));
  for(i=0;i<nirrep;i++){
    irrep[i]=(char*)calloc(5,sizeof(char));
  }
  for(i=0;i<nirrep;i++){
    fscanf(file,"%s",irrep[i]);    
  }
  fscanf(file,"%d",&norb);
  no=(int*)calloc(norb,sizeof(int));
  irreps=(char**)calloc(norb,sizeof(char*));
  for(i=0;i<norb;i++){
    irreps[i]=(char*)calloc(5,sizeof(char));
  }
  for(io=0;io<norb;io++){
    fscanf(file,"%d%s",&no[io],irreps[io]);
  }
  fscanf(file,"%d",&nb);
  ini=(int*)calloc(4*nb,sizeof(int));
  for(i=0;i<nb;i++){
    for(j=0;j<4;j++){
      fscanf(file,"%d",&ini[j+4*i]);
    }
  }
  fclose(file);
  hd=(hamilton*)malloc(nirrep*sizeof(hamilton));
  od=(orbitals*)malloc(norb*nirrep*sizeof(orbitals));
  for(k=0;k<nirrep;k++){
    if(id==root){
      lmerb=0;
      _genlm(symm,irrep[k],lmax);
      strcpy(fileel,"angular_");
      strcat(fileel,symm);
      strcat(fileel,"_");
      strcat(fileel,irrep[k]);
      strcat(fileel,".inp");
      file=fopen(fileel,"r");
      for(i=getc(file);i!=EOF;i=getc(file)){
	if(i=='\n'){
	  lmerb=lmerb+1;
	}
      }
      fclose(file);
    }    
    MPI_Bcast(&lmerb,1,MPI_INT,root,PETSC_COMM_WORLD);
    e_nl=(cmplx*)calloc(nbas*lmerb,sizeof(cmplx));    
    c_nl=(cmplx*)calloc(nbas*nbas*lmerb*lmerb,sizeof(cmplx));
    le=(int*)calloc(lmerb,sizeof(int));
    me=(int*)calloc(lmerb,sizeof(int));
    if(id==root){
      file=fopen(fileel,"r");
      for(i=0;i<lmerb;i++){
	fscanf(file,"%d%d%s",&le[i],&me[i],str);
      }
      fclose(file);
    }
    if(id==root){
      printf("Candelones con tostones %s %d\n",str,lmerb);
      fflush(stdout);
    }
    MPI_Bcast(le,lmerb,MPI_INT,root,PETSC_COMM_WORLD);
    MPI_Bcast(me,lmerb,MPI_INT,root,PETSC_COMM_WORLD);
    MPI_Bcast(&str,strlen(str)+1,MPI_CHAR,root,PETSC_COMM_WORLD);
    for(io=0;io<norb;io++){
      if(id==root){
	CORR._GETORB(fichero,symm,irreps[io],no[io],strs,&lmorb,&lo,&mo,&cilm);
      }
      MPI_Bcast(&lmorb,1,MPI_INT,root,PETSC_COMM_WORLD);
      MPI_Bcast(&strs,strlen(strs)+1,MPI_CHAR,root,PETSC_COMM_WORLD);
      if(id!=root){
	lo=(int*)calloc(lmorb,sizeof(int));
	mo=(int*)calloc(lmorb,sizeof(int));
	cilm=(cmplx*)calloc(nbas*lmorb,sizeof(cmplx));
      }
      MPI_Bcast(lo,lmorb,MPI_INT,root,PETSC_COMM_WORLD);
      MPI_Bcast(mo,lmorb,MPI_INT,root,PETSC_COMM_WORLD);
      MPI_Bcast(cilm,nbas*lmorb,MPI_DOUBLE_COMPLEX,root,PETSC_COMM_WORLD);
      
      //J
      if(id==root){
	start_time=MPI_Wtime();
      }
      _GETC3J(id,lmax,str,strs,lmerb,lmorb,le,me,lo,mo,&CKj,&CKi,&CJj,&CJi,1);
      if(id==root){
	end_time=MPI_Wtime();
	printf("tiempo CK % .8le\n",end_time-start_time);
      }
      
      J=(cmplx*)calloc(nbas*nbas*lmerb*lmerb,sizeof(cmplx));
      CORR._JO(id,numprocs,lmerb,le,lmorb,lo,cilm,CJj,CJi,J);
      free(CJi);
      free(CJj);
      
      K=(cmplx*)calloc(nbas*nbas*lmerb*lmerb,sizeof(cmplx));
      CORR._KO(id,numprocs,lmerb,le,lmorb,lo,cilm,CKj,CKi,K);
      free(CKi);
      free(CKj);
      
      _INITORB(io+norb*k,nbas,lmerb,strs,lmorb,lo,mo,cilm,J,K,&(od[0]));    
      CORR._JMK(id,numprocs,io+norb*k,&(od[0]));     

      free(K);
      free(J);
      free(lo);
      free(mo);
      free(cilm);      
    }

    if(id==root){
      start_time=MPI_Wtime();
    }
    _GETC3J(id,lamax,str,str,lmerb,lmerb,le,me,le,me,&CKj,&CKi,&CJj,&CJi,2);
    if(id==root){
      end_time=MPI_Wtime();
      printf("tiempo CK % .8le\n",end_time-start_time);
    }    
    NE=(cmplx*)calloc(nbas*nbas*lmerb*lmerb,sizeof(cmplx));
    CORR._NE(id,numprocs,lmerb,le,CKj,CKi,NE);
    free(CKi);    
    free(CKj);
    free(CJi);
    free(CJj);

    CORR._HH1(id,numprocs,k,norb,str,lmerb,le,me,&(od[0]),NE,c_nl,e_nl);
    _INITH(k,nbas,str,lmerb,le,me,c_nl,e_nl,&(hd[0]));    
    
    free(NE);
    free(le);
    free(me);
    free(e_nl);
    free(c_nl);
  }
  free(no);
  free(od);
  free(irrep);
  free(irreps);
  
  ///two-electron-H
  n=nb*(nb+1)/2;
  nmkl=(int*)calloc(4*n,sizeof(int));
  _mblock(nf,nb,nmkl);
  if(id==root){
    file=fopen("mblock.dat","w");
    fprintf(file,"%d\n",n);
    k=0;
    for(i=0;i<nb;i++){
      for(j=i;j<nb;j++){
	fprintf(file,"%5d %5d %5d %5d %5d\n",k+1,nmkl[0+4*k],nmkl[1+4*k],nmkl[2+4*k],nmkl[3+4*k]);
	k=k+1;
      }
    }
    fclose(file);
  }

  nvec=nb*nf*nf;
  nloc=(int)(floor(nvec/numprocs));
  if(id<nvec-numprocs*nloc){
    nloc=nloc+1;
  }
  d_nz=nloc;     
  o_nz=nvec-nloc;
  MatCreateAIJ(PETSC_COMM_WORLD,nloc,nloc,nvec,nvec,d_nz,PETSC_NULL,o_nz,PETSC_NULL,&AA);  
  k=0;
  for(il=0;il<nb;il++){
    for(jl=il;jl<nb;jl++){      
      if(id==root){
	AT=(double*)calloc(nf*nf*nf*nf,sizeof(double));
      }
      if(il==jl){
	type=0;
      }
      else{
	type=1;
      }
      int ijkl[4]={ini[0+4*il]-1,ini[1+4*il]-1,ini[0+4*jl]-1,ini[1+4*jl]-1};
      int ijkl0[4]={ini[2+4*il],ini[3+4*il],ini[2+4*jl],ini[3+4*jl]};
      CORR._HH2(id,numprocs,ijkl,ijkl0,nf,nf,&(hd[0]),AT,symm,type);
      
      if(id==root){
	for(i1=0;i1<nf*nf;i1++){
	  for(j1=0;j1<nf*nf;j1++){			
	    i=i1+nmkl[0+4*k]-1;
	    j=j1+nmkl[2+4*k]-1;
	    vAA=AT[j1+nf*nf*i1]+0.0*PETSC_i;
	    MatSetValue(AA,i,j,vAA,INSERT_VALUES);
	    if(il!=jl){
	      MatSetValue(AA,j,i,vAA,INSERT_VALUES);
	    }
	  }
	}
	free(AT);
      }
      k++;
    }
  }
  MatAssemblyBegin(AA,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(AA,MAT_FINAL_ASSEMBLY);
    
  EPSCreate(PETSC_COMM_WORLD,&eps);
  EPSSetOperators(eps,AA,NULL);
  EPSSetProblemType(eps,EPS_HEP);
  EPSSetFromOptions(eps);
  EPSSolve(eps);
  EPSGetIterationNumber(eps,&its);
  EPSGetConverged(eps,&nconv);
  EPSComputeError(eps,0,EPS_ERROR_RELATIVE,&error);
  if(id==root){
    printf("ERROR=% .8le ITS=%d\n",error,its);
  }

  VecCreateMPI(PETSC_COMM_WORLD,nloc,nvec,&x);
  VecDuplicate(x,&y);
  VecDuplicateVecs(x,nconv,&Q);  
  for(i=0;i<nconv;i++){
    //    EPSGetEigenpair(eps,i,&kr,&ki,NULL,NULL);
    EPSGetEigenpair(eps,i,&kr,&ki,x,y);
    VecCopy(x,Q[i]);
    VecNorm(Q[i],NORM_2,&norm);
    cval=1.0/norm;
    VecScale(Q[i],cval);

    VecNorm(Q[i],NORM_2,&norm);
    if(id==root){
      printf("%d % .8le % .8le % .8le\n",i,norm,PetscRealPart(kr),PetscImaginaryPart(ki));
    }
  }
  
  //  VecCreateMPI(PETSC_COMM_WORLD,nloc,nvec,&x);
  //  VecDuplicate(x,&y);
  VecZeroEntries(x);
  VecZeroEntries(y);
  EPSGetEigenpair(eps,0,&kr,&ki,x,y);
  VecNorm(x,NORM_2,&norm);
  cval=1.0/norm;
  VecScale(x,cval);

  VecNorm(x,NORM_2,&norm);
  if(id==root){
    printf("% .8le\n",norm);
  }
  PetscViewerBinaryOpen(PETSC_COMM_WORLD,solution,FILE_MODE_WRITE,&view);
  VecView(x,view);
  PetscViewerDestroy(&view);  
  EPSDestroy(&eps);

  PetscViewerBinaryOpen(PETSC_COMM_WORLD,"AA.m",FILE_MODE_WRITE,&view);
  MatView(AA,view);
  PetscViewerDestroy(&view);
  MatDestroy(&AA);  
  
  MatCreateAIJ(PETSC_COMM_WORLD,nloc,nloc,nvec,nvec,d_nz,PETSC_NULL,o_nz,PETSC_NULL,&BB);  
  k=0;
  for(il=0;il<nb;il++){
    for(jl=il;jl<nb;jl++){      
      if(id==root){
	AT=(double*)calloc(nf*nf*nf*nf,sizeof(double));
      }
      int ijkl[4]={ini[0+4*il]-1,ini[1+4*il]-1,ini[0+4*jl]-1,ini[1+4*jl]-1};
      int ijkl0[4]={ini[2+4*il],ini[3+4*il],ini[2+4*jl],ini[3+4*jl]};
      CORR._VVz(id,numprocs,ijkl,ijkl0,nf,nf,&(hd[0]),AT,&Vnuc);
      
      if(id==root){
	for(i1=0;i1<nf*nf;i1++){
	  for(j1=0;j1<nf*nf;j1++){			
	    i=i1+nmkl[0+4*k]-1;
	    j=j1+nmkl[2+4*k]-1;
	    if(il==jl && i1==j1 ){
	      vAA=AT[j1+nf*nf*i1]+Vnuc+0.0*PETSC_i;
	    }
	    else{
	      vAA=AT[j1+nf*nf*i1]+0.0*PETSC_i;
	    }
	    MatSetValue(BB,i,j,vAA,INSERT_VALUES);
	    if(il!=jl){
	      MatSetValue(BB,j,i,vAA,INSERT_VALUES);
	    }
	  }
	}
	free(AT);
      }
      k++;
    }
  }
  MatAssemblyBegin(BB,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(BB,MAT_FINAL_ASSEMBLY);

  for(i=0;i<nconv;i++){
    VecZeroEntries(y);
    MatMult(BB,x,y);    
    VecDot(Q[i],y,&ki);
    if(id==root){
      printf("%d % .8le\n",i,(double)PetscRealPart(ki));
    }    
  }

  VecDestroy(&x);
  VecDestroy(&y);
  VecDestroyVecs(nconv,&Q);

  PetscViewerBinaryOpen(PETSC_COMM_WORLD,"BB.m",FILE_MODE_WRITE,&view);
  MatView(BB,view);
  PetscViewerDestroy(&view);
  MatDestroy(&BB);
  
  free(ini);  
  free(nmkl);
  free(hd);
 
  
  SlepcFinalize();
  return 0;
}
