#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <gsl/gsl_vector.h>
#include <petsc.h>
#include <mpi.h>
#include "mkl.h"
#include "utils.h"
#include "class_grid.h"
#include "orb.h"
#include "sphere_lebedev_rule.h"
#include "sasha.h"
#include "bieleclist.h"

GRID::GRID(int nbas,int lamax,cmplx *Xi,cmplx *Wi,cmplx *TXX,cmplx *TIXX){ 
  int i;
  nbreak=nbas;
  nmlps=lamax;
  Xz=(cmplx*)calloc(nbreak,sizeof(cmplx));
  Wz=(cmplx*)calloc(nbreak,sizeof(cmplx));
  for(i=0;i<nbreak;i++){
    Xz[i]=Xi[i];
    Wz[i]=Wi[i];
  }
  KE=(cmplx*)calloc(nbreak*nbreak*nmlps,sizeof(cmplx));
  KEI=(cmplx*)calloc(nbreak*nbreak*nmlps,sizeof(cmplx));
  for(i=0;i<nbreak*nbreak*nmlps;i++){
    KE[i]=TXX[i];
    KEI[i]=TIXX[i];
  }
}

void GRID::_HH1(int id,int numprocs,int kk,int norb,char *str,int lmerb,int *le,int *me,orbitals *p,cmplx *NE,cmplx *c_nl,cmplx *e_nl){
  FILE *file;
  int i,j,i1,j1,li,lf,local_n,lmorb,io;
  int lwork,info=0,n;
  double *AT,*work,*w,start_time,end_time;
  cmplx beta,*A,z1,*B,*eo,E0;
  char sto[5],fichero[25]="angular_";
  char jobz='V',uplo='L';
  lmorb=lmerb;
  A=(cmplx*)calloc(nbreak*nbreak*lmorb*lmorb,sizeof(cmplx));
  local_n=(int)(floor(nbreak*lmorb/numprocs));
  if(id<nbreak*lmorb-numprocs*local_n){
    local_n=local_n+1;
  }
  for(i=0;i<local_n;i++){
    li=floor((i+id*local_n)/nbreak); 
    i1=i+id*local_n-nbreak*li;
    for(j=0;j<nbreak*lmorb;j++){
      lf=floor(j/nbreak);
      j1=j-nbreak*lf;
      if(li==lf){
      	A[j+lmorb*nbreak*(i+id*local_n)]=KE[j1+nbreak*i1+nbreak*nbreak*le[li]];
      }
      for(io=0;io<norb;io++){
      	A[j+lmorb*nbreak*(i+id*local_n)]=A[j+lmorb*nbreak*(i+id*local_n)]+2.0*p[io+norb*kk].J[j+lmorb*nbreak*(i+id*local_n)]-p[io+norb*kk].K[j+lmorb*nbreak*(i+id*local_n)];
      }
      A[j+lmorb*nbreak*(i+id*local_n)]=A[j+lmorb*nbreak*(i+id*local_n)]+NE[j+lmorb*nbreak*(i+id*local_n)];
    }
  }
  MPI_Allgather(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,A,local_n*nbreak*lmorb,MPI_DOUBLE_COMPLEX,PETSC_COMM_WORLD);

  eo=(cmplx*)calloc(norb,sizeof(cmplx));
  for(io=0;io<norb;io++){
    B=(cmplx*)calloc(lmorb*nbreak,sizeof(cmplx));
    for(i=0;i<local_n;i++){
      z1=0.0;
      for(j=0;j<lmorb*nbreak;j++){
	z1=z1+A[j+lmorb*nbreak*(i+id*local_n)]*p[io+norb*kk].cilm[j];
      }
      B[i+id*local_n]=z1;
    }
    MPI_Allgather(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,B,local_n,MPI_DOUBLE_COMPLEX,PETSC_COMM_WORLD);
    if(id==0){
      z1=0.0;
      for(i=0;i<lmorb*nbreak;i++){
	z1=z1+B[i]*conj(p[io+norb*kk].cilm[i]);
      }
      eo[io]=z1;
      printf("Energy= % .8le % .8le\n",creal(z1),cimag(z1));
    }
    free(B);
  }

  //core energy
  MPI_Bcast(eo,norb,MPI_DOUBLE_COMPLEX,0,PETSC_COMM_WORLD);
  E0=0.0;
  for(io=0;io<norb;io++){
    B=(cmplx*)calloc(lmorb*nbreak,sizeof(cmplx));
    for(i=0;i<local_n;i++){
      li=floor((i+id*local_n)/nbreak);
      i1=i+id*local_n-nbreak*li;
      z1=0.0;
      for(j=0;j<lmorb*nbreak;j++){
	lf=floor(j/nbreak);
	j1=j-nbreak*lf;
	if(li==lf){
	  z1=z1+(KE[j1+nbreak*i1+nbreak*nbreak*le[li]]+NE[j+lmorb*nbreak*(i+id*local_n)])*p[io+norb*kk].cilm[j];
	}
	else{
	  z1=z1+NE[j+lmorb*nbreak*(i+id*local_n)]*p[io+norb*kk].cilm[j];
	}
      }
      B[i+id*local_n]=z1;
    }
    MPI_Allgather(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,B,local_n,MPI_DOUBLE_COMPLEX,PETSC_COMM_WORLD);
    if(id==0){
      z1=0.0;
      for(i=0;i<lmorb*nbreak;i++){
	z1=z1+B[i]*conj(p[io+norb*kk].cilm[i]);
      }
      E0=z1+eo[io];
    }
    free(B);
  }
  free(eo);
  MPI_Bcast(&E0,1,MPI_DOUBLE_COMPLEX,0,PETSC_COMM_WORLD);
  if(id==0){
    printf("Core Energy=% .8le % .8le\n",creal(E0),cimag(E0));
  }
  //fin core energy
  
  //get the real orbitals version cmplx en otro fichero  
  if(id==0){
    start_time=MPI_Wtime();
    n=nbreak*lmorb;
    AT=(double*)calloc(n*n,sizeof(double));
    for(i=0;i<n*n;i++){
      AT[i]=creal(A[i]);
    }
    lwork=4*n;
    w=(double*)calloc(n,sizeof(double));
    work=(double*)calloc(lwork,sizeof(double));
    lwork=-1;
    dsyev_(&jobz,&uplo,&n,AT,&n,w,work,&lwork,&info);
    lwork=(int)(work[0]);
    free(work);
    work=(double*)calloc(lwork,sizeof(double));
    dsyev_(&jobz,&uplo,&n,AT,&n,w,work,&lwork,&info);        
    for(i=0;i<n;i++){
      e_nl[i]=w[i];
      printf("Eigenvalue %d % .8le\n",i,w[i]); 
      for(j=0;j<n;j++){
	c_nl[j+n*i]=AT[j+i*n];
      }
    }
    free(w);
    free(work);
    free(AT);
    // for(j=0;j<2;j++){
    //   strcpy(fichero,"orbital_di_");
    //   sprintf(sto,"%d",j);
    //   strcat(fichero,sto);
    //   strcat(fichero,"_");
    //   sprintf(sto,"%d",kk);
    //   strcat(fichero,sto);
    //   strcat(fichero,".inp");
    //   file=fopen(fichero,"w");
    //   for(i=0;i<nbreak;i++){    
    // 	beta=0.0;
    // 	for(li=0;li<lmorb;li++){
    // 	  beta=beta+c_nl[i+nbreak*li+n*j]*_Blm(str,le[li],me[li],0.0,0.0)/(cpow(Wz[i],0.5)*Xz[i]);
    // 	}
    // 	fprintf(file,"% .8le % .8le % .8le\n",creal(Xz[i]),creal(beta),cimag(beta));
    //   }
    //   fclose(file);
    end_time=MPI_Wtime();
    printf("tiempo Diag=% .8le\n",end_time-start_time);
    fflush(stdout);
  }
  free(A);
  
  MPI_Bcast(e_nl,nbreak*lmorb,MPI_DOUBLE_COMPLEX,0,PETSC_COMM_WORLD);
  MPI_Bcast(c_nl,nbreak*lmorb*nbreak*lmorb,MPI_DOUBLE_COMPLEX,0,PETSC_COMM_WORLD);
}

void GRID::_JMK(int id,int numprocs,int io,orbitals *p){
  int i,j,i1,j1,li,lf,k,local_n,lmorb;
  double max;
  cmplx z1,*A;
  lmorb=p[io].lmerb;
  if(id==0){
    gsl_vector *v1=gsl_vector_alloc(lmorb*nbreak*lmorb*nbreak);
    gsl_vector *v2=gsl_vector_alloc(lmorb*nbreak*lmorb*nbreak);
    for(i=0;i<lmorb*nbreak;i++){
      for(j=0;j<lmorb*nbreak;j++){
	k=j+lmorb*nbreak*i;
	gsl_vector_set(v1,k,fabs(creal(p[io].K[i+nbreak*lmorb*j])-creal(p[io].K[j+nbreak*lmorb*i])));
	gsl_vector_set(v2,k,fabs(creal(p[io].J[i+nbreak*lmorb*j])-creal(p[io].J[j+nbreak*lmorb*i])));
      }
    }
    max=gsl_vector_max(v2);
    printf("Max(Jij-Jji)=% .8le\n",max);
    max=gsl_vector_max(v1);
    printf("Max(Kij-Kji)=% .8le\n",max);
    gsl_vector_free(v1);
    gsl_vector_free(v2);
  }
}

cmplx GRID::_PSSN(int i,int j,int v){
  cmplx beta;
  beta=(2.0*v+1.0)*KEI[j+nbreak*i+nbreak*nbreak*v]/(Xz[j]*Xz[i]*cpow(Wz[i]*Wz[j],0.5));
  beta=beta+cpow(Xz[j]*Xz[i],v)/cpow(Xz[nbreak-1],2.0*v+1.0);
  return beta;
}

cmplx GRID::_LGRN(int i,cmplx x){
  int j;
  cmplx beta;
  beta=1.0;
  for(j=0;j<nbreak;j++){
    if(j!=i){
      beta=beta*(x-Xz[j])/(Xz[i]-Xz[j]);
    }
  }
  beta=beta/cpow(Wz[i],0.5);  
  return beta;
}

void GRID::_KO(int id,int numprocs,int lmerb,int *le,int lmorb,int *lo,cmplx *cilm,cmplx *C3j,cmplx *C3i,cmplx *K){
  int i,j,k,l,nsph,i1,i2,j1,j2,l1,l2,m1,m2,li,lf,local_n;
  int lmax,nchannels,*la,*ma;
  double start_time,end_time;
  cmplx beta,z1,z2,z3,z4,z5;
  lmax=(int)(floor(nmlps/2));
  nsph=(lmax-1)*(lmax+1)+1;
  la=(int*)calloc(nsph,sizeof(int));
  ma=(int*)calloc(nsph,sizeof(int));
  k=0;
  for(j=0;j<lmax;j++){
    for(i=0;i<=j;i++){
      if(i==0){
	la[k]=j;
	ma[k]=i;
	k++;
      }
      else{
	la[k]=j;
	ma[k]=i;
	k++;
	la[k]=j;
	ma[k]=-i;
	k++;
      }
    }
  }
  local_n=(int)(floor(nbreak*lmerb/numprocs));
  if(id<nbreak*lmerb-numprocs*local_n){
    local_n=local_n+1;
  }
  if(id==0){
    start_time=MPI_Wtime();
  }
  for(i=0;i<local_n;i++){
    li=floor((i+id*local_n)/nbreak); 
    i1=i+id*local_n-nbreak*li; 
    for(j=0;j<nbreak*lmerb;j++){ 
      lf=floor(j/nbreak);
      j1=j-nbreak*lf;     
      z1=0.0+0.0*I;
      for(l=0;l<nsph;l++){
	z2=4.0*M_PI*pow(-1,ma[l])*_PSSN(i1,j1,la[l])/(2.0*la[l]+1.0);
	for(l1=0;l1<lmorb;l1++){
	  if(la[l]>=abs(lo[l1]-le[li]) && la[l]<=lo[l1]+le[li]){ //esto es nuevo
	    z3=C3j[l+nsph*(l1+lmorb*li)];
	    for(l2=0;l2<lmorb;l2++){
	      if(la[l]>=abs(lo[l2]-le[lf]) && la[l]<=lo[l2]+le[lf]){ //esto es nuevo
		z4=C3i[l+nsph*(l2+lmorb*lf)];
		z5=cilm[i1+nbreak*l1]*cilm[j1+nbreak*l2];
		z1=z1+z2*z3*z4*z5;
	      }
	    }
	  }
	}
      }
      K[j+lmerb*nbreak*(i+id*local_n)]=z1;      
    }
  }
  MPI_Allgather(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,K,local_n*nbreak*lmerb,MPI_DOUBLE_COMPLEX,PETSC_COMM_WORLD);
  if(id==0){
    end_time=MPI_Wtime();
    printf("% .8le\n",end_time-start_time);
    fflush(stdout);
  }  
  free(la);
  free(ma);
}

void GRID::_JO(int id,int numprocs,int lmerb,int *le,int lmorb,int *lo,cmplx *cilm,cmplx *C3j,cmplx *C3i,cmplx *J){
  int i,j,k,l,nsph,i1,i2,j1,j2,l1,l2,m1,m2,li,lf,local_n;
  int lmax,nchannels,*la,*ma;
  double start_time,end_time;
  cmplx beta,z1,z2,z3,z4,z5;
  lmax=(int)(floor(nmlps/2));
  nsph=(lmax-1)*(lmax+1)+1;
  la=(int*)calloc(nsph,sizeof(int));
  ma=(int*)calloc(nsph,sizeof(int));
  k=0;
  for(j=0;j<lmax;j++){
    for(i=0;i<=j;i++){
      if(i==0){
	la[k]=j;
	ma[k]=i;
	k++;
      }
      else{
	la[k]=j;
	ma[k]=i;
	k++;
	la[k]=j;
	ma[k]=-i;
	k++;
      }
    }
  }
  
  local_n=(int)(floor(nbreak*lmerb/numprocs));
  if(id<nbreak*lmerb-numprocs*local_n){
    local_n=local_n+1;
  }
  if(id==0){
    start_time=MPI_Wtime();
  }
  for(i=0;i<local_n;i++){
    li=floor((i+id*local_n)/nbreak); 
    i1=i+id*local_n-nbreak*li; 
    for(j=0;j<nbreak*lmerb;j++){ 
      lf=floor(j/nbreak);
      j1=j-nbreak*lf;      
      z1=0.0+0.0*I;
      if(i1==j1){
	for(k=0;k<nbreak;k++){
	  for(l=0;l<nsph;l++){
	    if(la[l]>=abs(le[li]-le[lf]) && la[l]<=le[li]+le[lf]){ //esto es nuevo
	      z2=4.0*M_PI*pow(-1,ma[l])*_PSSN(i1,k,la[l])/(2.0*la[l]+1.0);
	      z3=C3j[l+nsph*(lf+lmerb*li)];
	      for(l1=0;l1<lmorb;l1++){
		for(l2=0;l2<lmorb;l2++){
		  if(la[l]>=abs(lo[l1]-lo[l2]) && la[l]<=lo[l1]+lo[l2]){ //esto es nuevo
		    z4=C3i[l+nsph*(l2+lmorb*l1)];
		    z5=cilm[k+nbreak*l1]*cilm[k+nbreak*l2];
		    z1=z1+z2*z3*z4*z5;
		  }
		}
	      }
	    }
      	  }
      	}
      }
      J[j+lmerb*nbreak*(i+id*local_n)]=z1;      
    }
  }
  MPI_Allgather(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,J,local_n*nbreak*lmerb,MPI_DOUBLE_COMPLEX,PETSC_COMM_WORLD);
  if(id==0){
    end_time=MPI_Wtime();
    printf("% .8le\n",end_time-start_time);
    fflush(stdout);
  }  
  free(la);
  free(ma);
}

void GRID::_NE(int id,int numprocs,int lmerb,int *le,cmplx *C3j,cmplx *C3i,cmplx *NE){
  FILE *file;
  int i,j,k,l,m,nsph,i1,i2,j1,j2,l1,l2,m1,m2,li,lf,local_n;
  int lmax,nchannels,*la,*ma,ncent,*Z;
  double start_time,end_time,**center_xyz,*R,theta,phi;
  cmplx beta,z1,z2,z3,z4,z5,***In,yi,*Ia;
  char str[50];
  const double AnsgAu=1.889726124565062;
  lmax=nmlps;
  nsph=(lmax-1)*(lmax+1)+1;
  la=(int*)calloc(nsph,sizeof(int));
  ma=(int*)calloc(nsph,sizeof(int));
  k=0;
  for(j=0;j<lmax;j++){
    for(i=0;i<=j;i++){
      if(i==0){
	la[k]=j;
	ma[k]=i;
	k++;
      }
      else{
	la[k]=j;
	ma[k]=i;
	k++;
	la[k]=j;
	ma[k]=-i;
	k++;
      }
    }
  }
  
  file=fopen("read.inp","r");
  fscanf(file,"%d",&ncent);
  fclose(file);
  Z=(int*)calloc(ncent,sizeof(double));
  center_xyz=(double**)calloc(ncent,sizeof(double));
  for(i=0;i<ncent;i++){
    center_xyz[i]=(double*)calloc(3,sizeof(double));
  }
  R=(double*)calloc(ncent,sizeof(double));
  file=fopen("center.inp","r");
  for(i=0;i<ncent;i++){
    fscanf(file,"%s%d%d",str,&k,&Z[i]);
    for(j=0;j<3;j++){
      fscanf(file,"%le",&center_xyz[i][j]);
      center_xyz[i][j]=center_xyz[i][j]*AnsgAu;
    }
    R[i]=sqrt(pow(center_xyz[i][0],2)+pow(center_xyz[i][1],2)+pow(center_xyz[i][2],2));
  }
  fclose(file);

  In=(cmplx***)calloc(ncent,sizeof(cmplx));
  for(i=0;i<ncent;i++){
    In[i]=(cmplx**)calloc(nbreak,sizeof(cmplx));
    for(j=0;j<nbreak;j++){
      In[i][j]=(cmplx*)calloc(lmax,sizeof(cmplx));
    }
  }
  Ia=(cmplx*)calloc(ncent*nbreak,sizeof(cmplx));

  // local_n=(int)(floor(nbreak/numprocs));
  // if(id<nbreak-numprocs*local_n){
  //   local_n=local_n+1;
  // }
  // for(m=0;m<local_n;m++){
  //   for(i=0;i<ncent;i++){
  //     Ia[i+ncent*(m+id*local_n)]=_LGRN(m+id*local_n,R[i])/R[i];
  //   }
  // }
  // MPI_Allgather(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,Ia,local_n*ncent,MPI_DOUBLE_COMPLEX,PETSC_COMM_WORLD);
  for(m=0;m<nbreak;m++){
    for(i=0;i<ncent;i++){
      Ia[i+ncent*m]=_LGRN(m,R[i])/R[i];
    }
  }
  
  for(i=0;i<ncent;i++){
    for(l=0;l<lmax;l++){
      for(k=0;k<nbreak;k++){
	yi=0.0;
	for(m=0;m<nbreak;m++){
	  z1=Ia[i+ncent*m];
	  z2=(2.0*l+1.0)*KEI[m+nbreak*k+nbreak*nbreak*l]/(Xz[k]*cpow(Wz[k],0.5));
	  yi=yi+z1*z2;
	}
	z3=cpow(R[i]*Xz[k],l)/cpow(Xz[nbreak-1],2.0*l+1.0);
	yi=yi+z3;
	In[i][k][l]=-Z[i]*yi;
      }
    }
  }
  free(Ia);
  
  local_n=(int)(floor(nbreak*lmerb/numprocs));
  if(id<nbreak*lmerb-numprocs*local_n){
    local_n=local_n+1;
  }
  if(id==0){
    start_time=MPI_Wtime();
  }
  for(i=0;i<local_n;i++){
    li=floor((i+id*local_n)/nbreak); 
    i1=i+id*local_n-nbreak*li;
    for(j=0;j<nbreak*lmerb;j++){ 
      lf=floor(j/nbreak);
      j1=j-nbreak*lf;
      z1=0.0+0.0*I;
      if(i1==j1){
	for(k=0;k<ncent;k++){
	  theta=acos(center_xyz[k][2]/R[k]);
	  phi=atan2(center_xyz[k][1],center_xyz[k][0]);
	  if(id==0 && i1==0 && j1==0){
	    printf("Center %d % .8le % .8le\n",k,theta,phi);
	  }
	  if(k==0){
	    if(li==lf){
	      z1=-Z[k]/Xz[i1];
	    }
	  }
	  else{
	    for(l=0;l<nsph;l++){
	      if(la[l]>=abs(le[li]-le[lf]) && la[l]<=le[li]+le[lf]){ //esto es nuevo
		z2=4.0*M_PI*In[k][i1][la[l]]*C3j[l+nsph*(lf+lmerb*li)]/(2.0*la[l]+1.0);
		z1=z1+z2*conj(_Ylm(la[l],ma[l],theta,phi));
	      }
	    }
	  }
	}
      }
      NE[j+lmerb*nbreak*(i+id*local_n)]=z1;      
    }
  }
  MPI_Allgather(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,NE,local_n*nbreak*lmerb,MPI_DOUBLE_COMPLEX,PETSC_COMM_WORLD);
  if(id==0){
    end_time=MPI_Wtime();
    printf("% .8le\n",end_time-start_time);
    fflush(stdout);
  }  
  free(la);
  free(ma);
  free(Z);
  free(In);
  free(R);
  free(center_xyz);
}
 

// void GRID::_NE(int id,int numprocs,int lmerb,int *le,cmplx *C3j,cmplx *C3i,cmplx *NE){
//   FILE *file;
//   int i,j,k,l,nsph,i1,i2,j1,j2,l1,l2,m1,m2,li,lf,local_n;
//   int lmax,nchannels,*la,*ma,ncent,*Z;
//   double start_time,end_time,**center_xyz,yi,xi,*R,theta,phi;
//   cmplx beta,z1,z2,z3,z4,z5,***In;
//   char str[50];
//   const double AnsgAu=1.889726124565062;
//   lmax=nmlps;
//   nsph=(lmax-1)*(lmax+1)+1;
//   la=(int*)calloc(nsph,sizeof(int));
//   ma=(int*)calloc(nsph,sizeof(int));
//   k=0;
//   for(j=0;j<lmax;j++){
//     for(i=0;i<=j;i++){
//       if(i==0){
// 	la[k]=j;
// 	ma[k]=i;
// 	k++;
//       }
//       else{
// 	la[k]=j;
// 	ma[k]=i;
// 	k++;
// 	la[k]=j;
// 	ma[k]=-i;
// 	k++;
//       }
//     }
//   }
  
//   file=fopen("read.inp","r");
//   fscanf(file,"%d",&ncent);
//   fclose(file);
//   Z=(int*)calloc(ncent,sizeof(double));
//   center_xyz=(double**)calloc(ncent,sizeof(double));
//   for(i=0;i<ncent;i++){
//     center_xyz[i]=(double*)calloc(3,sizeof(double));
//   }
//   file=fopen("center.inp","r");
//   for(i=0;i<ncent;i++){
//     fscanf(file,"%s%d%d",str,&k,&Z[i]);
//     for(j=0;j<3;j++){
//       fscanf(file,"%le",&center_xyz[i][j]);
//       center_xyz[i][j]=center_xyz[i][j]*AnsgAu;
//     }
//   }
//   fclose(file);

//   R=(double*)calloc(ncent,sizeof(double));
//   In=(cmplx***)calloc(ncent,sizeof(cmplx));
//   for(i=0;i<ncent;i++){
//     In[i]=(cmplx**)calloc(nbreak,sizeof(cmplx));
//     for(j=0;j<nbreak;j++){
//       In[i][j]=(cmplx*)calloc(lmax,sizeof(cmplx));
//     }
//   }

//   for(i=0;i<ncent;i++){
//     R[i]=sqrt(pow(center_xyz[i][0],2)+pow(center_xyz[i][1],2)+pow(center_xyz[i][2],2));
//     for(l=0;l<lmax;l++){
//       for(k=0;k<nbreak;k++){
// 	xi=creal(Xz[k]);   //change double to cmplx
// 	if(xi<R[i]){
// 	  yi=pow(xi/R[i],l)*(1.0/R[i]);
// 	}
// 	else{
// 	  yi=pow(R[i]/xi,l)*(1.0/xi);
// 	}
// 	In[i][k][l]=-Z[i]*yi;
//       }
//     }
//   }
  
//   local_n=(int)(floor(nbreak*lmerb/numprocs));
//   if(id<nbreak*lmerb-numprocs*local_n){
//     local_n=local_n+1;
//   }
//   if(id==0){
//     start_time=MPI_Wtime();
//   }
//   for(i=0;i<local_n;i++){
//     li=floor((i+id*local_n)/nbreak); 
//     i1=i+id*local_n-nbreak*li;
//     for(j=0;j<nbreak*lmerb;j++){ 
//       lf=floor(j/nbreak);
//       j1=j-nbreak*lf;
//       z1=0.0+0.0*I;
//       if(i1==j1){
// 	for(k=0;k<ncent;k++){
// 	  theta=acos(center_xyz[k][2]/R[k]);
// 	  phi=atan2(center_xyz[k][1],center_xyz[k][0]);
// 	  if(id==0 && i1==0 && j1==0){
// 	    printf("Center %d % .8le % .8le\n",k,theta,phi);
// 	  }
// 	  if(k==0){
// 	    if(li==lf){
// 	      z1=-Z[k]/Xz[i1];
// 	    }
// 	  }
// 	  else{
// 	    for(l=0;l<nsph;l++){
// 	      if(la[l]>=abs(le[li]-le[lf]) && la[l]<=le[li]+le[lf]){ //esto es nuevo
// 		z2=4.0*M_PI*In[k][i1][la[l]]*C3j[l+nsph*(lf+lmerb*li)]/(2.0*la[l]+1.0);
// 		z1=z1+z2*conj(_Ylm(la[l],ma[l],theta,phi));
// 	      }
// 	    }
// 	  }
// 	}
//       }
//       NE[j+lmerb*nbreak*(i+id*local_n)]=z1;      
//     }
//   }
//   MPI_Allgather(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,NE,local_n*nbreak*lmerb,MPI_DOUBLE_COMPLEX,PETSC_COMM_WORLD);
//   if(id==0){
//     end_time=MPI_Wtime();
//     printf("% .8le\n",end_time-start_time);
//     fflush(stdout);
//   }  
//   free(la);
//   free(ma);
//   free(Z);
//   free(In);
//   free(R);
//   free(center_xyz);
// }

void GRID::_GETORB(char *orbfile,char *symm,char *irrep,int lo,char *str1,int *lmorb,int **li,int **mi,cmplx **cilm){
  FILE *file;
  int i,j,k,l,ngrid,v,nchannels;
  double *orb,*x,*y,*z,*xi,*yi,*zi,*wi,theta,phi;
  cmplx beta;
  char str[10],sto[5],fichero[25]="angular_";
  const int n=2354; //el que mejor funciona es 194
  nchannels=0;
  l=(int)(floor(nmlps/2));
  _genlm(symm,irrep,l);
  strcat(fichero,symm);  
  strcat(fichero,"_");
  strcat(fichero,irrep);
  strcat(fichero,".inp");
  file=fopen(fichero,"r");
  for(i=getc(file);i!=EOF;i=getc(file)){
    if(i=='\n'){ 
      nchannels=nchannels+1;
    }
  }
  fclose(file);
  *lmorb=nchannels;
  printf("%d\n",nchannels);
  *li=(int*)calloc(nchannels,sizeof(int));
  *mi=(int*)calloc(nchannels,sizeof(int));
  file=fopen(fichero,"r");
  for(i=0;i<nchannels;i++){
    fscanf(file,"%d%d%s",&(*li)[i],&(*mi)[i],str);
  }
  fclose(file);
  strcpy(str1,str);
  
  xi=(double*)calloc(n,sizeof(double));
  yi=(double*)calloc(n,sizeof(double));
  zi=(double*)calloc(n,sizeof(double));
  wi=(double*)calloc(n,sizeof(double));
  ld2354(xi,yi,zi,wi);
  ngrid=n*nbreak;
  orb=(double*)calloc(ngrid,sizeof(double));
  x=(double*)calloc(ngrid,sizeof(double));
  y=(double*)calloc(ngrid,sizeof(double));
  z=(double*)calloc(ngrid,sizeof(double));
  for(i=0;i<nbreak;i++){
    for(j=0;j<n;j++){
      x[j+n*i]=xi[j]*creal(Xz[i]);
      y[j+n*i]=yi[j]*creal(Xz[i]);
      z[j+n*i]=zi[j]*creal(Xz[i]);
    }
  }
  _getorb(orbfile,lo,ngrid,x,y,z,orb);

  *cilm=(cmplx*)calloc(nbreak*nchannels,sizeof(cmplx));
  for(l=0;l<nchannels;l++){
    for(i=0;i<nbreak;i++){
      beta=0.0+0.0*I;
      for(j=0;j<n;j++){
        theta=acos(zi[j]);
        phi=atan2(yi[j],xi[j]);
	if(strcmp(symm,"ylm")==0){
	  beta=beta+cpow(Wz[i],0.5)*wi[j]*orb[j+n*i]*conj(_Ylm((*li)[l],(*mi)[l],theta,phi))*Xz[i];
	}
	else{
	  beta=beta+cpow(Wz[i],0.5)*wi[j]*orb[j+n*i]*_Blm(str,(*li)[l],(*mi)[l],theta,phi)*Xz[i];
	}
      }
      (*cilm)[i+nbreak*l]=4.0*M_PI*beta;
    }
  }  
  beta=0.0+0.0*I;
  for(i=0;i<nbreak;i++){
    for(l=0;l<nchannels;l++){
      beta=beta+((*cilm)[i+nbreak*l])*conj((*cilm)[i+nbreak*l]);
    }
  }
  printf("%s %d%s Norm=% .8le % .8le\n",symm,lo+1,irrep,creal(beta),cimag(beta));
  fflush(stdout);
  for(l=0;l<nchannels;l++){
    for(i=0;i<nbreak;i++){
      (*cilm)[i+nbreak*l]=(*cilm)[i+nbreak*l]/cpow(beta,0.5);
    }
  }
  beta=0.0+0.0*I;
  for(i=0;i<nbreak;i++){
    for(l=0;l<nchannels;l++){
      beta=beta+((*cilm)[i+nbreak*l])*conj((*cilm)[i+nbreak*l]);
    }
  }
  printf("%s %d%s Norm=% .8le % .8le\n",symm,lo+1,irrep,creal(beta),cimag(beta));
  fflush(stdout);

  strcpy(fichero,"orbital_");
  sprintf(sto,"%d",lo);
  strcat(fichero,sto);
  strcat(fichero,".inp");
  file=fopen(fichero,"w");
  for(i=0;i<nbreak;i++){    
    beta=0.0;
    for(l=0;l<nchannels;l++){
      beta=beta+(*cilm)[i+nbreak*l]*_Blm(str,(*li)[l],(*mi)[l],0.0,0.0)/(cpow(Wz[i],0.5)*Xz[i]);
    }
    fprintf(file,"% .8le % .8le % .8le\n",creal(Xz[i]),creal(beta),cimag(beta));
  }
  fclose(file);
  
  free(xi);
  free(yi);
  free(zi);
  free(wi);
  free(x);
  free(y);
  free(z);
  free(orb);
}

void GRID::_HH2(int id,int numprocs,int nmkl[],int nmkl0[],int nf,int mf,hamilton *p,double *AT,char *symm,int type){
  FILE *file;
  int i,j,m1,n1,m2,n2,local_n,lmorb,io;
  int l,l1,l2,lmax,nsph,*la,*ma,k,k1,k2,lmerb1,lmerb2,nbas,*le1,*me1,*le2,*me2,n;
  int li1,li2,lf1,lf2,i1,i2,j1,j2,lmerb3,lmerb4,*le3,*me3,*le4,*me4,k3,k4;
  cmplx *CKj,*CKi,*CJj,*CJi,*B,*B1,*B12,z1,z2,z3,z4,z5;
  cmplx *CKj2,*CKi2,*CJj2,*CJi2;
  double start_time,end_time;
  char str[5];  
  int nbie,*ijkl,lnri;
  int i0,j0,k0,l0,n01,n03,m02,m04;
  
  i0=nmkl[0];
  j0=nmkl[1];
  k0=nmkl[2];
  l0=nmkl[3];
  n01=nmkl0[0];
  m02=nmkl0[1];
  n03=nmkl0[2];
  m04=nmkl0[3];
  if(id==0){
    printf("%d %d %d %d\n",i0,j0,k0,l0);
    printf("%d %d %d %d\n",n01,m02,n03,m04);
    fflush(stdout);
  }

  io=_dpr(symm,i0,j0);
  io=_dpr(symm,io,k0);
  io=_dpr(symm,io,l0);
  if(io==0){  
    lmerb1=p[i0].lmerb;
    le1=(int*)calloc(lmerb1,sizeof(int));
    me1=(int*)calloc(lmerb1,sizeof(int));
    for(l=0;l<lmerb1;l++){
      le1[l]=p[i0].le[l];
      me1[l]=p[i0].me[l];
    }
    lmerb2=p[j0].lmerb;
    le2=(int*)calloc(lmerb2,sizeof(int));
    me2=(int*)calloc(lmerb2,sizeof(int));
    for(l=0;l<lmerb2;l++){
      le2[l]=p[j0].le[l];
      me2[l]=p[j0].me[l];
    }
    lmerb3=p[k0].lmerb;
    le3=(int*)calloc(lmerb3,sizeof(int));
    me3=(int*)calloc(lmerb3,sizeof(int));
    for(l=0;l<lmerb3;l++){
      le3[l]=p[k0].le[l];
      me3[l]=p[k0].me[l];
    }  
    lmerb4=p[l0].lmerb;
    le4=(int*)calloc(lmerb4,sizeof(int));
    me4=(int*)calloc(lmerb4,sizeof(int));
    for(l=0;l<lmerb4;l++){
      le4[l]=p[l0].le[l];
      me4[l]=p[l0].me[l];
    }
    
    //bielectronicas
    //    lmax=(int)(floor(nmlps/2));
    lmax=nmlps;
    nsph=(lmax-1)*(lmax+1)+1;
    la=(int*)calloc(nsph,sizeof(int));
    ma=(int*)calloc(nsph,sizeof(int));
    k=0;
    for(j=0;j<lmax;j++){
      for(i=0;i<=j;i++){
	if(i==0){
	  la[k]=j;
	  ma[k]=i;
	  k++;
	}
	else{
	  la[k]=j;
	  ma[k]=i;
	  k++;
	  la[k]=j;
	  ma[k]=-i;
	  k++;
	}
      }
    }
    if(id==0){
      start_time=MPI_Wtime();
    }
    _GETC3J(id,lmax,p[i0].str,p[k0].str,lmerb1,lmerb3,le1,me1,le3,me3,&CKj,&CKi,&CJj,&CJi,1);
    free(CJj);
    free(CJi);
    _GETC3J(id,lmax,p[j0].str,p[l0].str,lmerb2,lmerb4,le2,me2,le4,me4,&CKj2,&CKi2,&CJj2,&CJi2,1);
    free(CJj2);
    free(CJi2);    
    if(id==0){
      end_time=MPI_Wtime();
      printf("tiempo CK % .8le\n",end_time-start_time);
      fflush(stdout);
    }
    
    if(id==0){
      start_time=MPI_Wtime();
    }  
    n=nbreak*lmerb2*lmerb4;
    local_n=(int)(floor(nbreak*lmerb1/numprocs));
    if(id<nbreak*lmerb1-numprocs*local_n){
      local_n=local_n+1;
    }
    B=(cmplx*)calloc(local_n*lmerb1*lmerb3*n,sizeof(cmplx));
    for(i2=0;i2<local_n;i2++){	  
      li1=floor((i2+id*local_n)/nbreak);
      i1=i2+id*local_n-nbreak*li1;	
      for(lf1=0;lf1<lmerb3;lf1++){
	i=lf1+lmerb3*i2;
	for(j1=0;j1<nbreak;j1++){
	  for(li2=0;li2<lmerb2;li2++){
	    for(lf2=0;lf2<lmerb4;lf2++){
	      j=lf2+lmerb4*li2+lmerb2*lmerb4*j1;
	      z1=0.0+0.0*I;
	      for(l=0;l<nsph;l++){
		if(la[l]>=abs(le1[li1]-le3[lf1]) && la[l]<=le1[li1]+le3[lf1] && la[l]>=abs(le2[li2]-le4[lf2]) && la[l]<=le2[li2]+le4[lf2]){ 
		  z2=4.0*M_PI*pow(-1,ma[l])*_PSSN(i1,j1,la[l])/(2.0*la[l]+1.0);
		  z3=CKi[l+nsph*(lf1+lmerb3*li1)];
		  z4=CKj2[l+nsph*(lf2+lmerb4*li2)];
		  z1=z1+z2*z3*z4;
		}
	      }
	      B[j+n*i]=z1;
	    }
	  }
	}
      }
    }
    
    if(id==0){
      end_time=MPI_Wtime();
      printf("tiempo electron-electron=% .8le\n",end_time-start_time);
      fflush(stdout);
    }
    free(la);
    free(ma);
    free(le1);
    free(me1);
    free(le2);
    free(me2);
    free(le3);
    free(me3);
    free(le4);
    free(me4);  
    free(CKj);
    free(CKi);
    free(CKj2);
    free(CKi2);

    // nbie=nf*nf*nf*nf;
    // ijkl=(int*)calloc(4*nbie,sizeof(int));
    // _bieleclist_2(nf,ijkl);
    if(i0==j0 && i0==k0 && i0==l0){
      nbie=nf*(nf+1)*(nf*nf+nf+2)/8;
      ijkl=(int*)calloc(4*nbie,sizeof(int));
      _bieleclist_8(nf,ijkl);
    }
    else if(i0==k0 && j0==l0 && i0!=j0 && k0!=l0){
      nbie=nf*nf*(nf+1)*(nf+1)/4;
      ijkl=(int*)calloc(4*nbie,sizeof(int));
      _bieleclist_4(nf,ijkl);
    }
    else{
      nbie=nf*nf*(nf*nf+1)/2;
      ijkl=(int*)calloc(4*nbie,sizeof(int));
      _bieleclist_2(nf,ijkl);
    }
    B12=(cmplx*)calloc(nbie,sizeof(cmplx));
    if(id==0){
      start_time=MPI_Wtime();
    }
    k1=nbreak*lmerb1;
    k2=nbreak*lmerb2;
    k3=nbreak*lmerb3;
    k4=nbreak*lmerb4;
    B1=(cmplx*)calloc(local_n*lmerb1*lmerb3,sizeof(cmplx));
    for(lnri=0;lnri<nbie;lnri++){
      n1=ijkl[0+4*lnri];
      n2=ijkl[1+4*lnri];
      m1=ijkl[2+4*lnri];
      m2=ijkl[3+4*lnri];
      
      z2=0.0+0.0*I;
      for(i2=0;i2<local_n;i2++){
	li1=floor((i2+id*local_n)/nbreak);
	i1=i2+id*local_n-nbreak*li1;	  
	for(lf1=0;lf1<lmerb3;lf1++){
	  i=lf1+lmerb1*i2;
	  z1=0.0+0.0*I;
	  for(j1=0;j1<nbreak;j1++){
	    for(li2=0;li2<lmerb2;li2++){
	      for(lf2=0;lf2<lmerb4;lf2++){
		j=lf2+lmerb4*li2+lmerb2*lmerb4*j1;
		z1=z1+B[j+n*i]*p[j0].c_nl[j1+nbreak*li2+k2*(m1+m02)]*p[l0].c_nl[j1+nbreak*lf2+k4*(m2+m04)];
	      }
	    }
	  }
	  B1[i]=z1;
	}
      }
      
      z1=0.0+0.0*I;
      for(i2=0;i2<local_n;i2++){
	li1=floor((i2+id*local_n)/nbreak);
	i1=i2+id*local_n-nbreak*li1;      
	for(lf1=0;lf1<lmerb3;lf1++){
	  i=lf1+lmerb1*i2;
	  z1=z1+B1[i]*p[i0].c_nl[i1+nbreak*li1+k1*(n1+n01)]*p[k0].c_nl[i1+nbreak*lf1+k3*(n2+n03)];
	}
      }
      
      MPI_Allreduce(&z1,&z2,1,MPI_DOUBLE_COMPLEX,MPI_SUM,PETSC_COMM_WORLD);    
      B12[lnri]=z2;    
    }
    free(B1);
    free(B);
    if(id==0){
      end_time=MPI_Wtime();
      printf("tiempo four index transf=% .8le\n",end_time-start_time);
      fflush(stdout);
    }    
    
    if(id==0){
      for(lnri=0;lnri<nbie;lnri++){
	n1=ijkl[0+4*lnri];
	n2=ijkl[1+4*lnri];
	m1=ijkl[2+4*lnri];
	m2=ijkl[3+4*lnri];
	
	if(i0==j0 && i0==k0 && i0==l0){
	  AT[m2+nf*n2+nf*nf*(m1+nf*n1)]=creal(B12[lnri]);
	  AT[m1+nf*n2+nf*nf*(m2+nf*n1)]=creal(B12[lnri]);
	  AT[m2+nf*n1+nf*nf*(m1+nf*n2)]=creal(B12[lnri]);
	  AT[m1+nf*n1+nf*nf*(m2+nf*n2)]=creal(B12[lnri]);
	  AT[n2+nf*m2+nf*nf*(n1+nf*m1)]=creal(B12[lnri]);
	  AT[n1+nf*m2+nf*nf*(n2+nf*m1)]=creal(B12[lnri]);
	  AT[n2+nf*m1+nf*nf*(n1+nf*m2)]=creal(B12[lnri]);
	  AT[n1+nf*m1+nf*nf*(n2+nf*m2)]=creal(B12[lnri]);
	}
	else if(i0==k0 && j0==l0 && i0!=j0 && k0!=l0){
	  AT[m2+nf*n2+nf*nf*(m1+nf*n1)]=creal(B12[lnri]);
	  AT[m1+nf*n2+nf*nf*(m2+nf*n1)]=creal(B12[lnri]);
	  AT[m2+nf*n1+nf*nf*(m1+nf*n2)]=creal(B12[lnri]);
	  AT[m1+nf*n1+nf*nf*(m2+nf*n2)]=creal(B12[lnri]);
	}	
	else{
	  AT[m2+nf*n2+nf*nf*(m1+nf*n1)]=creal(B12[lnri]);
	  AT[n2+nf*m2+nf*nf*(n1+nf*m1)]=creal(B12[lnri]);
	}
	
	if(cimag(B12[lnri])>1.0e-16){
	  printf("Mierda la has cagao % .8le\n",cimag(B12[lnri]));
	}
	
      }
    }
    free(B12);
    free(ijkl);
  }

  if(id==0){
    if(type==0){      
      for(n1=0;n1<nf;n1++){
      	for(m1=0;m1<nf;m1++){
	  AT[m1+nf*n1+nf*nf*(m1+nf*n1)]=AT[m1+nf*n1+nf*nf*(m1+nf*n1)]+creal(p[i0].e_nl[n1+n01]+p[j0].e_nl[m1+m02]);
      	}
      }
    }
  }  
}


void GRID::_VVz(int id,int numprocs,int nmkl[],int nmkl0[],int nf,int mf,hamilton *p,double *AT,double *Vnuc){
  FILE *file;
  int i,j,m1,n1,m2,n2,local_n,lmorb,io,ncent,*Z;
  int l,l1,l2,lmax,nsph,*la,*ma,k,k1,k2,lmerb1,lmerb2,nbas,*le1,*me1,*le2,*me2,n;
  int li1,li2,lf1,lf2,i1,i2,j1,j2,lmerb3,lmerb4,*le3,*me3,*le4,*me4,k3,k4;
  cmplx *CKj,*CKi,*CJj,*CJi,*B,*B1,*B12,z1,z2,z3,z4,z5;
  cmplx *CKj2,*CKi2,*CJj2,*CJi2;
  double start_time,end_time,**center_xyz,*R;
  char str[5];  
  int nbie,*ijkl,lnri;
  int i0,j0,k0,l0,n01,n03,m02,m04;
  const double AnsgAu=1.889726124565062;
  
  i0=nmkl[0];
  j0=nmkl[1];
  k0=nmkl[2];
  l0=nmkl[3];
  n01=nmkl0[0];
  m02=nmkl0[1];
  n03=nmkl0[2];
  m04=nmkl0[3];
  if(id==0){
    printf("%d %d %d %d\n",i0,j0,k0,l0);
    printf("%d %d %d %d\n",n01,m02,n03,m04);
    fflush(stdout);
  }

  lmerb1=p[i0].lmerb;
  le1=(int*)calloc(lmerb1,sizeof(int));
  me1=(int*)calloc(lmerb1,sizeof(int));
  for(l=0;l<lmerb1;l++){
    le1[l]=p[i0].le[l];
    me1[l]=p[i0].me[l];
  }
  lmerb2=p[j0].lmerb;
  le2=(int*)calloc(lmerb2,sizeof(int));
  me2=(int*)calloc(lmerb2,sizeof(int));
  for(l=0;l<lmerb2;l++){
    le2[l]=p[j0].le[l];
    me2[l]=p[j0].me[l];
  }
  lmerb3=p[k0].lmerb;
  le3=(int*)calloc(lmerb3,sizeof(int));
  me3=(int*)calloc(lmerb3,sizeof(int));
  for(l=0;l<lmerb3;l++){
    le3[l]=p[k0].le[l];
    me3[l]=p[k0].me[l];
  }  
  lmerb4=p[l0].lmerb;
  le4=(int*)calloc(lmerb4,sizeof(int));
  me4=(int*)calloc(lmerb4,sizeof(int));
  for(l=0;l<lmerb4;l++){
    le4[l]=p[l0].le[l];
    me4[l]=p[l0].me[l];
  }

  lmax=2;
  nsph=(lmax-1)*(lmax+1)+1;
  la=(int*)calloc(nsph,sizeof(int));
  ma=(int*)calloc(nsph,sizeof(int));
  k=0;
  for(j=0;j<lmax;j++){
    for(i=0;i<=j;i++){
      if(i==0){
	la[k]=j;
	ma[k]=i;
	k++;
      }
      else{
	la[k]=j;
	ma[k]=i;
	k++;
	la[k]=j;
	ma[k]=-i;
	k++;
      }
    }
  }
  
  if(id==0){
    start_time=MPI_Wtime();
  }
  _GETC3J(id,lmax,p[i0].str,p[k0].str,lmerb1,lmerb3,le1,me1,le3,me3,&CKj,&CKi,&CJj,&CJi,1);
  free(CJj);
  free(CJi);
  _GETC3J(id,lmax,p[j0].str,p[l0].str,lmerb2,lmerb4,le2,me2,le4,me4,&CKj2,&CKi2,&CJj2,&CJi2,1);
  free(CJj2);
  free(CJi2);    
  if(id==0){
    end_time=MPI_Wtime();
    printf("tiempo CK % .8le\n",end_time-start_time);
    fflush(stdout);
  }
  

  file=fopen("read.inp","r");
  fscanf(file,"%d",&ncent);
  fclose(file);
  Z=(int*)calloc(ncent,sizeof(double));
  center_xyz=(double**)calloc(ncent,sizeof(double));
  for(i=0;i<ncent;i++){
    center_xyz[i]=(double*)calloc(3,sizeof(double));
  }
  R=(double*)calloc(ncent,sizeof(double));
  file=fopen("center.inp","r");
  for(i=0;i<ncent;i++){
    fscanf(file,"%s%d%d",str,&k,&Z[i]);
    for(j=0;j<3;j++){
      fscanf(file,"%le",&center_xyz[i][j]);
      center_xyz[i][j]=center_xyz[i][j]*AnsgAu;
    }
    R[i]=sqrt(pow(center_xyz[i][0],2)+pow(center_xyz[i][1],2)+pow(center_xyz[i][2],2));
  }
  fclose(file);

  *Vnuc=0.0;
  for(i=0;i<ncent;i++){
    *Vnuc=*Vnuc+Z[i]*R[i];
  }
  free(Z);
  free(R);
  free(center_xyz);

  // nbie=nf*nf*nf*nf;
  // ijkl=(int*)calloc(4*nbie,sizeof(int));
  // _bieleclist_0(nf,ijkl);
  if(i0==j0 && i0==k0 && i0==l0){
    nbie=nf*(nf+1)*(nf*nf+nf+2)/8;
    ijkl=(int*)calloc(4*nbie,sizeof(int));
    _bieleclist_8(nf,ijkl);
  }
  else if(i0==k0 && j0==l0 && i0!=j0 && k0!=l0){
    nbie=nf*nf*(nf+1)*(nf+1)/4;
    ijkl=(int*)calloc(4*nbie,sizeof(int));
    _bieleclist_4(nf,ijkl);
  }
  else{
    nbie=nf*nf*(nf*nf+1)/2;
    ijkl=(int*)calloc(4*nbie,sizeof(int));
    _bieleclist_2(nf,ijkl);
  }
  B12=(cmplx*)calloc(nbie,sizeof(cmplx));
  if(id==0){
    start_time=MPI_Wtime();
  }
  k1=nbreak*lmerb1;
  k2=nbreak*lmerb2;
  k3=nbreak*lmerb3;
  k4=nbreak*lmerb4;

  for(lnri=0;lnri<nbie;lnri++){
    n1=ijkl[0+4*lnri];
    n2=ijkl[1+4*lnri];
    m1=ijkl[2+4*lnri];
    m2=ijkl[3+4*lnri];

    local_n=(int)(floor(nbreak*lmerb2/numprocs));
    if(id<nbreak*lmerb2-numprocs*local_n){
      local_n=local_n+1;
    }
    z3=0.0+0.0*I;
    z1=0.0+0.0*I;
    if(i0==k0 && n1+n01==n2+n03){

      for(i2=0;i2<local_n;i2++){
	li2=floor((i2+id*local_n)/nbreak);
	i1=i2+id*local_n-nbreak*li2;
	for(lf2=0;lf2<lmerb4;lf2++){	    
	  for(l=0;l<nsph;l++){
	    if((la[l]==1 && ma[l]==0) && la[l]==abs(le2[li2]-le4[lf2])){
	      z1=z1+Xz[i1]*p[j0].c_nl[i1+nbreak*li2+k2*(m1+m02)]*p[l0].c_nl[i1+nbreak*lf2+k4*(m2+m04)]*CKj2[l+nsph*(lf2+lmerb4*li2)];
	    }
	  }
	}
      }
    }
    MPI_Allreduce(&z1,&z3,1,MPI_DOUBLE_COMPLEX,MPI_SUM,PETSC_COMM_WORLD);

    local_n=(int)(floor(nbreak*lmerb1/numprocs));
    if(id<nbreak*lmerb1-numprocs*local_n){
      local_n=local_n+1;
    }
    z5=0.0+0.0*I;
    z1=0.0+0.0*I;
    if(j0==l0 && m1+m02==m2+m04){
      for(i2=0;i2<local_n;i2++){
	li1=floor((i2+id*local_n)/nbreak);
	i1=i2+id*local_n-nbreak*li1;
	for(lf1=0;lf1<lmerb1;lf1++){      
	  for(l=0;l<nsph;l++){
	    if((la[l]==1 && ma[l]==0) && la[l]==abs(le1[li1]-le3[lf1])){
	      z1=z1+Xz[i1]*p[i0].c_nl[i1+nbreak*li1+k1*(n1+n01)]*p[k0].c_nl[i1+nbreak*lf1+k3*(n2+n03)]*CKj[l+nsph*(lf1+lmerb3*li1)];
	    }
	  }
	}
      }
    }
    MPI_Allreduce(&z1,&z5,1,MPI_DOUBLE_COMPLEX,MPI_SUM,PETSC_COMM_WORLD);

    B12[lnri]=-sqrt(4.0*M_PI/3.0)*(z3+z5);
  }
  free(la);
  free(ma);
  free(le1);
  free(me1);
  free(le2);
  free(me2);
  free(le3);
  free(me3);
  free(le4);
  free(me4);  
  free(CKj);
  free(CKi);
  free(CKj2);
  free(CKi2);
  
  if(id==0){
    end_time=MPI_Wtime();
    printf("tiempo four two-index dipole transf=% .8le\n",end_time-start_time);
    fflush(stdout);
  }    

  
  if(id==0){    
    for(lnri=0;lnri<nbie;lnri++){
      n1=ijkl[0+4*lnri];
      n2=ijkl[1+4*lnri];
      m1=ijkl[2+4*lnri];
      m2=ijkl[3+4*lnri];
      
      if(i0==j0 && i0==k0 && i0==l0){
        AT[m2+nf*n2+nf*nf*(m1+nf*n1)]=creal(B12[lnri]);
        AT[m1+nf*n2+nf*nf*(m2+nf*n1)]=creal(B12[lnri]);
        AT[m2+nf*n1+nf*nf*(m1+nf*n2)]=creal(B12[lnri]);
        AT[m1+nf*n1+nf*nf*(m2+nf*n2)]=creal(B12[lnri]);
        AT[n2+nf*m2+nf*nf*(n1+nf*m1)]=creal(B12[lnri]);
        AT[n1+nf*m2+nf*nf*(n2+nf*m1)]=creal(B12[lnri]);
        AT[n2+nf*m1+nf*nf*(n1+nf*m2)]=creal(B12[lnri]);
        AT[n1+nf*m1+nf*nf*(n2+nf*m2)]=creal(B12[lnri]);
      }
      else if(i0==k0 && j0==l0 && i0!=j0 && k0!=l0){
        AT[m2+nf*n2+nf*nf*(m1+nf*n1)]=creal(B12[lnri]);
        AT[m1+nf*n2+nf*nf*(m2+nf*n1)]=creal(B12[lnri]);
        AT[m2+nf*n1+nf*nf*(m1+nf*n2)]=creal(B12[lnri]);
        AT[m1+nf*n1+nf*nf*(m2+nf*n2)]=creal(B12[lnri]);
      }	
      else{      
	AT[m2+nf*n2+nf*nf*(m1+nf*n1)]=creal(B12[lnri]);
      	AT[n2+nf*m2+nf*nf*(n1+nf*m1)]=creal(B12[lnri]);
      }
      
      if(cimag(B12[lnri])>1.0e-16){
	printf("Mierda la has cagao % .8le % .8le\n",creal(B12[lnri]),cimag(B12[lnri]));
      }      
    }
  }
  free(B12);
  free(ijkl);
  
}


