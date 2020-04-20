#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_vector.h>
#include <mpi.h>
#include "utils.h"
#include "class_grid.h"
#include "orb.h"
#include "sphere_lebedev_rule.h"
#include "sasha.h"

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

// void GRID::_HH1(double Zcharge,int id,int numprocs,int norb,orbitals *p,cmplx *NE,cmplx *H1){
//   int i,j,i1,j1,li,lf,local_n,lmorb,io;
//   double z1,*A,*B;
//   gsl_matrix_view m;
//   gsl_vector *eval;
//   gsl_matrix *evec;
//   gsl_vector_view evec_i;
//   gsl_eigen_symmv_workspace *w;
//   lmorb=p[0].lmerb;
//   A=(double*)calloc(nbreak*nbreak*lmorb*lmorb,sizeof(double));
//   local_n=(int)(floor(nbreak*lmorb/numprocs));
//   if(id<nbreak*lmorb-numprocs*local_n){
//     local_n=local_n+1;
//   }
//   for(i=0;i<local_n;i++){
//     li=floor((i+id*local_n)/nbreak); 
//     i1=i+id*local_n-nbreak*li;
//     for(j=0;j<nbreak*lmorb;j++){
//       lf=floor(j/nbreak);
//       j1=j-nbreak*lf;
//       if(li==lf){
// 	A[j+lmorb*nbreak*(i+id*local_n)]=creal(KE[j1+nbreak*i1+nbreak*nbreak*p[0].lo[li]]);
// 	if(i1==j1){
// 	  A[j+lmorb*nbreak*(i+id*local_n)]=A[j+lmorb*nbreak*(i+id*local_n)]-Zcharge/creal(Xz[i1]);
// 	}
//       }
//       for(io=0;io<norb;io++){
// 	A[j+lmorb*nbreak*(i+id*local_n)]=A[j+lmorb*nbreak*(i+id*local_n)]+creal(2.0*p[io].J[j+lmorb*nbreak*(i+id*local_n)]-p[io].K[j+lmorb*nbreak*(i+id*local_n)]);
//       }
//     }
//   }
//   MPI_Allgather(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,A,local_n*nbreak*lmorb,MPI_DOUBLE,MPI_COMM_WORLD);

//   for(io=0;io<norb;io++){
//     B=(double*)calloc(lmorb*nbreak,sizeof(double));
//     for(i=0;i<local_n;i++){
//       z1=0.0;
//       for(j=0;j<lmorb*nbreak;j++){
// 	z1=z1+A[j+lmorb*nbreak*(i+id*local_n)]*creal(p[io].cilm[j]);
//       }
//       B[i+id*local_n]=z1;
//     }
//     MPI_Allgather(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,B,local_n,MPI_DOUBLE,MPI_COMM_WORLD);
//     if(id==0){
//       z1=0.0;
//       for(i=0;i<lmorb*nbreak;i++){
// 	z1=z1+B[i]*creal(p[io].cilm[i]);
//       }
//       printf("io %d Energy=% .8le\n",io,z1);
//     }
//     free(B);
//   }
  
//   if(id==0){
//     m=gsl_matrix_view_array(A,nbreak*lmorb,nbreak*lmorb);
//     eval=gsl_vector_alloc(nbreak*lmorb);
//     evec=gsl_matrix_alloc(nbreak*lmorb,nbreak*lmorb);
//     w=gsl_eigen_symmv_alloc(nbreak*lmorb);
//     gsl_eigen_symmv(&m.matrix,eval,evec,w);
//     gsl_eigen_symmv_free(w);
//     gsl_eigen_symmv_sort(eval,evec,GSL_EIGEN_SORT_VAL_ASC);
//     for(i=0;i<4;i++){
//       z1=gsl_vector_get(eval,i);
//       printf("Eigenvalue %d % .8le\n",i,z1);
//     }
//     gsl_vector_free(eval);
//     gsl_matrix_free(evec);
//   }
//   free(A);
// }

void GRID::_HH1(double Zcharge,int id,int numprocs,int norb,orbitals *p,cmplx *NE,cmplx *H1){
  int i,j,i1,j1,li,lf,local_n,lmorb,io;
  double z1,*A,*B;
  gsl_matrix_view m;
  gsl_vector *eval;
  gsl_matrix *evec;
  gsl_vector_view evec_i;
  gsl_eigen_symmv_workspace *w;
  lmorb=p[0].lmerb;
  A=(double*)calloc(nbreak*nbreak*lmorb*lmorb,sizeof(double));
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
	A[j+lmorb*nbreak*(i+id*local_n)]=creal(KE[j1+nbreak*i1+nbreak*nbreak*p[2].lo[li]]);
      }
      for(io=0;io<norb;io++){
	A[j+lmorb*nbreak*(i+id*local_n)]=A[j+lmorb*nbreak*(i+id*local_n)]+2.0*creal(p[io].J[j+lmorb*nbreak*(i+id*local_n)])-creal(p[io].K[j+lmorb*nbreak*(i+id*local_n)]);
      }
      A[j+lmorb*nbreak*(i+id*local_n)]=A[j+lmorb*nbreak*(i+id*local_n)]+creal(NE[j+lmorb*nbreak*(i+id*local_n)]);
    }
  }
  MPI_Allgather(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,A,local_n*nbreak*lmorb,MPI_DOUBLE,MPI_COMM_WORLD);

  for(io=0;io<norb;io++){
    B=(double*)calloc(lmorb*nbreak,sizeof(double));
    for(i=0;i<local_n;i++){
      z1=0.0;
      for(j=0;j<lmorb*nbreak;j++){
	z1=z1+A[j+lmorb*nbreak*(i+id*local_n)]*creal(p[io].cilm[j]);
      }
      B[i+id*local_n]=z1;
    }
    MPI_Allgather(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,B,local_n,MPI_DOUBLE,MPI_COMM_WORLD);
    if(id==0){
      z1=0.0;
      for(i=0;i<lmorb*nbreak;i++){
	z1=z1+B[i]*creal(p[io].cilm[i]);
      }
      printf("Energy=% .8le\n",z1);
    }
    free(B);
  }
    
  if(id==0){
    m=gsl_matrix_view_array(A,nbreak*lmorb,nbreak*lmorb);
    eval=gsl_vector_alloc(nbreak*lmorb);
    evec=gsl_matrix_alloc(nbreak*lmorb,nbreak*lmorb);
    w=gsl_eigen_symmv_alloc(nbreak*lmorb);
    gsl_eigen_symmv(&m.matrix,eval,evec,w);
    gsl_eigen_symmv_free(w);
    gsl_eigen_symmv_sort(eval,evec,GSL_EIGEN_SORT_VAL_ASC);
    for(i=0;i<4;i++){
      z1=gsl_vector_get(eval,i);
      printf("Eigenvalue %d % .8le\n",i,z1);
    }
    gsl_vector_free(eval);
    gsl_matrix_free(evec);
  }
  free(A);
}

void GRID::_JMK(int id,int numprocs,int io,orbitals *p){
  int i,j,i1,j1,li,lf,k,local_n,lmorb;
  double max;
  cmplx z1,*A;
  lmorb=p[io].lmerb;
  // A=(cmplx*)calloc(lmorb*nbreak,sizeof(cmplx));
  // local_n=(int)(floor(nbreak*lmorb/numprocs));
  // if(id<nbreak*lmorb-numprocs*local_n){
  //   local_n=local_n+1;
  // }
  // for(i=0;i<local_n;i++){
  //   z1=0.0+0.0*I;
  //   for(j=0;j<lmorb*nbreak;j++){
  //     z1=z1+(p[io].J[j+lmorb*nbreak*(i+id*local_n)]-p[io].K[j+lmorb*nbreak*(i+id*local_n)])*p[io].cilm[j];
  //   }
  //   A[i+id*local_n]=z1;
  // }
  // MPI_Allgather(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,A,local_n,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD);
  if(id==0){
    // gsl_vector *v=gsl_vector_alloc(lmorb*nbreak);
    // for(i=0;i<lmorb*nbreak;i++){
    //   gsl_vector_set(v,i,creal(A[i]));      
    // }
    // max=gsl_vector_max(v);
    // gsl_vector_free(v);
    // printf("Max((J-K)chi)=% .8le\n",max);
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
  //  free(A);
}

cmplx GRID::_PSSN(int i,int j,int v){
  cmplx beta;
  beta=(2.0*v+1.0)*KEI[j+nbreak*i+nbreak*nbreak*v]/(Xz[j]*Xz[i]*cpow(Wz[i]*Wz[j],0.5));
  beta=beta+cpow(Xz[j]*Xz[i],v)/cpow(Xz[nbreak-1],2.0*v+1.0);
  return beta;
}

void GRID::_KO(int id,int numprocs,int lmerb,int *le,int lmorb,int *lo,cmplx *cilm,double *C3j,double *C3i,cmplx *K){
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
  MPI_Allgather(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,K,local_n*nbreak*lmerb,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD);
  if(id==0){
    end_time=MPI_Wtime();
    printf("% .8le\n",end_time-start_time);
    fflush(stdout);
  }  
  free(la);
  free(ma);
}

void GRID::_JO(int id,int numprocs,int lmerb,int *le,int lmorb,int *lo,cmplx *cilm,double *C3j,double *C3i,cmplx *J){
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
  MPI_Allgather(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,J,local_n*nbreak*lmerb,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD);
  if(id==0){
    end_time=MPI_Wtime();
    printf("% .8le\n",end_time-start_time);
    fflush(stdout);
  }  
  free(la);
  free(ma);
}


void GRID::_NE(int id,int numprocs,int lmerb,int *le,double *C3j,double *C3i,cmplx *NE){
  FILE *file;
  int i,j,k,l,nsph,i1,i2,j1,j2,l1,l2,m1,m2,li,lf,local_n;
  int lmax,nchannels,*la,*ma,ncent,*Z;
  double start_time,end_time,**center_xyz,yi,xi,*R,theta,phi;
  cmplx beta,z1,z2,z3,z4,z5,***In;
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
  file=fopen("center.inp","r");
  for(i=0;i<ncent;i++){
    fscanf(file,"%s%d%d",str,&k,&Z[i]);
    for(j=0;j<3;j++){
      fscanf(file,"%le",&center_xyz[i][j]);
      center_xyz[i][j]=center_xyz[i][j]*AnsgAu;
    }
  }
  fclose(file);

  R=(double*)calloc(ncent,sizeof(double));
  In=(cmplx***)calloc(ncent,sizeof(cmplx));
  for(i=0;i<ncent;i++){
    In[i]=(cmplx**)calloc(nbreak,sizeof(cmplx));
    for(j=0;j<nbreak;j++){
      In[i][j]=(cmplx*)calloc(lmax,sizeof(cmplx));
    }
  }
  //only for diatomics in works for H2O
  for(i=0;i<ncent;i++){
    R[i]=sqrt(pow(center_xyz[i][0],2)+pow(center_xyz[i][1],2)+pow(center_xyz[i][2],2));
    for(l=0;l<lmax;l++){
      for(k=0;k<nbreak;k++){
	xi=creal(Xz[k]);
	if(xi<R[i]){
	  yi=pow(xi/R[i],l)*(1.0/R[i]);
	}
	else{
	  yi=pow(R[i]/xi,l)*(1.0/xi);
	}
	In[i][k][l]=-Z[i]*yi;
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
  MPI_Allgather(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,NE,local_n*nbreak*lmerb,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD);
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

  // /////esto es para pintar orbitales en el eje z make sense only for diatomics
  // strcpy(fichero,"orbital_");
  // sprintf(sto,"%d",lo);
  // strcat(fichero,sto);
  // strcat(fichero,".inp");
  // file=fopen(fichero,"w");
  // for(i=nbreak-1;i>=0;i--){
  //   for(j=0;j<n;j++){
  //     theta=acos(zi[j]);
  //     phi=atan2(yi[j],xi[j]);
  //     beta=0.0;
  //     for(l=0;l<nchannels;l++){
  // 	beta=beta+(*cilm)[i+nbreak*l]*_Blm(str,(*li)[l],(*mi)[l],theta,phi)/(sqrt(creal(Wz[i]))*creal(Xz[i]));
  //   }
  //     fprintf(file,"% .8le % .8le % .8le\n",creal(Xz[i])*sin(theta)*cos(phi),creal(Xz[i])*cos(theta),creal(beta));
  //   }
  // }
  // for(i=0;i<nbreak;i++){
  //   for(j=0;j<n;j++){
  //     theta=acos(zi[j]);
  //     phi=atan2(yi[j],xi[j]);
  //     beta=0.0;
  //     for(l=0;l<nchannels;l++){
  // 	beta=beta+(*cilm)[i+nbreak*l]*_Blm(str,(*li)[l],(*mi)[l],theta,phi)/(sqrt(creal(Wz[i]))*creal(Xz[i]));
  //     }
  //     fprintf(file,"% .8le % .8le % .8le\n",creal(Xz[i])*sin(theta)*cos(phi),creal(Xz[i])*cos(theta),creal(beta));
  //   }
  // }  
  // fclose(file);
  
  strcpy(fichero,"orbital_");
  sprintf(sto,"%d",lo);
  strcat(fichero,sto);
  strcat(fichero,".inp");
  file=fopen(fichero,"w");
  for(i=nbreak-1;i>=0;i--){    
    beta=0.0;
    for(l=0;l<nchannels;l++){
      beta=beta+(*cilm)[i+nbreak*l]*_Blm(str,(*li)[l],(*mi)[l],M_PI,0.0)/(cpow(Wz[i],0.5)*Xz[i]);
    }
    fprintf(file,"% .8le % .8le % .8le\n",-creal(Xz[i]),creal(beta),cimag(beta));
  }
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


