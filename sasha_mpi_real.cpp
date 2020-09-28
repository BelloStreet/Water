#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <mpi.h>
#include "utils.h"
#include "sphere_lebedev_rule.h"
#include <gsl/gsl_sf_legendre.h>
#include "sasha.h"

double _Blm(char *type,int l,int m,double theta,double phi){
  int m1;
  double y,c1;
  m1=abs(m);
  if(m!=0){
    y=sqrt(2.0)*gsl_sf_legendre_sphPlm(l,m1,cos(theta));
  }
  else{
    y=gsl_sf_legendre_sphPlm(l,m1,cos(theta));
  }
  if(strcmp(type,"c")==0){
    c1=y*cos(m1*phi);
  }
  if(strcmp(type,"s")==0){
    c1=y*sin(m1*phi);
  }
  return c1;
}


void _genlm(char *symm,char *irrep,int l){
  FILE *file;
  int i,m,m0,l0,j;
  char cs[2],fichero[25]="angular_";
  strcat(fichero,symm);
  strcat(fichero,"_");
  strcat(fichero,irrep);
  strcat(fichero,".inp");
  //C2v
  if(strcmp(symm,"C2v")==0){
    if(strcmp(irrep,"a1")==0){
      m0=0;
      strcpy(cs,"c");
    }
    if(strcmp(irrep,"a2")==0){
      m0=2;
      strcpy(cs,"s");
    }
    if(strcmp(irrep,"b1")==0){
      m0=1;
      strcpy(cs,"c");
    }
    if(strcmp(irrep,"b2")==0){
      m0=1;
      strcpy(cs,"s");
    }
    file=fopen(fichero,"w");
    for(i=0;i<l;i++){
      if(m0<=i){
	fprintf(file,"%d %d %s\n",i,m0,cs);
      }
      // for(m=m0;m<=i;m=m+2){
      //  	fprintf(file,"%d %d %s\n",i,m,cs);
      // }
    }
    fclose(file);
  }
  ///D2h
  if(strcmp(symm,"D2h")==0){
    if(strcmp(irrep,"ag")==0){
      l0=0;
      m0=0;
      strcpy(cs,"c");
    }
    if(strcmp(irrep,"b1g")==0){
      l0=0;
      m0=2; //      m0=0;
      strcpy(cs,"s");
    }
    if(strcmp(irrep,"b2g")==0){
      l0=2;
      m0=1;
      strcpy(cs,"c");
    }
    if(strcmp(irrep,"b3g")==0){
      l0=2;
      m0=1;
      strcpy(cs,"s");
    }
    if(strcmp(irrep,"au")==0){
      l0=1;
      m0=2; //      m0=0;
      strcpy(cs,"s");
    }
    if(strcmp(irrep,"b1u")==0){
      l0=1;
      m0=0;
      strcpy(cs,"c");
    }
    if(strcmp(irrep,"b2u")==0){
      l0=1;
      m0=1;
      strcpy(cs,"s");
    }
    if(strcmp(irrep,"b3u")==0){
      l0=1;
      m0=1;
      strcpy(cs,"c");
    }
    file=fopen(fichero,"w");
    for(i=l0;i<l;i=i+2){
      for(m=m0;m<=i;m=m+2){
      	fprintf(file,"%d %d %s\n",i,m,cs);
      }
    }
    fclose(file);
  }
  //complex shperical harmonics
  if(strcmp(symm,"ylm")==0 && strcmp(irrep,"ylm")==0){
    strcpy(cs,"y");
    file=fopen(fichero,"w");
    for(j=0;j<l;j++){
      for(i=-j;i<=j;i++){
	fprintf(file,"%d %d %s\n",j,i,cs);
      }
    }
    fclose(file);
  }
  
}

void _C3jBlm(char *str1,char *str2,int j1,int j2,int j3,int m1,int m2,int m3,cmplx *zeta,cmplx *beta){
  int i,j,k,l;
  double y,z,*xi,*yi,*zi,*wi,phi,theta;
  const int n=2354;//n=974;
  xi=(double*)calloc(n,sizeof(double));
  yi=(double*)calloc(n,sizeof(double));
  zi=(double*)calloc(n,sizeof(double));
  wi=(double*)calloc(n,sizeof(double));
  ld2354(xi,yi,zi,wi);
  *zeta=0.0+0.0*I;
  *beta=0.0+0.0*I;
  for(j=0;j<n;j++){
    theta=acos(zi[j]);
    phi=atan2(yi[j],xi[j]);
    if(strcmp(str1,"y")==0 && strcmp(str2,"y")==0){
      *zeta=*zeta+wi[j]*conj(_Ylm(j1,m1,theta,phi))*_Ylm(j2,m2,theta,phi)*_Ylm(j3,m3,theta,phi);
      *beta=*beta+wi[j]*conj(_Ylm(j1,m1,theta,phi))*_Ylm(j2,m2,theta,phi)*_Ylm(j3,-m3,theta,phi);
    }
    else{
      if(m3>0){
	*zeta=*zeta+wi[j]*_Blm(str1,j1,m1,theta,phi)*_Blm(str2,j2,m2,theta,phi)*_Blm("c",j3,m3,theta,phi);
	*beta=*beta+wi[j]*_Blm(str1,j1,m1,theta,phi)*_Blm(str2,j2,m2,theta,phi)*_Blm("s",j3,-m3,theta,phi);
      }
      if(m3<0){
	*zeta=*zeta+wi[j]*_Blm(str1,j1,m1,theta,phi)*_Blm(str2,j2,m2,theta,phi)*_Blm("s",j3,m3,theta,phi);
	*beta=*beta+wi[j]*_Blm(str1,j1,m1,theta,phi)*_Blm(str2,j2,m2,theta,phi)*_Blm("c",j3,-m3,theta,phi);

      }
      if(m3==0){
	*zeta=*zeta+wi[j]*_Blm(str1,j1,m1,theta,phi)*_Blm(str2,j2,m2,theta,phi)*_Blm("c",j3,m3,theta,phi);
        *beta=*beta+wi[j]*_Blm(str1,j1,m1,theta,phi)*_Blm(str2,j2,m2,theta,phi)*_Blm("c",j3,m3,theta,phi);
      }
    }
  }
  *zeta=4.0*M_PI*(*zeta);
  *beta=4.0*M_PI*(*beta);
  free(xi);
  free(yi);
  free(zi);
  free(wi);
}

void _GETC3J(int id,int lmax,char *str1,char *str2, int lmorb1,int lmorb2,int *lo1,int *mo1,int *lo2,int *mo2,cmplx **C3Kj,cmplx **C3Ki,cmplx **C3Jj,cmplx **C3Ji,int type){
  int i,j,k,lv,mv,nsph,*la,*ma,tag1=9,tag2=99;
  cmplx zeta,beta,*c3i,*c3j;
  MPI_Status status1,status2;
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
  *C3Kj=(cmplx*)calloc(lmorb1*lmorb2*nsph,sizeof(cmplx));
  *C3Ki=(cmplx*)calloc(lmorb1*lmorb2*nsph,sizeof(cmplx));
  *C3Jj=(cmplx*)calloc(lmorb1*lmorb1*nsph,sizeof(cmplx));
  *C3Ji=(cmplx*)calloc(lmorb2*lmorb2*nsph,sizeof(cmplx));

  c3j=(cmplx*)calloc(lmorb2*nsph,sizeof(cmplx));
  c3i=(cmplx*)calloc(lmorb2*nsph,sizeof(cmplx));    
  if(id<lmorb1){
    for(j=0;j<lmorb2;j++){      
      k=0;
      for(lv=0;lv<lmax;lv++){
	if(la[k]>=abs(lo1[id]-lo2[j]) && la[k]<=lo1[id]+lo2[j]){
	  for(mv=0;mv<=lv;mv++){  
	    if(mv==0){
	      _C3jBlm(str1,str2,lo1[id],lo2[j],la[k],mo1[id],mo2[j],ma[k],&zeta,&beta);
	      c3j[k+nsph*j]=zeta;
	      c3i[k+nsph*j]=beta;
	      k++;
	    }
	    else{
	      _C3jBlm(str1,str2,lo1[id],lo2[j],la[k],mo1[id],mo2[j],ma[k],&zeta,&beta);
	      c3j[k+nsph*j]=zeta;
	      c3i[k+nsph*j]=beta;
	      k++;
	      c3j[k+nsph*j]=beta;
	      c3i[k+nsph*j]=zeta;
	      k++;
	    }	    
	  }
	}
	else{
	  for(mv=0;mv<=lv;mv++){  
	    if(mv==0){
	      k++;
	    }
	    else{
	      k=k+2;
	    }
	  }
	}
      }
    }
  }
  if(id<lmorb1){
    MPI_Send(&(c3j[0]),lmorb2*nsph,MPI_DOUBLE_COMPLEX,57,tag1,MPI_COMM_WORLD);
    MPI_Send(&(c3i[0]),lmorb2*nsph,MPI_DOUBLE_COMPLEX,57,tag2,MPI_COMM_WORLD);
  }
  if(id==57){
    for(i=0;i<lmorb1;i++){
      MPI_Recv(&(c3j[0]),lmorb2*nsph,MPI_DOUBLE_COMPLEX,i,tag1,MPI_COMM_WORLD,&status1);
      MPI_Recv(&(c3i[0]),lmorb2*nsph,MPI_DOUBLE_COMPLEX,i,tag2,MPI_COMM_WORLD,&status2);
      for(j=0;j<lmorb2;j++){
	for(k=0;k<nsph;k++){
	  (*C3Kj)[k+nsph*(j+lmorb2*i)]=c3j[k+nsph*j];
	  (*C3Ki)[k+nsph*(j+lmorb2*i)]=c3i[k+nsph*j];
	}
      }
    }
  }
  free(c3j);
  free(c3i);  
  
  c3j=(cmplx*)calloc(lmorb1*nsph,sizeof(cmplx));
  if(id<lmorb1){
    for(j=0;j<lmorb1;j++){
      k=0;
      for(lv=0;lv<lmax;lv++){
	if(la[k]>=abs(lo1[id]-lo1[j]) && la[k]<=lo1[id]+lo1[j]){
	  for(mv=0;mv<=lv;mv++){  
	    if(mv==0){
	      _C3jBlm(str1,str1,lo1[id],lo1[j],la[k],mo1[id],mo1[j],ma[k],&zeta,&beta);
	      c3j[k+nsph*j]=zeta;
	      k++;
	    }
	    else{
	      _C3jBlm(str1,str1,lo1[id],lo1[j],la[k],mo1[id],mo1[j],ma[k],&zeta,&beta);
	      c3j[k+nsph*j]=zeta;
	      k++;
	      c3j[k+nsph*j]=beta;
	      k++;
	    }	    
	  }
	}
	else{
	  for(mv=0;mv<=lv;mv++){  
	    if(mv==0){
	      k++;
	    }
	    else{
	      k=k+2;
	    }	    
	  }
	}
      }
    }
  }
  if(id<lmorb1){
    MPI_Send(&(c3j[0]),lmorb1*nsph,MPI_DOUBLE_COMPLEX,57,tag1,MPI_COMM_WORLD);
  }
  if(id==57){
    for(i=0;i<lmorb1;i++){
      MPI_Recv(&(c3j[0]),lmorb1*nsph,MPI_DOUBLE_COMPLEX,i,tag1,MPI_COMM_WORLD,&status1);
      for(j=0;j<lmorb1;j++){
	for(k=0;k<nsph;k++){
	  (*C3Jj)[k+nsph*(j+lmorb1*i)]=c3j[k+nsph*j];
	}
      }
    }
  }
  free(c3j);

  c3i=(cmplx*)calloc(lmorb2*nsph,sizeof(cmplx));
  if(id<lmorb2){
    for(j=0;j<lmorb2;j++){
      k=0;
      for(lv=0;lv<lmax;lv++){
	if(la[k]>=abs(lo2[id]-lo2[j]) && la[k]<=lo2[id]+lo2[j]){
	  for(mv=0;mv<=lv;mv++){  	  
	    if(mv==0){
	      _C3jBlm(str2,str2,lo2[id],lo2[j],la[k],mo2[id],mo2[j],ma[k],&zeta,&beta);
	      c3i[k+nsph*j]=beta;
	      k++;
	    }
	    else{
	      _C3jBlm(str2,str2,lo2[id],lo2[j],la[k],mo2[id],mo2[j],ma[k],&zeta,&beta);
	      c3i[k+nsph*j]=beta;
	      k++;
	      c3i[k+nsph*j]=zeta;
	      k++;
	    }	    
	  }
	}
	else{
	  for(mv=0;mv<=lv;mv++){  	  
	    if(mv==0){
	      k++;
	    }
	    else{
	      k=k+2;
	    }	    
	  }
	}
      }
    }
  }
  if(id<lmorb2){
    MPI_Send(&(c3i[0]),lmorb2*nsph,MPI_DOUBLE_COMPLEX,57,tag2,MPI_COMM_WORLD);
  }
  if(id==57){
    for(i=0;i<lmorb2;i++){
      MPI_Recv(&(c3i[0]),lmorb2*nsph,MPI_DOUBLE_COMPLEX,i,tag2,MPI_COMM_WORLD,&status2);
      for(j=0;j<lmorb2;j++){
	for(k=0;k<nsph;k++){
	  (*C3Ji)[k+nsph*(j+lmorb2*i)]=c3i[k+nsph*j];
	}
      }
    }
  }
  free(c3i); 

  MPI_Bcast(&((*C3Kj)[0]),lmorb1*lmorb2*nsph,MPI_DOUBLE_COMPLEX,57,MPI_COMM_WORLD);
  MPI_Bcast(&((*C3Ki)[0]),lmorb1*lmorb2*nsph,MPI_DOUBLE_COMPLEX,57,MPI_COMM_WORLD);
  MPI_Bcast(&((*C3Ji)[0]),lmorb2*lmorb2*nsph,MPI_DOUBLE_COMPLEX,57,MPI_COMM_WORLD);
  MPI_Bcast(&((*C3Jj)[0]),lmorb1*lmorb1*nsph,MPI_DOUBLE_COMPLEX,57,MPI_COMM_WORLD);
  
  free(la);
  free(ma);  
}

