#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <complex.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_integration.h>
using namespace std;
#include "sphere_lebedev_rule.h"
double f(double x,double y,double z);
double B(int nbreak,int k,double a,double b,double xi,int i);
_Complex double Ylm(int l,int m,double theta,double phi);
int main(void){
  FILE *file;
  int n=146,i,j,N=100,l,m,nchannels=16,v;
  double z1,theta,phi,x1,y1;
  double alpha,*x,*y,*z,*w,*Xi,*Wi,xi,wi;
  int li[16]={0,1,1,1,2,2,2,2,2,3,3,3,3,3,3,3};
  int mi[16]={0,-1,0,1,-2,-1,0,1,2,-3,-2,-1,0,1,2,3};
  _Complex double **cilm,beta;
  gsl_integration_glfixed_table *t;

  Xi=(double*)calloc(N,sizeof(double));
  Wi=(double*)calloc(N,sizeof(double));
  t=gsl_integration_glfixed_table_alloc(N);
  for(i=0;i<N;i++){
    gsl_integration_glfixed_point(0.0,20.0,i,&xi,&wi,t);
    Xi[i]=xi;
    Wi[i]=wi;
  }
  gsl_integration_glfixed_table_free(t);
  
  x=(double*)calloc(n,sizeof(double));
  y=(double*)calloc(n,sizeof(double));
  z=(double*)calloc(n,sizeof(double));
  w=(double*)calloc(n,sizeof(double));
  ld0146(x,y,z,w);

  //  alpha=0.0;
  //   for(i=0;i<n;i++){
  //     alpha=alpha+w[i]*f(x[i],y[i],z[i]);
  //   }
  //   printf("% .8le\n",4.0*M_PI*alpha);  
  //integral de la funcion en coordenadas esfericas
  alpha=0.0;
  for(i=0;i<N;i++){
    for(j=0;j<n;j++){
      alpha=alpha+Wi[i]*w[j]*f(Xi[i]*x[j],Xi[i]*y[j],Xi[i]*z[j])*pow(Xi[i],2.0);
    }
  }
  printf("% .8le\n",4.0*M_PI*alpha);
  free(Xi);
  free(Wi);

  file=fopen("quadrature.dat","r");
  fscanf(file,"%d",&N);
  Xi=(double*)calloc(N,sizeof(double));
  Wi=(double*)calloc(N,sizeof(double));
  for(i=0;i<N;i++){
    fscanf(file,"%le%le%le%le",&Xi[i],&z1,&Wi[i],&z1);
  }
  fclose(file);

  //expansion in spherical coordinates
  cilm=(_Complex double**)calloc(N,sizeof(_Complex double));
  for(i=0;i<N;i++){
    cilm[i]=(_Complex double*)calloc(nchannels,sizeof(_Complex double));
  }
  for(i=0;i<N;i++){
    for(l=0;l<nchannels;l++){
      beta=0.0;
      for(j=0;j<n;j++){
	theta=acos(z[j]);
	phi=atan2(y[j],x[j]);
	beta=beta+sqrt(Wi[i])*w[j]*f(Xi[i]*x[j],Xi[i]*y[j],Xi[i]*z[j])*conj(Ylm(li[l],mi[l],theta,phi));
      }
      cilm[i][l]=4.0*M_PI*beta;
    }
  }

  beta=0.0;
  theta=0.4*M_PI;
  phi=0.7*M_PI;
  z1=Xi[15]*cos(theta);
  y1=Xi[15]*sin(theta)*sin(phi);
  x1=Xi[15]*sin(theta)*cos(phi);
  for(l=0;l<nchannels;l++){
    beta=beta+cilm[15][l]*Ylm(li[l],mi[l],theta,phi)/(sqrt(Wi[15]));
  }
  printf("% .8le % .8le % .8le\n",creal(beta),cimag(beta),f(x1,y1,z1));

  beta=0.0;
  theta=0.65*M_PI;
  phi=0.35*M_PI;
  z1=Xi[4]*cos(theta);
  y1=Xi[4]*sin(theta)*sin(phi);
  x1=Xi[4]*sin(theta)*cos(phi);
  for(l=0;l<nchannels;l++){
    beta=beta+cilm[4][l]*Ylm(li[l],mi[l],theta,phi)/(sqrt(Wi[4]));
  }
  printf("% .8le % .8le % .8le\n",creal(beta),cimag(beta),f(x1,y1,z1));

  free(cilm);
  free(x);
  free(y);
  free(z);
  free(w);
  free(Xi);
  free(Wi);
  return 0;
}

double f(double x,double y,double z){
  //  return 1.0+x+pow(y,2.0)+pow(x,2.0)*y+pow(x,4.0)+pow(y,5.0)+pow(x*y*z,2.0);
  //  return sqrt(pow(x,2.0)+pow(y,2.0)+pow(z,2.0));
  //  return exp(pow(pow(x,2.0)+pow(y,2.0)+pow(z,2.0),1.5));
  return exp(-sqrt(pow(x,2.0)+pow(y,2.0)+pow(z,2.0)))*x*y*z;
}

_Complex double Ylm(int l,int m,double theta,double phi){
  int m1;
  double y;
  _Complex double c1,c2,c3;
  m1=abs(m);
  y=pow(-1.0,m)*gsl_sf_legendre_sphPlm(l,m1,cos(theta));
  c1=cos(m1*phi)+sin(m1*phi)*I;
  if(m>=0){
    c3=c1*y;
  }
  else{
    c3=pow(-1.0,m)*conj(c1*y);
  }
  return c3;
}


