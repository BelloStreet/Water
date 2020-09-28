#include "utils.h"
#include <string.h> 
#include <stdio.h> 
#include <stdlib.h> 
#include <math.h>
#include <complex.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_legendre.h>

double _fcsplines(double z,struct datos *p){
  int n;
  double yi=0.0,*x,*y;
  n=p->nn;
  x=p->xx;
  y=p->yy;
  if(z>=x[0] && z<=x[n-1]){
    gsl_interp_accel *acc=gsl_interp_accel_alloc();
    gsl_spline *spline=gsl_spline_alloc(gsl_interp_cspline,n);
    gsl_spline_init(spline,x,y,n);
    yi=gsl_spline_eval(spline,z,acc);
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
  }
  return yi;
}

cmplx _Ylm(int l,int m,double theta,double phi){
  int m1;
  double y;
  cmplx c1,c2,c3;
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


void _INITORB(int io,int nbas,int lmerb,char *str,int lmorb,int *lo,int *mo,cmplx *cilm,cmplx *J,cmplx *K,orbitals *p){
  int i;
  p[io].lmorb=lmorb;
  p[io].lmerb=lmerb;
  for(i=0;i<nbas*nbas*lmerb*lmerb;i++){
    p[io].J[i]=J[i];
    p[io].K[i]=K[i];
  }
  for(i=0;i<nbas*lmorb;i++){
    p[io].cilm[i]=cilm[i];
  }
  for(i=0;i<lmorb;i++){
    p[io].lo[i]=lo[i];
    p[io].mo[i]=mo[i];
  }
  strcpy(p[io].str,str); 
}

void _INITH(int io,int nbas,char *str,int lmerb,int *le,int *me,cmplx *c_nl,cmplx *e_nl,hamilton *p){
  int i;
  p[io].nbas=nbas;
  p[io].lmerb=lmerb;
  for(i=0;i<nbas*nbas*lmerb*lmerb;i++){
    p[io].c_nl[i]=c_nl[i];
  }
  for(i=0;i<nbas*lmerb;i++){
    p[io].e_nl[i]=e_nl[i];
  }
  for(i=0;i<lmerb;i++){
    p[io].le[i]=le[i];
    p[io].me[i]=me[i];
  }
  strcpy(p[io].str,str);
}

