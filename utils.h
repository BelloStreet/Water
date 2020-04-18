#ifndef HEADER_UTILS_H
#define HEADER_UTILS_H

#define maxorb 121 //72
#define maxdiv 116

typedef _Complex double cmplx;

typedef struct datos{
  int nn;
  double *xx;
  double *yy;
}datos;

typedef struct obitals{
  cmplx J[maxdiv*maxdiv*maxorb*maxorb];
  cmplx K[maxdiv*maxdiv*maxorb*maxorb];
  cmplx cilm[maxorb*maxdiv];//*cilm;
  int lo[maxorb];//*lo;
  int mo[maxorb];//*mo;
  int lmorb;
  int lmerb;
  char str[5];//*str;
}orbitals;


extern double _fcsplines(double z,struct datos *p);
extern cmplx _Ylm(int l,int m,double theta,double phi);

extern void _INITORB(int io,int nbas,int lmerb,char *str,int lmorb,int *lo,int *mo,cmplx *cilm,cmplx *J,cmplx *K,orbitals *p);

void _ORTORB(int nbas,int io,int il,orbitals *p);
#endif
