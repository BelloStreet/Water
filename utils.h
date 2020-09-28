#ifndef HEADER_UTILS_H
#define HEADER_UTILS_H

#define maxorb 31 //72
#define maxdiv 200

typedef _Complex double cmplx;

typedef struct datos{
  int nn;
  double *xx;
  double *yy;
}datos;

typedef struct obitals{
  cmplx J[maxdiv*maxdiv*maxorb*maxorb];
  cmplx K[maxdiv*maxdiv*maxorb*maxorb];
  cmplx cilm[maxorb*maxdiv];
  int lo[maxorb];
  int mo[maxorb];
  int lmorb;
  int lmerb;
  char str[5];
}orbitals;

typedef struct hamilton{
  cmplx c_nl[maxdiv*maxdiv*maxorb*maxorb];
  cmplx e_nl[maxdiv*maxorb];
  int le[maxorb];
  int me[maxorb];
  int lmerb;
  int nbas;
  char str[5];
}hamilton;



extern double _fcsplines(double z,struct datos *p);
extern cmplx _Ylm(int l,int m,double theta,double phi);

extern void _INITORB(int io,int nbas,int lmerb,char *str,int lmorb,int *lo,int *mo,cmplx *cilm,cmplx *J,cmplx *K,orbitals *p);

extern void _INITH(int io,int nbas,char *str,int lmerb,int *le,int *me,cmplx *c_nl,cmplx *e_nl,hamilton *p);

#endif
