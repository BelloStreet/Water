#ifndef HEADER_SASHA_H
#define HEADER_SASHA_H

extern cmplx _Ylm(int l,int m,double theta,double phi);
extern double _Blm(char *type,int l,int m,double theta,double phi);
extern void _genlm(char *symm,char *irrep,int l);
extern void _C3jBlm(char *str1,char *str2,int j1,int j2,int j3,int m1,int m2,int m3,cmplx *zeta,cmplx *beta);
extern void _GETC3J(int id,int lmax,char *str1,char *str2, int lmorb1,int lmorb2,int *lo1,int *mo1,int *lo2,int *mo2,double **C3Kj,double **C3Ki,double **C3Jj,double **C3Ji,int type);


#endif
