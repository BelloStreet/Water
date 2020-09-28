#ifndef HEADER_CLASS_GRID_H
#define HEADER_CLASS_GRID_H

typedef _Complex double cmplx;

class GRID{
 private:
  int nbreak;
  int nmlps;
  cmplx *Xz;
  cmplx *Wz;
  cmplx *KE;
  cmplx *KEI;
 public:
  void _HH1(int id,int numprocs,int kk,int norb,char *str,int lmerb,int *le,int *me,orbitals *p,cmplx *NE,cmplx *c_nl,cmplx *e_nl);
  void _HH2(int id,int numprocs,int nmkl[],int nmkl0[],int nf,int mf,hamilton *p,double *AT,char *symm,int type);
  void _VVz(int id,int numprocs,int nmkl[],int nmkl0[],int nf,int mf,hamilton *p,double *AT,double *Vnuc);
  void _GETORB(char *orbfile,char *symm,char *irrep,int lo,char *str,int *lmorb,int **li,int **mi,cmplx **cilm);
  cmplx _PSSN(int i,int j,int v);
  cmplx _LGRN(int i,cmplx x);
  void _KO(int id,int numprocs,int lmerb,int *le,int lmorb,int *lo,cmplx *cilm,cmplx *C3j,cmplx *C3i,cmplx *K);
  void _JO(int id,int numprocs,int lmerb,int *le,int lmorb,int *lo,cmplx *cilm,cmplx *C3j,cmplx *C3i,cmplx *J);
  void _JMK(int id,int numprocs,int io,orbitals *p);
  void _NE(int id,int numprocs,int lmerb,int *le,cmplx *C3j,cmplx *C3i,cmplx *NE);
  GRID(int nbas,int lamax,cmplx *Xi,cmplx *Wi,cmplx *TXX,cmplx *TIXX);
};


#endif
