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
  void _HH1(double Zcharge,int id,int numprocs,int norb,orbitals *p,cmplx *NE,cmplx *H1);
  void _GETORB(char *orbfile,char *symm,char *irrep,int lo,char *str,int *lmorb,int **li,int **mi,cmplx **cilm);
  cmplx _PSSN(int i,int j,int v);
  void _KO(int id,int numprocs,int lmerb,int *le,int lmorb,int *lo,cmplx *cilm,cmplx *C3j,cmplx *C3i,cmplx *K);
  void _JO(int id,int numprocs,int lmerb,int *le,int lmorb,int *lo,cmplx *cilm,cmplx *C3j,cmplx *C3i,cmplx *J);
  void _JMK(int id,int numprocs,int io,orbitals *p);
  void _NE(int id,int numprocs,int lmerb,int *le,cmplx *C3j,cmplx *C3i,cmplx *NE);
  GRID(int nbas,int lamax,cmplx *Xi,cmplx *Wi,cmplx *TXX,cmplx *TIXX);
};


#endif
