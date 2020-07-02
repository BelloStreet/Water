#include <complex.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_spline.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "class_grid.h"
#include "orb.h"
#include "sasha.h"
#include "utils.h"
int main(int argc, char **argv) {
  FILE *file;
  int i, j, norb, l, tag = 99;
  int id = 0, numprocs, root = 0, nbas, lamax, lmax, io, *lo = NULL, *mo = NULL,
      lmorb;
  int *la, *ma, *le, *me, *no, lmerb = 0;
  double z1, z2, z3, z4, start_time, end_time;
  cmplx *Xi, *Wi, *TXX, *TIXX, *H1, *cilm = NULL, *K, *J, *NE, *CKj, *CKi, *CJj,
                                    *CJi;
  char command[150], fichero[50], fileel[50] = "angular_";
  char symm[5] = "ylm", irrep[5] = "ylm", str[2] = "y", **irreps, strs[2] = "c";
  orbitals od[5];
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &id);
  printf("id=%d numprocs=%d\n", id, numprocs);
  fflush(stdout);
  file = fopen("quadrature.dat", "r");
  fscanf(file, "%d", &lamax);
  fscanf(file, "%d", &nbas);
  Xi = (cmplx *)calloc(nbas, sizeof(cmplx));
  Wi = (cmplx *)calloc(nbas, sizeof(cmplx));
  for (i = 0; i < nbas; i++) {
    fscanf(file, "%le%le%le%le", &z1, &z2, &z3, &z4);
    Xi[i] = z1 + z2 * I;
    Wi[i] = z3 + z4 * I;
  }
  fclose(file);

  lmax = floor(lamax / 2);
  TXX = (cmplx *)calloc(nbas * nbas * lamax, sizeof(cmplx));
  TIXX = (cmplx *)calloc(nbas * nbas * lamax, sizeof(cmplx));
  file = fopen("kinetic.dat", "r");
  for (i = 0; i < nbas * nbas * lamax; i++) {
    fscanf(file, "%le%le", &z1, &z2);
    TXX[i] = z1 + z2 * I;
  }
  for (i = 0; i < nbas * nbas * lamax; i++) {
    fscanf(file, "%le%le", &z1, &z2);
    TIXX[i] = z1 + z2 * I;
  }
  fclose(file);
  GRID CORR(nbas, lamax, Xi, Wi, TXX, TIXX);
  free(Xi);
  free(Wi);
  free(TXX);
  free(TIXX);

  H1 = (cmplx *)calloc(nbas * nbas, sizeof(cmplx));

  file = fopen("orb.inp", "r");
  fscanf(file, "%s", fichero);
  fscanf(file, "%s%s", symm, irrep);
  fscanf(file, "%d", &norb);
  no = (int *)calloc(norb, sizeof(int));
  irreps = (char **)calloc(norb, sizeof(char *));
  for (i = 0; i < norb; i++) {
    irreps[i] = (char *)calloc(5, sizeof(char));
  }
  for (io = 0; io < norb; io++) {
    fscanf(file, "%d%s", &no[io], irreps[io]);
  }
  fclose(file);
  if (id == root) {
    _genlm(symm, irrep, lmax);
    strcat(fileel, symm);
    strcat(fileel, "_");
    strcat(fileel, irrep);
    strcat(fileel, ".inp");
    file = fopen(fileel, "r");
    for (i = getc(file); i != EOF; i = getc(file)) {
      if (i == '\n') {
        lmerb = lmerb + 1;
      }
    }
    fclose(file);
  }
  MPI_Bcast(&lmerb, 1, MPI_INT, root, MPI_COMM_WORLD);
  le = (int *)calloc(lmerb, sizeof(int));
  me = (int *)calloc(lmerb, sizeof(int));
  if (id == root) {
    file = fopen(fileel, "r");
    for (i = 0; i < lmerb; i++) {
      fscanf(file, "%d%d%s", &le[i], &me[i], str);
    }
    fclose(file);
  }
  MPI_Bcast(le, lmerb, MPI_INT, root, MPI_COMM_WORLD);
  MPI_Bcast(me, lmerb, MPI_INT, root, MPI_COMM_WORLD);
  MPI_Bcast(&str, strlen(str) + 1, MPI_CHAR, root, MPI_COMM_WORLD);
  if (id == root + 1) {
    printf("Candelones con tostones %s\n", str);
  }

  for (io = 0; io < norb; io++) {
    if (id == root) {
      CORR._GETORB(fichero, symm, irreps[io], no[io], strs, &lmorb, &lo, &mo,
                   &cilm);
    }
    MPI_Bcast(&lmorb, 1, MPI_INT, root, MPI_COMM_WORLD);
    MPI_Bcast(&strs, strlen(strs) + 1, MPI_CHAR, root, MPI_COMM_WORLD);
    if (id != root) {
      lo = (int *)calloc(lmorb, sizeof(int));
      mo = (int *)calloc(lmorb, sizeof(int));
      cilm = (cmplx *)calloc(nbas * lmorb, sizeof(cmplx));
    }
    MPI_Bcast(lo, lmorb, MPI_INT, root, MPI_COMM_WORLD);
    MPI_Bcast(mo, lmorb, MPI_INT, root, MPI_COMM_WORLD);
    MPI_Bcast(cilm, nbas * lmorb, MPI_DOUBLE_COMPLEX, root, MPI_COMM_WORLD);

    // J
    if (id == root) {
      start_time = MPI_Wtime();
    }
    _GETC3J(id, lmax, str, strs, lmerb, lmorb, le, me, lo, mo, &CKj, &CKi, &CJj,
            &CJi, 1);
    if (id == root) {
      end_time = MPI_Wtime();
      printf("tiempo CK % .8le\n", end_time - start_time);
    }

    J = (cmplx *)calloc(nbas * nbas * lmerb * lmerb, sizeof(cmplx));
    CORR._JO(id, numprocs, lmerb, le, lmorb, lo, cilm, CJj, CJi, J);
    free(CJi);
    free(CJj);

    K = (cmplx *)calloc(nbas * nbas * lmerb * lmerb, sizeof(cmplx));
    CORR._KO(id, numprocs, lmerb, le, lmorb, lo, cilm, CKj, CKi, K);
    free(CKi);
    free(CKj);

    _INITORB(io, nbas, lmerb, strs, lmorb, lo, mo, cilm, J, K, &(od[0]));
    CORR._JMK(id, numprocs, io, &(od[0]));

    free(K);
    free(J);
    free(lo);
    free(mo);
    free(cilm);
  }
  if (id == root) {
    start_time = MPI_Wtime();
  }
  if (id == root) {
    _ORTORB(nbas, 0, 1, &(od[0]));
    _ORTORB(nbas, 0, 3, &(od[0]));
    _ORTORB(nbas, 1, 3, &(od[0]));
  }
  _GETC3J(id, lamax, str, str, lmerb, lmerb, le, me, le, me, &CKj, &CKi, &CJj,
          &CJi, 2);
  if (id == root) {
    end_time = MPI_Wtime();
    printf("tiempo CK % .8le\n", end_time - start_time);
  }
  NE = (cmplx *)calloc(nbas * nbas * lmerb * lmerb, sizeof(cmplx));
  CORR._NE(id, numprocs, lmerb, le, CKj, CKi, NE);
  free(CKi);
  free(CKj);

  CORR._HH1(4.0, id, numprocs, norb, &(od[0]), NE, H1);
  // two electron hamiltonian
  free(NE);
  free(le);
  free(me);
  free(H1);

  MPI_Finalize();
  return 0;
}
