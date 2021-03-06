#include <complex.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_vector.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "class_grid.h"
#include "orb.h"
#include "sasha.h"
#include "sphere_lebedev_rule.h"
#include "utils.h"

GRID::GRID(int nbas, int lamax, cmplx *Xi, cmplx *Wi, cmplx *TXX, cmplx *TIXX) {
  int i;
  nbreak = nbas;
  nmlps = lamax;
  Xz = (cmplx *)calloc(nbreak, sizeof(cmplx));
  Wz = (cmplx *)calloc(nbreak, sizeof(cmplx));
  for (i = 0; i < nbreak; i++) {
    Xz[i] = Xi[i];
    Wz[i] = Wi[i];
  }
  KE = (cmplx *)calloc(nbreak * nbreak * nmlps, sizeof(cmplx));
  KEI = (cmplx *)calloc(nbreak * nbreak * nmlps, sizeof(cmplx));
  for (i = 0; i < nbreak * nbreak * nmlps; i++) {
    KE[i] = TXX[i];
    KEI[i] = TIXX[i];
  }
}

void GRID::_HH1(double Zcharge, int id, int numprocs, int lmorb, int *lo,
                int *mo, cmplx *K, cmplx *J, cmplx *H1) {
  int i, j, i1, j1, li, lf, local_n;
  double z1, *A;
  gsl_matrix_view m;
  gsl_vector *eval;
  gsl_matrix *evec;
  gsl_vector_view evec_i;
  gsl_eigen_symmv_workspace *w;
  A = (double *)calloc(nbreak * nbreak * lmorb * lmorb, sizeof(double));
  local_n = (int)(floor(nbreak * lmorb / numprocs));
  if (id < nbreak * lmorb - numprocs * local_n) {
    local_n = local_n + 1;
  }

  for (i = 0; i < local_n; i++) {
    li = floor((i + id * local_n) / nbreak);  // calculate local indexes
    i1 = i + id * local_n - nbreak * li;
    for (j = 0; j < nbreak * lmorb; j++) {
      lf = floor(j / nbreak);
      j1 = j - nbreak * lf;
      if (li == lf) {
        A[j + lmorb * nbreak * (i + id * local_n)] =
            creal(KE[j1 + nbreak * i1 + nbreak * nbreak * lo[li]]);
        if (i1 == j1) {
          A[j + lmorb * nbreak * (i + id * local_n)] =
              A[j + lmorb * nbreak * (i + id * local_n)] -
              Zcharge / creal(Xz[i1]);
        }
      }
      A[j + lmorb * nbreak * (i + id * local_n)] =
          A[j + lmorb * nbreak * (i + id * local_n)] +
          creal(2.0 * J[j + lmorb * nbreak * (i + id * local_n)] -
                K[j + lmorb * nbreak * (i + id * local_n)]);
    }
  }
  MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, A, local_n * nbreak * lmorb,
                MPI_DOUBLE, MPI_COMM_WORLD);
  if (id == 0) {
    m = gsl_matrix_view_array(A, nbreak * lmorb, nbreak * lmorb);
    eval = gsl_vector_alloc(nbreak * lmorb);
    evec = gsl_matrix_alloc(nbreak * lmorb, nbreak * lmorb);
    w = gsl_eigen_symmv_alloc(nbreak * lmorb);
    gsl_eigen_symmv(&m.matrix, eval, evec, w);
    gsl_eigen_symmv_free(w);
    gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_ASC);
    for (i = 0; i < nbreak * lmorb; i++) {
      z1 = gsl_vector_get(eval, i);
      printf("%d % .8le\n", i, z1);
    }
    gsl_vector_free(eval);
    gsl_matrix_free(evec);
  }
  free(A);
}

void GRID::_JMK(int id, int numprocs, int lmorb, cmplx *cilm, cmplx *K,
                cmplx *J) {
  int i, j, i1, j1, li, lf, local_n;
  cmplx z1, *A;
  A = (cmplx *)calloc(lmorb * nbreak, sizeof(cmplx));
  local_n = (int)(floor(nbreak * lmorb / numprocs));
  if (id < nbreak * lmorb - numprocs * local_n) {
    local_n = local_n + 1;
  }
  for (i = 0; i < local_n; i++) {
    z1 = 0.0;
    for (j = 0; j < lmorb * nbreak; j++) {
      z1 = z1 + (J[j + lmorb * nbreak * (i + id * local_n)] -
                 K[j + lmorb * nbreak * (i + id * local_n)]) *
                    cilm[j];
    }
    A[i + id * local_n] = z1;
  }
  MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, A, local_n,
                MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD);
  if (id == 0) {
    gsl_vector *v = gsl_vector_alloc(lmorb * nbreak);
    for (i = 0; i < lmorb * nbreak; i++) {
      gsl_vector_set(v, i, creal(A[i]));
    }
    z1 = gsl_vector_max(v);
    printf("Max=% .8le\n", z1);
    for (i = 0; i < lmorb * nbreak; i++) {
      for (j = 0; j < lmorb * nbreak; j = j + 60) {
        printf(
            "% .8le % .8le % .8le % .8le\n", creal(K[i + nbreak * lmorb * j]),
            creal(K[j + nbreak * lmorb * i]), creal(J[i + nbreak * lmorb * j]),
            creal(J[j + nbreak * lmorb * i]));
      }
    }
  }
  free(A);
}

cmplx GRID::_PSSN(int i, int j, int v) {
  cmplx beta;
  beta = (2.0 * v + 1.0) * KEI[j + nbreak * i + nbreak * nbreak * v] /
         (Xz[j] * Xz[i] * cpow(Wz[i] * Wz[j], 0.5));
  beta = beta + cpow(Xz[j] * Xz[i], v) / cpow(Xz[nbreak - 1], 2.0 * v + 1.0);
  return beta;
}

void GRID::_KO(int id, int numprocs, char *str, int lmorb, int *lo, int *mo,
               cmplx *cilm, cmplx *K) {
  int i, j, k, l, nsph, i1, i2, j1, j2, l1, l2, m1, m2, li, lf, local_n;
  int lmax, nchannels, *la, *ma;
  double *C3j, *C3i, start_time, end_time;
  cmplx beta, z1, z2, z3, z4, z5;
  //  nsph=nmlps*(2*nmlps+1);
  nsph = (nmlps - 1) * (nmlps + 1) + 1;
  C3j = (double *)calloc(lmorb * lmorb * nsph, sizeof(double));
  C3i = (double *)calloc(lmorb * lmorb * nsph, sizeof(double));
  la = (int *)calloc(nsph, sizeof(int));
  ma = (int *)calloc(nsph, sizeof(int));
  k = 0;
  for (j = 0; j < nmlps; j++) {
    for (i = -j; i <= j; i++) {
      la[k] = j;
      ma[k] = i;
      k++;
    }
  }
  //  strcpy(str,str);
  for (i = 0; i < lmorb; i++) {
    for (j = 0; j < lmorb; j++) {
      for (k = 0; k < nsph; k++) {
        C3j[k + nsph * (j + lmorb * i)] =
            _C3jBlm(str, str, lo[i], lo[j], la[k], mo[i], mo[j], ma[k]);
        C3i[k + nsph * (j + lmorb * i)] =
            _C3jBlm(str, str, lo[i], lo[j], la[k], mo[i], mo[j], -ma[k]);
      }
    }
  }
  local_n = (int)(floor(nbreak * lmorb / numprocs));
  if (id < nbreak * lmorb - numprocs * local_n) {
    local_n = local_n + 1;
  }
  if (id == 0) {
    start_time = MPI_Wtime();
  }
  for (i = 0; i < local_n; i++) {
    li = floor((i + id * local_n) / nbreak);  // calculate local indexes
    i1 = i + id * local_n - nbreak * li;
    for (j = 0; j < nbreak * lmorb; j++) {  // make the upper diagonal matrix
      lf = floor(j / nbreak);
      j1 = j - nbreak * lf;
      z1 = 0.0;
      for (l = 0; l < nsph; l++) {
        z2 = 4.0 * M_PI * pow(-1, ma[l]) * _PSSN(i1, j1, la[l]) /
             (2.0 * la[l] + 1.0);
        for (l1 = 0; l1 < lmorb; l1++) {
          z3 = C3j[l + nsph * (l1 + lmorb * li)];
          for (l2 = 0; l2 < lmorb; l2++) {
            z4 = C3i[l + nsph * (l2 + lmorb * lf)];
            z5 = cilm[i1 + nbreak * l1] * cilm[j1 + nbreak * l2];
            z1 = z1 + z2 * z3 * z4 * z5;
          }
        }
      }
      K[j + lmorb * nbreak * (i + id * local_n)] = z1;
    }
  }
  MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, K, local_n * nbreak * lmorb,
                MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD);
  if (id == 0) {
    end_time = MPI_Wtime();
    printf("% .8le\n", end_time - start_time);
    fflush(stdout);
  }
  free(la);
  free(ma);
}

void GRID::_JO(int id, int numprocs, char *str, int lmorb, int *lo, int *mo,
               cmplx *cilm, cmplx *J) {
  int i, j, k, l, nsph, i1, i2, j1, j2, l1, l2, m1, m2, li, lf, local_n;
  int lmax, nchannels, *la, *ma;
  double *C3j, *C3i, start_time, end_time;
  cmplx beta, z1, z2, z3, z4, z5;
  //  nsph=nmlps*(2*nmlps+1);
  nsph = (nmlps - 1) * (nmlps + 1) + 1;
  C3j = (double *)calloc(lmorb * lmorb * nsph, sizeof(double));
  C3i = (double *)calloc(lmorb * lmorb * nsph, sizeof(double));
  la = (int *)calloc(nsph, sizeof(int));
  ma = (int *)calloc(nsph, sizeof(int));
  k = 0;
  for (j = 0; j < nmlps; j++) {
    for (i = -j; i <= j; i++) {
      la[k] = j;
      ma[k] = i;
      k++;
    }
  }
  for (i = 0; i < lmorb; i++) {
    for (j = 0; j < lmorb; j++) {
      for (k = 0; k < nsph; k++) {
        C3j[k + nsph * (j + lmorb * i)] =
            _C3jBlm(str, str, lo[i], lo[j], la[k], mo[i], mo[j], ma[k]);
        C3i[k + nsph * (j + lmorb * i)] =
            _C3jBlm(str, str, lo[i], lo[j], la[k], mo[i], mo[j], -ma[k]);
      }
    }
  }
  local_n = (int)(floor(nbreak * lmorb / numprocs));
  if (id < nbreak * lmorb - numprocs * local_n) {
    local_n = local_n + 1;
  }
  if (id == 0) {
    start_time = MPI_Wtime();
  }
  for (i = 0; i < local_n; i++) {
    li = floor((i + id * local_n) / nbreak);  // calculate local indexes
    i1 = i + id * local_n - nbreak * li;
    for (j = 0; j < nbreak * lmorb; j++) {  // make the upper diagonal matrix
      lf = floor(j / nbreak);
      j1 = j - nbreak * lf;
      z1 = 0.0;
      if (i1 == j1) {
        for (k = 0; k < nbreak; k++) {
          for (l = 0; l < nsph; l++) {
            z2 = 4.0 * M_PI * pow(-1, ma[l]) * _PSSN(i1, k, la[l]) /
                 (2.0 * la[l] + 1.0);
            z3 = C3j[l + nsph * (li + lmorb * lf)];
            for (l1 = 0; l1 < lmorb; l1++) {
              for (l2 = 0; l2 < lmorb; l2++) {
                z4 = C3i[l + nsph * (l1 + lmorb * l2)];
                z5 = cilm[k + nbreak * l1] * cilm[k + nbreak * l2];
                z1 = z1 + z2 * z3 * z4 * z5;
              }
            }
          }
        }
      }
      J[j + lmorb * nbreak * (i + id * local_n)] = z1;
    }
  }
  MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, J, local_n * nbreak * lmorb,
                MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD);
  if (id == 0) {
    end_time = MPI_Wtime();
    printf("% .8le\n", end_time - start_time);
    fflush(stdout);
  }
  free(la);
  free(ma);
}

void GRID::_GETORB(char *orbfile, char *symm, char *irrep, int lo, char *str1,
                   int *lmorb, int **li, int **mi, cmplx **cilm) {
  FILE *file;
  int i, j, k, l, ngrid, v, nchannels;
  double *orb, *x, *y, *z, *xi, *yi, *zi, *wi, theta, phi;
  cmplx beta;
  char str[10], fichero[25] = "angular_";
  const int n = 110;  // el que mejor funciona es 194
  nchannels = 0;
  l = (int)(ceil(nmlps / 2));
  _genlm(symm, irrep, l);
  strcat(fichero, symm);
  strcat(fichero, "_");
  strcat(fichero, irrep);
  strcat(fichero, ".inp");
  file = fopen(fichero, "r");
  for (i = getc(file); i != EOF; i = getc(file)) {
    if (i == '\n') {
      nchannels = nchannels + 1;
    }
  }
  fclose(file);
  *lmorb = nchannels;
  printf("%d\n", nchannels);
  *li = (int *)calloc(nchannels, sizeof(int));
  *mi = (int *)calloc(nchannels, sizeof(int));
  file = fopen(fichero, "r");
  for (i = 0; i < nchannels; i++) {
    fscanf(file, "%d%d%s", &(*li)[i], &(*mi)[i], str);
  }
  fclose(file);
  strcpy(str1, str);

  xi = (double *)calloc(n, sizeof(double));
  yi = (double *)calloc(n, sizeof(double));
  zi = (double *)calloc(n, sizeof(double));
  wi = (double *)calloc(n, sizeof(double));
  ld0110(xi, yi, zi, wi);
  ngrid = n * nbreak;
  orb = (double *)calloc(ngrid, sizeof(double));
  x = (double *)calloc(ngrid, sizeof(double));
  y = (double *)calloc(ngrid, sizeof(double));
  z = (double *)calloc(ngrid, sizeof(double));
  for (i = 0; i < nbreak; i++) {
    for (j = 0; j < n; j++) {
      x[j + n * i] = xi[j] * creal(Xz[i]);
      y[j + n * i] = yi[j] * creal(Xz[i]);
      z[j + n * i] = zi[j] * creal(Xz[i]);
    }
  }
  _getorb(orbfile, lo, ngrid, x, y, z, orb);

  *cilm = (cmplx *)calloc(nbreak * nchannels, sizeof(cmplx));
  for (l = 0; l < nchannels; l++) {
    for (i = 0; i < nbreak; i++) {
      beta = 0.0;
      for (j = 0; j < n; j++) {
        theta = acos(zi[j]);
        phi = atan2(yi[j], xi[j]);
        if (strcmp(symm, "ylm") == 0) {
          beta = beta + sqrt(creal(Wz[i])) * wi[j] * orb[j + n * i] *
                            conj(_Ylm((*li)[l], (*mi)[l], theta, phi)) *
                            creal(Xz[i]);
        } else {
          beta = beta + sqrt(creal(Wz[i])) * wi[j] * orb[j + n * i] *
                            _Blm(str, (*li)[l], (*mi)[l], theta, phi) *
                            creal(Xz[i]);
        }
      }
      (*cilm)[i + nbreak * l] = 4.0 * M_PI * beta;
    }
  }
  beta = 0.0;
  for (i = 0; i < nbreak; i++) {
    for (l = 0; l < nchannels; l++) {
      beta = beta + ((*cilm)[i + nbreak * l]) * conj((*cilm)[i + nbreak * l]);
    }
  }
  printf("%s %d%s Norm=% .8le % .8le\n", symm, lo + 1, irrep, creal(beta),
         cimag(beta));

  // /////esto es para pintar orbitales en el eje z make sense only for
  // diatomics for(i=nbreak-1;i>=0;i--){
  //   beta=0.0;
  //   for(l=0;l<nchannels;l++){
  //     //
  //     beta=beta+(*cilm)[i+nbreak*l]*_Ylm((*li)[l],(*mi)[l],M_PI,0.0)/(sqrt(creal(Wz[i]))*creal(Xz[i]));
  //     beta=beta+(*cilm)[i+nbreak*l]*_Blm(str,(*li)[l],(*mi)[l],M_PI,0.0)/(sqrt(creal(Wz[i]))*creal(Xz[i]));
  //   }
  //   printf("% .8le % .8le % .8le\n",-creal(Xz[i]),creal(beta),cimag(beta));
  // }
  // for(i=0;i<nbreak;i++){
  //   beta=0.0;
  //   for(l=0;l<nchannels;l++){
  //     //
  //     beta=beta+(*cilm)[i+nbreak*l]*_Ylm((*li)[l],(*mi)[l],0.0,0.0)/(sqrt(creal(Wz[i]))*creal(Xz[i]));
  //     beta=beta+(*cilm)[i+nbreak*l]*_Blm(str,(*li)[l],(*mi)[l],0.0,0.0)/(sqrt(creal(Wz[i]))*creal(Xz[i]));
  //   }
  //   printf("% .8le % .8le % .8le\n",creal(Xz[i]),creal(beta),cimag(beta));
  // }

  free(xi);
  free(yi);
  free(zi);
  free(wi);
  free(x);
  free(y);
  free(z);
  free(orb);
}

