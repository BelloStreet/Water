#include "orb.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void _getorb(char *fichero, int l, int ngrid, double *xi, double *yi,
             double *zi, double *mo_orb) {
  FILE *file;
  int i, j, k, v, count = 0, ncent = 0, *pos, nlines, icent, n, nbasis, off,
                  n_orb, n_occu_orb;
  int *ns_icent, *np_icent, *nd_icent, *nf_icent, xgrid, ygrid, zgrid;
  int *aux_s_icent, *aux_p_icent, *aux_d_icent, *aux_f_icent;
  int **nprim_s_icent, **nprim_p_icent, **nprim_d_icent, **nprim_f_icent;
  double val1, val2, val, a, **center_xyz, **mo_occu, *occu, fact, S_xyz, x, y,
      z, *we;
  double **coeff_s, **expon_s, **coeff_p, **expon_p, **coeff_d, **expon_d,
      **coeff_f, **expon_f;
  char str[50], str1[200], str2[200], command[500];
  const double AnsgAu = 1.889726124565062;
  strcpy(command, "./orb.sh ");
  strcat(command, fichero);
  system(command);
  file = fopen("read.inp", "r");
  fscanf(file, "%d", &ncent);
  printf("Number of centers detected=%d\n", ncent);
  pos = (int *)calloc(ncent, sizeof(int));
  for (i = 0; i < ncent; i++) {
    fscanf(file, "%d", &pos[i]);
  }
  fscanf(file, "%d", &nlines);
  fscanf(file, "%d", &n_orb);
  fscanf(file, "%d", &n_occu_orb);
  fclose(file);

  center_xyz = (double **)calloc(ncent, sizeof(double));
  for (i = 0; i < ncent; i++) {
    center_xyz[i] = (double *)calloc(3, sizeof(double));
  }
  file = fopen("center.inp", "r");
  for (i = 0; i < ncent; i++) {
    fscanf(file, "%s%d%d", str, &k, &k);
    for (j = 0; j < 3; j++) {
      fscanf(file, "%le", &center_xyz[i][j]);
      center_xyz[i][j] = center_xyz[i][j] * AnsgAu;
    }
    printf("Position of Nuclei %d x=% .8le y=% .8le z=% .8le\n", i + 1,
           center_xyz[i][0], center_xyz[i][1], center_xyz[i][2]);
  }
  fclose(file);

  ns_icent = (int *)calloc(ncent, sizeof(int));
  np_icent = (int *)calloc(ncent, sizeof(int));
  nd_icent = (int *)calloc(ncent, sizeof(int));
  nf_icent = (int *)calloc(ncent, sizeof(int));
  count = 1;
  nbasis = 0;
  file = fopen("gaussian.inp", "r");
  while (count != nlines) {
    for (i = 0; i < ncent; i++) {
      if (count == pos[i]) {
        fscanf(file, "%d%d", &icent, &n);
      }
    }
    fgets(str, 200, file);
    if (strncmp(str, " s  ", 4) == 0) {
      ns_icent[icent - 1] = ns_icent[icent - 1] + 1;
      nbasis = nbasis + 1;
    }
    if (strncmp(str, " p  ", 4) == 0) {
      np_icent[icent - 1] = np_icent[icent - 1] + 1;
      nbasis = nbasis + 3;
    }
    if (strncmp(str, " d  ", 4) == 0) {
      nd_icent[icent - 1] = nd_icent[icent - 1] + 1;
      nbasis = nbasis + 6;
    }
    if (strncmp(str, " f  ", 4) == 0) {
      nf_icent[icent - 1] = nf_icent[icent - 1] + 1;
      nbasis = nbasis + 10;
    }
    count++;
  }
  fclose(file);
  free(pos);

  nprim_s_icent = (int **)calloc(ncent, sizeof(int *));
  nprim_p_icent = (int **)calloc(ncent, sizeof(int *));
  nprim_d_icent = (int **)calloc(ncent, sizeof(int *));
  nprim_f_icent = (int **)calloc(ncent, sizeof(int *));
  for (i = 0; i < ncent; i++) {
    nprim_s_icent[i] = (int *)calloc(ns_icent[i], sizeof(int));
    nprim_p_icent[i] = (int *)calloc(np_icent[i], sizeof(int));
    nprim_d_icent[i] = (int *)calloc(nd_icent[i], sizeof(int));
    nprim_f_icent[i] = (int *)calloc(nf_icent[i], sizeof(int));
  }
  aux_s_icent = (int *)calloc(ncent, sizeof(int));
  aux_p_icent = (int *)calloc(ncent, sizeof(int));
  aux_d_icent = (int *)calloc(ncent, sizeof(int));
  aux_f_icent = (int *)calloc(ncent, sizeof(int));

  file = fopen("gaussian.inp", "r");
  for (i = 0; i < ncent; i++) {
    fscanf(file, "%d%d", &icent, &n);
    for (j = 0; j < ns_icent[icent - 1]; j++) {
      fscanf(file, "%s%d%le", str, &k, &a);
      nprim_s_icent[icent - 1][j] = k;
      aux_s_icent[icent - 1] = aux_s_icent[icent - 1] + k;
      for (k = 0; k < nprim_s_icent[icent - 1][j]; k++) {
        fscanf(file, "%s%s", str, str);
      }
    }
    for (j = 0; j < np_icent[icent - 1]; j++) {
      fscanf(file, "%s%d%le", str, &k, &a);
      nprim_p_icent[icent - 1][j] = k;
      aux_p_icent[icent - 1] = aux_p_icent[icent - 1] + k;
      for (k = 0; k < nprim_p_icent[icent - 1][j]; k++) {
        fscanf(file, "%s%s", str, str);
      }
    }
    for (j = 0; j < nd_icent[icent - 1]; j++) {
      fscanf(file, "%s%d%le", str, &k, &a);
      nprim_d_icent[icent - 1][j] = k;
      aux_d_icent[icent - 1] = aux_d_icent[icent - 1] + k;
      for (k = 0; k < nprim_d_icent[icent - 1][j]; k++) {
        fscanf(file, "%s%s", str, str);
      }
    }
    for (j = 0; j < nf_icent[icent - 1]; j++) {
      fscanf(file, "%s%d%le", str, &k, &a);
      nprim_f_icent[icent - 1][j] = k;
      aux_f_icent[icent - 1] = aux_f_icent[icent - 1] + k;
      for (k = 0; k < nprim_f_icent[icent - 1][j]; k++) {
        fscanf(file, "%s%s", str, str);
      }
    }
  }
  fclose(file);

  coeff_s = (double **)calloc(ncent, sizeof(double));
  coeff_p = (double **)calloc(ncent, sizeof(double));
  coeff_d = (double **)calloc(ncent, sizeof(double));
  coeff_f = (double **)calloc(ncent, sizeof(double));
  expon_s = (double **)calloc(ncent, sizeof(double));
  expon_p = (double **)calloc(ncent, sizeof(double));
  expon_d = (double **)calloc(ncent, sizeof(double));
  expon_f = (double **)calloc(ncent, sizeof(double));
  for (i = 0; i < ncent; i++) {
    coeff_s[i] = (double *)calloc(aux_s_icent[i], sizeof(double));
    coeff_p[i] = (double *)calloc(aux_p_icent[i], sizeof(double));
    coeff_d[i] = (double *)calloc(aux_d_icent[i], sizeof(double));
    coeff_f[i] = (double *)calloc(aux_f_icent[i], sizeof(double));
    expon_s[i] = (double *)calloc(aux_s_icent[i], sizeof(double));
    expon_p[i] = (double *)calloc(aux_p_icent[i], sizeof(double));
    expon_d[i] = (double *)calloc(aux_d_icent[i], sizeof(double));
    expon_f[i] = (double *)calloc(aux_f_icent[i], sizeof(double));
  }
  file = fopen("gaussian.inp", "r");
  for (i = 0; i < ncent; i++) {
    fscanf(file, "%d%d", &icent, &n);
    off = 0;
    for (j = 0; j < ns_icent[icent - 1]; j++) {
      fscanf(file, "%s%d%le", str, &k, &a);
      for (k = 0; k < nprim_s_icent[icent - 1][j]; k++) {
        fscanf(file, "%s%s", str1, str2);
        _CHANGE_STR(str1, &val1);
        _CHANGE_STR(str2, &val2);
        expon_s[icent - 1][off] = val1;
        fact = pow(2.0 * val1 / M_PI, 0.75);
        coeff_s[icent - 1][off] = fact * val2;
        off++;
      }
    }
    off = 0;
    for (j = 0; j < np_icent[icent - 1]; j++) {
      fscanf(file, "%s%d%le", str, &k, &a);
      for (k = 0; k < nprim_p_icent[icent - 1][j]; k++) {
        fscanf(file, "%s%s", str1, str2);
        _CHANGE_STR(str1, &val1);
        _CHANGE_STR(str2, &val2);
        expon_p[icent - 1][off] = val1;
        fact = pow(pow(2.0, 7.0) * pow(val1, 5.0) / pow(M_PI, 3.0), 0.25);
        coeff_p[icent - 1][off] = fact * val2;
        off++;
      }
    }
    off = 0;
    for (j = 0; j < nd_icent[icent - 1]; j++) {
      fscanf(file, "%s%d%le", str, &k, &a);
      for (k = 0; k < nprim_d_icent[icent - 1][j]; k++) {
        fscanf(file, "%s%s", str1, str2);
        _CHANGE_STR(str1, &val1);
        _CHANGE_STR(str2, &val2);
        expon_d[icent - 1][off] = val1;
        fact = pow(pow(2.0, 11.0) * pow(val1, 7.0) / pow(M_PI, 3.0), 0.25);
        coeff_d[icent - 1][off] = fact * val2;
        off++;
      }
    }
    off = 0;
    for (j = 0; j < nf_icent[icent - 1]; j++) {
      fscanf(file, "%s%d%le", str, &k, &a);
      for (k = 0; k < nprim_f_icent[icent - 1][j]; k++) {
        fscanf(file, "%s%s", str1, str2);
        _CHANGE_STR(str1, &val1);
        _CHANGE_STR(str2, &val2);
        expon_f[icent - 1][off] = val1;
        fact = pow(pow(2.0, 15.0) * pow(val1, 9.0) / pow(M_PI, 3.0), 0.25);
        coeff_f[icent - 1][off] = fact * val2;
        off++;
      }
    }
  }
  fclose(file);
  free(aux_s_icent);
  free(aux_p_icent);
  free(aux_d_icent);
  free(aux_f_icent);

  printf("Number of orbitals detected=%d\n", n_orb);
  printf("Number of occupied orbitals detected=%d\n", n_occu_orb);
  occu = (double *)calloc(n_occu_orb, sizeof(double));
  mo_occu = (double **)calloc(n_occu_orb, sizeof(double));
  for (i = 0; i < n_occu_orb; i++) {
    mo_occu[i] = (double *)calloc(nbasis, sizeof(double));
  }
  count = 0;
  file = fopen("molecule.inp", "r");
  for (i = 0; i < n_orb; i++) {
    fscanf(file, "%s%s", str, str);
    fscanf(file, "%s%s", str, str);
    fscanf(file, "%s%s", str, str);
    fscanf(file, "%s%le", str, &a);
    if (a != 0.0) {
      occu[count] = a;
      for (j = 0; j < nbasis; j++) {
        fscanf(file, "%d%s", &k, str);
        _CHANGE_STR(str, &val1);
        mo_occu[count][j] = val1;
      }
      count++;
    } else {
      for (j = 0; j < nbasis; j++) {
        fscanf(file, "%d%s", &k, str);
      }
    }
  }
  fclose(file);

  we = (double *)calloc(4, sizeof(double));
  for (i = 0; i < 4; i++) {
    we[i] = 1.0;
    for (j = 2 * i - 1; j > 0; j = j - 2) {
      we[i] = we[i] * sqrt(1.0 / j);
    }
  }
  for (xgrid = 0; xgrid < ngrid; xgrid++) {
    val = 0.0;
    v = 0;
    for (i = 0; i < ncent; i++) {
      x = xi[xgrid] - center_xyz[i][0];
      y = yi[xgrid] - center_xyz[i][1];
      z = zi[xgrid] - center_xyz[i][2];
      off = 0;
      for (j = 0; j < ns_icent[i]; j++) {
        S_xyz = 0.0;
        for (k = 0; k < nprim_s_icent[i][j]; k++) {
          val1 = expon_s[i][off];
          val2 = coeff_s[i][off];
          S_xyz = S_xyz + val2 * exp(-val1 * (x * x + y * y + z * z));
          off++;
        }
        val = val + S_xyz * mo_occu[l][v];
        v++;
      }
      off = 0;
      for (j = 0; j < np_icent[i]; j++) {
        S_xyz = 0.0;
        for (k = 0; k < nprim_p_icent[i][j]; k++) {
          val1 = expon_p[i][off];
          val2 = coeff_p[i][off];
          S_xyz = S_xyz + val2 * exp(-val1 * (x * x + y * y + z * z));
          off++;
        }
        val = val + S_xyz * x * mo_occu[l][v];
        v++;
        val = val + S_xyz * y * mo_occu[l][v];
        v++;
        val = val + S_xyz * z * mo_occu[l][v];
        v++;
      }
      off = 0;
      for (j = 0; j < nd_icent[i]; j++) {
        S_xyz = 0.0;
        for (k = 0; k < nprim_d_icent[i][j]; k++) {
          val1 = expon_d[i][off];
          val2 = coeff_d[i][off];
          S_xyz = S_xyz + val2 * exp(-val1 * (x * x + y * y + z * z));
          off++;
        }
        val = val + S_xyz * pow(x, 2) * mo_occu[l][v] * we[2];
        v++;
        val = val + S_xyz * pow(y, 2) * mo_occu[l][v] * we[2];
        v++;
        val = val + S_xyz * pow(z, 2) * mo_occu[l][v] * we[2];
        v++;
        val = val + S_xyz * x * y * mo_occu[l][v];
        v++;
        val = val + S_xyz * x * z * mo_occu[l][v];
        v++;
        val = val + S_xyz * y * z * mo_occu[l][v];
        v++;
      }
      off = 0;
      for (j = 0; j < nf_icent[i]; j++) {
        S_xyz = 0.0;
        for (k = 0; k < nprim_f_icent[i][j]; k++) {
          val1 = expon_f[i][off];
          val2 = coeff_f[i][off];
          S_xyz = S_xyz + val2 * exp(-val1 * (x * x + y * y + z * z));
          off++;
        }
        val = val + S_xyz * pow(x, 3) * mo_occu[l][v] * we[3];
        v++;
        val = val + S_xyz * pow(y, 3) * mo_occu[l][v] * we[3];
        v++;
        val = val + S_xyz * pow(z, 3) * mo_occu[l][v] * we[3];
        v++;
        val = val + S_xyz * pow(y, 2) * x * mo_occu[l][v] * we[2];
        v++;
        val = val + S_xyz * pow(x, 2) * y * mo_occu[l][v] * we[2];
        v++;
        val = val + S_xyz * pow(x, 2) * z * mo_occu[l][v] * we[2];
        v++;
        val = val + S_xyz * pow(z, 2) * x * mo_occu[l][v] * we[2];
        v++;
        val = val + S_xyz * pow(z, 2) * y * mo_occu[l][v] * we[2];
        v++;
        val = val + S_xyz * pow(y, 2) * z * mo_occu[l][v] * we[2];
        v++;
        val = val + S_xyz * x * y * z * mo_occu[l][v];
        v++;
      }
    }
    mo_orb[xgrid] = val;
  }

  free(we);
  free(ns_icent);
  free(np_icent);
  free(nd_icent);
  free(nf_icent);
  free(nprim_s_icent);
  free(nprim_p_icent);
  free(nprim_d_icent);
  free(nprim_f_icent);
  free(expon_s);
  free(coeff_s);
  free(expon_p);
  free(coeff_p);
  free(expon_d);
  free(coeff_d);
  free(expon_f);
  free(coeff_f);
  free(occu);
  free(mo_occu);
  free(center_xyz);
}

void _CHANGE_STR(char *str, double *val) {  // ok
  int i, len;
  char *end;
  len = strlen(str);
  for (i = 0; i < len; i++) {
    if (str[i] == 'D') {
      str[i] = 'E';
    }
  }
  *val = strtod(str, &end);  //*val=atof(str);
}
