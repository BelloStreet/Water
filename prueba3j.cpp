#include <complex.h>
#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_sf_legendre.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "sphere_lebedev_rule.h"
typedef _Complex double cmplx;
cmplx _Ylm(int l, int m, double theta, double phi);
double _Blm(char *type, int l, int m, double theta, double phi);
int main(void) {
  int i, j, k, l, j1 = 2, j2 = 1, j3 = 1, m1 = 1, m2 = 0, m3 = -1;
  double y, z, *xi, *yi, *zi, *wi, phi, theta;
  cmplx beta, z1, z2, z3, z4;
  char str1[3] = "s", str2[3] = "c";
  const int n = 974;
  xi = (double *)calloc(n, sizeof(double));
  yi = (double *)calloc(n, sizeof(double));
  zi = (double *)calloc(n, sizeof(double));
  wi = (double *)calloc(n, sizeof(double));
  ld0974(xi, yi, zi, wi);
  beta = 0.0;
  for (j = 0; j < n; j++) {
    theta = acos(zi[j]);
    phi = atan2(yi[j], xi[j]);
    //    beta=beta+wi[j]*_Ylm(j1,m1,theta,phi)*_Ylm(j2,m2,theta,phi)*_Ylm(j3,m3,theta,phi);
    beta = beta + wi[j] * _Blm(str1, j1, m1, theta, phi) *
                      _Blm(str2, j2, m2, theta, phi) * _Ylm(j3, m3, theta, phi);
  }
  beta = 4.0 * M_PI * beta;
  //  y=gsl_sf_coupling_3j(2*j1,2*j2,2*j3,2*m1,2*m2,2*m3)*gsl_sf_coupling_3j(2*j1,2*j2,2*j3,0,0,0);
  //  y=y*sqrt((2.0*j1+1.0)*(2.0*j2+1.0)*(2.0*j3+1.0)/(4.0*M_PI));

  z1 = gsl_sf_coupling_3j(2 * j1, 2 * j2, 2 * j3, 2 * m1, 2 * m2, 2 * m3) *
       gsl_sf_coupling_3j(2 * j1, 2 * j2, 2 * j3, 0, 0, 0);
  z1 = z1 * sqrt((2.0 * j1 + 1.0) * (2.0 * j2 + 1.0) * (2.0 * j3 + 1.0) /
                 (4.0 * M_PI));
  z2 = gsl_sf_coupling_3j(2 * j1, 2 * j2, 2 * j3, 2 * m1, -2 * m2, 2 * m3) *
       gsl_sf_coupling_3j(2 * j1, 2 * j2, 2 * j3, 0, 0, 0);
  z2 = z2 * sqrt((2.0 * j1 + 1.0) * (2.0 * j2 + 1.0) * (2.0 * j3 + 1.0) /
                 (4.0 * M_PI));
  z3 = gsl_sf_coupling_3j(2 * j1, 2 * j2, 2 * j3, -2 * m1, 2 * m2, 2 * m3) *
       gsl_sf_coupling_3j(2 * j1, 2 * j2, 2 * j3, 0, 0, 0);
  z3 = z3 * sqrt((2.0 * j1 + 1.0) * (2.0 * j2 + 1.0) * (2.0 * j3 + 1.0) /
                 (4.0 * M_PI));
  z4 = gsl_sf_coupling_3j(2 * j1, 2 * j2, 2 * j3, -2 * m1, -2 * m2, 2 * m3) *
       gsl_sf_coupling_3j(2 * j1, 2 * j2, 2 * j3, 0, 0, 0);
  z4 = z4 * sqrt((2.0 * j1 + 1.0) * (2.0 * j2 + 1.0) * (2.0 * j3 + 1.0) /
                 (4.0 * M_PI));

  y = 0.5 * creal(z1 + z2 - z3 - z4);

  printf("% .8le % .8le % .8le % .8le\n", creal(beta), cimag(beta),
         0.5 * creal(z1 + z2 - z3 - z4), 0.5 * cimag(z1 + z2 - z3 - z4));

  free(xi);
  free(yi);
  free(zi);
  free(wi);
  return 0;
}

cmplx _Ylm(int l, int m, double theta, double phi) {
  int m1;
  double y;
  cmplx c1, c2, c3;
  m1 = abs(m);
  y = gsl_sf_legendre_sphPlm(l, m1, cos(theta));
  c1 = cos(m1 * phi) + sin(m1 * phi) * I;
  if (m >= 0) {
    c3 = c1 * y;
  } else {
    c3 = pow(-1.0, m) * conj(c1 * y);
  }
  return c3;
}

double _Blm(char *type, int l, int m, double theta, double phi) {
  int m1;
  double y, c1;
  m1 = abs(m);
  if (m != 0) {
    y = sqrt(2.0) * gsl_sf_legendre_sphPlm(l, m1, cos(theta));
  } else {
    y = gsl_sf_legendre_sphPlm(l, m1, cos(theta));
  }
  if (strcmp(type, "c") == 0) {
    c1 = y * cos(m1 * phi);
  }
  if (strcmp(type, "s") == 0) {
    c1 = y * sin(m1 * phi);
  }
  return c1;
}
