#include "../include/Koperator.hpp"
#include <gsl/gsl_sf_legendre.h>
#include <iostream>
#include <mkl_lapacke.h>
#include <mpi.h>

Koperator::Koperator(std::shared_ptr<AngularGrid> a_angular_grid,
                     std::shared_ptr<MOPartialWaveRepresentation> a_ket_orbital,
                     std::shared_ptr<MOPartialWaveRepresentation> a_bra_orbital,
                     std::shared_ptr<Toperator> a_T) {}

Koperator::~Koperator() {}

std::complex<double> Koperator::getK(int index) const {
  return m_dvr_rep[index];
}

void Koperator::C3jBlm(std::shared_ptr<AngularGrid> a_angular_grid,
                       const int &L, const int &M, const std::string &type1,
                       const int &l1, const int &m1, const std::string &type2,
                       const int &l2, const int &m2,
                       std::complex<double> &a_zeta,
                       std::complex<double> &a_beta) {
  for (int i = 0; i < a_angular_grid->getAngularOrder(); i++) {
    double theta = acos(a_angular_grid->getZ(i));
    double phi = atan2(a_angular_grid->getY(i), a_angular_grid->getX(i));
    a_zeta += a_angular_grid->getAngularWeight(i) *
              Blm(type1, l1, m1, theta, phi) * Blm(type2, l2, m2, theta, phi) *
              Ylm(L, M, theta, phi);
    a_beta += a_angular_grid->getAngularWeight(i) *
              Blm(type1, l1, m1, theta, phi) * Blm(type2, l2, m2, theta, phi) *
              Ylm(L, -M, theta, phi);
  }
}

std::complex<double> Koperator::Ylm(const int l, const int m,
                                    const double theta, const double phi) {
  int m1 = abs(m);
  double y = gsl_sf_legendre_sphPlm(l, m1, cos(theta));
  std::complex<double> c1 = std::polar(1.0, m1 * phi);
  std::complex<double> c3;
  if (m >= 0) {
    c3 = c1 * y;
  } else {
    c3 = pow(-1.0, m) * conj(c1 * y);
  }
  return c3;
}

double Koperator::Blm(const std::string type, const int l, const int m,
                      const double theta, const double phi) {
  double value = 0.0, y = 0.0;
  int m1 = abs(m);
  if (m != 0) {
    y = sqrt(2.0) * gsl_sf_legendre_sphPlm(l, m1, cos(theta));
  } else {
    y = gsl_sf_legendre_sphPlm(l, m1, cos(theta));
  }
  if (type == "cosine") {
    value = y * cos(m1 * phi);
  } else if (type == "sine") {
    value = y * sin(m1 * phi);
  }
  return value;
}
