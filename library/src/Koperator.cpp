#include "../include/Koperator.hpp"
#include <gsl/gsl_sf_legendre.h>
#include <iostream>
#include <mkl_lapacke.h>
#include <mpi.h>

Koperator::Koperator(
    const size_t &Nbas, std::shared_ptr<AngularGrid> a_angular_grid,
    std::shared_ptr<MOPartialWaveRepresentation> a_ket_bra_orbital,
    std::shared_ptr<MOPartialWaveRepresentation> a_J_orbital,
    std::shared_ptr<Toperator> a_T) {

  size_t grid_Lmax = a_angular_grid->getLmax();
  size_t num_spherical_harmonics = (grid_Lmax - 1) * (grid_Lmax + 1) + 1;

  /* Building the angular coefficients */
  size_t ket_bra_num_channels = a_ket_bra_orbital->getNumChannels();
  size_t J_orbital_num_channels = a_J_orbital->getNumChannels();
  auto CJ_1 = std::make_unique<std::complex<double>[]>(
      ket_bra_num_channels * J_orbital_num_channels * num_spherical_harmonics);
  auto CJ_2 = std::make_unique<std::complex<double>[]>(
      ket_bra_num_channels * J_orbital_num_channels * num_spherical_harmonics);
  for (int i = 0; i < ket_bra_num_channels; ++i) {
    for (int j = 0; j < J_orbital_num_channels; ++j) {
      int k = 0;
      for (int lv = 0; lv < grid_Lmax; ++lv) {
        if (a_angular_grid->getL(k) >=
                abs(a_ket_bra_orbital->getL(i) - a_J_orbital->getL(j)) &&
            a_angular_grid->getL(k) <=
                a_ket_bra_orbital->getL(i) + a_J_orbital->getL(j)) {
          for (int mv = 0; mv <= lv; ++mv) {
            std::complex<double> zeta(0, 0), beta(0, 0);
            if (mv == 0) {
              C3jBlm(a_angular_grid, a_angular_grid->getL(k),
                     a_angular_grid->getM(k), a_ket_bra_orbital->getType(),
                     a_ket_bra_orbital->getL(i), a_ket_bra_orbital->getM(i),
                     a_J_orbital->getType(), a_J_orbital->getL(j),
                     a_J_orbital->getM(j), zeta, beta);
              CJ_1[i * ket_bra_num_channels * num_spherical_harmonics +
                   j * num_spherical_harmonics + k] = zeta;
              CJ_2[i * J_orbital_num_channels * num_spherical_harmonics +
                   j * num_spherical_harmonics + k] = beta;
              k++;
            } else {
              C3jBlm(a_angular_grid, a_angular_grid->getL(k),
                     a_angular_grid->getM(k), a_ket_bra_orbital->getType(),
                     a_ket_bra_orbital->getL(i), a_ket_bra_orbital->getM(i),
                     a_J_orbital->getType(), a_J_orbital->getL(j),
                     a_J_orbital->getM(j), zeta, beta);
              CJ_1[i * ket_bra_num_channels * num_spherical_harmonics +
                   j * num_spherical_harmonics + k] = zeta;
              CJ_2[i * J_orbital_num_channels * num_spherical_harmonics +
                   j * num_spherical_harmonics + k] = beta;
              k++;
              CJ_1[i * ket_bra_num_channels * num_spherical_harmonics +
                   j * num_spherical_harmonics + k] = beta;
              CJ_2[i * J_orbital_num_channels * num_spherical_harmonics +
                   j * num_spherical_harmonics + k] = zeta;
              k++;
            }
          }
        } else {
          for (int mv = 0; mv <= lv; ++mv) {
            if (mv == 0) {
              k++;
            } else {
              k = k + 2;
            }
          }
        }
      }
    }
  }

  m_dvr_rep = std::make_unique<std::complex<double>[]>(
      Nbas * Nbas * ket_bra_num_channels * ket_bra_num_channels);
  for (int i = 0; i < Nbas * ket_bra_num_channels; ++i) {
    int ket_r = i / Nbas;
    int ket_l = i % Nbas;
    for (int j = 0; j < Nbas * ket_bra_num_channels; ++j) {
      int bra_r = j / Nbas;
      int bra_l = j % Nbas;
      std::complex<double> tmp_J(1, 0);
      for (int L = 0; L < num_spherical_harmonics; ++L) {
        std::complex<double> tmp_1 =
            4.0 * M_PI * pow(-1, a_angular_grid->getM(L)) *
            a_T->getTIXX(a_angular_grid->getL(L) * Nbas * Nbas + bra_r * Nbas +
                         ket_r) /
            (2.0 * a_angular_grid->getL(L) + 1.0);
        for (int l1 = 0; l1 < J_orbital_num_channels; ++l1) {
          if (a_angular_grid->getL(L) >=
                  abs(a_ket_bra_orbital->getL(ket_l) - a_J_orbital->getL(l1)) &&
              a_angular_grid->getL(L) <=
                  a_ket_bra_orbital->getL(ket_l) + a_J_orbital->getL(l1)) {
            std::complex<double> tmp_2 =
                CJ_1[ket_l * J_orbital_num_channels * num_spherical_harmonics +
                     l1 * num_spherical_harmonics + L];
            for (int l2 = 0; l2 < J_orbital_num_channels; ++l2) {
              if (a_angular_grid->getL(L) >=
                      abs(a_ket_bra_orbital->getL(bra_l) -
                          a_J_orbital->getL(l2)) &&
                  a_angular_grid->getL(L) <=
                      a_ket_bra_orbital->getL(bra_l) + a_J_orbital->getL(l2)) {
                tmp_J =
                    tmp_1 * tmp_2 *
                    CJ_2[l1 * J_orbital_num_channels * num_spherical_harmonics +
                         l2 * num_spherical_harmonics + L];
                a_J_orbital->getPartialWaveRep(l1 * Nbas + ket_r) *
                    a_J_orbital->getPartialWaveRep(l2 * Nbas + bra_r);
              }
            }
          }
        }
      }
      m_dvr_rep[i * Nbas * ket_bra_num_channels + j] = tmp_J;
    }
  }
}

Koperator::Koperator(
    const int &a_numprocs, const std::array<int, 2> &a_proccessor_xy,
    const size_t &Nbas, std::shared_ptr<AngularGrid> a_angular_grid,
    std::shared_ptr<MOPartialWaveRepresentation> a_ket_bra_orbital,
    std::shared_ptr<MOPartialWaveRepresentation> a_J_orbital,
    std::shared_ptr<Toperator> a_T) {

  size_t grid_Lmax = a_angular_grid->getLmax();
  size_t num_spherical_harmonics = (grid_Lmax - 1) * (grid_Lmax + 1) + 1;

  /* Building the angular coefficients */
  size_t ket_bra_num_channels = a_ket_bra_orbital->getNumChannels();
  size_t J_orbital_num_channels = a_J_orbital->getNumChannels();
  auto CJ_1 = std::make_unique<std::complex<double>[]>(
      ket_bra_num_channels * J_orbital_num_channels * num_spherical_harmonics);
  auto CJ_2 = std::make_unique<std::complex<double>[]>(
      ket_bra_num_channels * J_orbital_num_channels * num_spherical_harmonics);
  for (int i = 0; i < ket_bra_num_channels; ++i) {
    for (int j = 0; j < J_orbital_num_channels; ++j) {
      int k = 0;
      for (int lv = 0; lv < grid_Lmax; ++lv) {
        if (a_angular_grid->getL(k) >=
                abs(a_ket_bra_orbital->getL(i) - a_J_orbital->getL(j)) &&
            a_angular_grid->getL(k) <=
                a_ket_bra_orbital->getL(i) + a_J_orbital->getL(j)) {
          for (int mv = 0; mv <= lv; ++mv) {
            std::complex<double> zeta(0, 0), beta(0, 0);
            if (mv == 0) {
              C3jBlm(a_angular_grid, a_angular_grid->getL(k),
                     a_angular_grid->getM(k), a_ket_bra_orbital->getType(),
                     a_ket_bra_orbital->getL(i), a_ket_bra_orbital->getM(i),
                     a_J_orbital->getType(), a_J_orbital->getL(j),
                     a_J_orbital->getM(j), zeta, beta);
              CJ_1[i * ket_bra_num_channels * num_spherical_harmonics +
                   j * num_spherical_harmonics + k] = zeta;
              CJ_2[i * J_orbital_num_channels * num_spherical_harmonics +
                   j * num_spherical_harmonics + k] = beta;
              k++;
            } else {
              C3jBlm(a_angular_grid, a_angular_grid->getL(k),
                     a_angular_grid->getM(k), a_ket_bra_orbital->getType(),
                     a_ket_bra_orbital->getL(i), a_ket_bra_orbital->getM(i),
                     a_J_orbital->getType(), a_J_orbital->getL(j),
                     a_J_orbital->getM(j), zeta, beta);
              CJ_1[i * ket_bra_num_channels * num_spherical_harmonics +
                   j * num_spherical_harmonics + k] = zeta;
              CJ_2[i * J_orbital_num_channels * num_spherical_harmonics +
                   j * num_spherical_harmonics + k] = beta;
              k++;
              CJ_1[i * ket_bra_num_channels * num_spherical_harmonics +
                   j * num_spherical_harmonics + k] = beta;
              CJ_2[i * J_orbital_num_channels * num_spherical_harmonics +
                   j * num_spherical_harmonics + k] = zeta;
              k++;
            }
          }
        } else {
          for (int mv = 0; mv <= lv; ++mv) {
            if (mv == 0) {
              k++;
            } else {
              k = k + 2;
            }
          }
        }
      }
    }
  }

  for (int i = 0; i < Nbas * ket_bra_num_channels; ++i) {
    int ket_r = i / Nbas;
    int ket_l = i % Nbas;
    for (int j = 0; j < Nbas * ket_bra_num_channels; ++j) {
      int bra_r = j / Nbas;
      int bra_l = j % Nbas;
      std::complex<double> tmp_J(1, 0);
      for (int L = 0; L < num_spherical_harmonics; ++L) {
        std::complex<double> tmp_1 =
            4.0 * M_PI * pow(-1, a_angular_grid->getM(L)) *
            a_T->getTIXX(a_angular_grid->getL(L) * Nbas * Nbas + bra_r * Nbas +
                         ket_r) /
            (2.0 * a_angular_grid->getL(L) + 1.0);
        for (int l1 = 0; l1 < J_orbital_num_channels; ++l1) {
          if (a_angular_grid->getL(L) >=
                  abs(a_ket_bra_orbital->getL(ket_l) - a_J_orbital->getL(l1)) &&
              a_angular_grid->getL(L) <=
                  a_ket_bra_orbital->getL(ket_l) + a_J_orbital->getL(l1)) {
            std::complex<double> tmp_2 =
                CJ_1[ket_l * ket_bra_num_channels * num_spherical_harmonics +
                     bra_l * num_spherical_harmonics + L];
            for (int l2 = 0; l2 < J_orbital_num_channels; ++l2) {
              if (a_angular_grid->getL(L) >=
                      abs(a_ket_bra_orbital->getL(bra_l) -
                          a_J_orbital->getL(l2)) &&
                  a_angular_grid->getL(L) <=
                      a_ket_bra_orbital->getL(bra_l) + a_J_orbital->getL(l2)) {
                tmp_J =
                    tmp_1 * tmp_2 *
                    CJ_2[l1 * J_orbital_num_channels * num_spherical_harmonics +
                         l2 * num_spherical_harmonics + L];
                a_ket_bra_orbital->getPartialWaveRep(l1 * Nbas + ket_r) *
                    a_ket_bra_orbital->getPartialWaveRep(l2 * Nbas + bra_r);
              }
            }
          }
        }
      }
      m_dvr_rep[i * Nbas * ket_bra_num_channels + j] = tmp_J;
    }
  }
}

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
