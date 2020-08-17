#include "../include/Nuclearoperator.hpp"
#include <gsl/gsl_sf_legendre.h>
#include <mkl_lapacke.h>
#include <mpi.h>

Nuclearoperator::Nuclearoperator(
    const size_t &a_Nbas, std::shared_ptr<AngularGrid> a_angular_grid,
    std::shared_ptr<MOPartialWaveRepresentation> a_ket_bra_orbital,
    std::shared_ptr<MOPartialWaveRepresentation> a_J_orbital,
    std::shared_ptr<Toperator> a_T) {

  size_t grid_Lmax = a_angular_grid->getLmax();
  size_t num_spherical_harmonics = (grid_Lmax - 1) * (grid_Lmax + 1) + 1;
  int test_smallerLmax = floor(grid_Lmax / 2);

  /* Building the angular coefficients */
  size_t ket_bra_num_channels = a_ket_bra_orbital->getNumChannels();
  auto CJ_1 = std::make_unique<std::complex<double>[]>(
      ket_bra_num_channels * ket_bra_num_channels * num_spherical_harmonics);

  for (int i = 0; i < ket_bra_num_channels; ++i) {
    for (int j = 0; j < ket_bra_num_channels; ++j) {
      int k = 0;
      for (int lv = 0; lv < grid_Lmax; ++lv) {
        if (a_angular_grid->getL(k) >=
                abs(a_ket_bra_orbital->getL(i) - a_ket_bra_orbital->getL(j)) &&
            a_angular_grid->getL(k) <=
                a_ket_bra_orbital->getL(i) + a_ket_bra_orbital->getL(j)) {
          for (int mv = 0; mv <= lv; ++mv) {
            std::complex<double> zeta(0, 0), beta(0, 0);
            if (mv == 0) {
              C3jBlm(a_angular_grid, a_angular_grid->getL(k),
                     a_angular_grid->getM(k), a_ket_bra_orbital->getType(),
                     a_ket_bra_orbital->getL(i), a_ket_bra_orbital->getM(i),
                     a_ket_bra_orbital->getType(), a_ket_bra_orbital->getL(j),
                     a_ket_bra_orbital->getM(j), zeta, beta);
              CJ_1[i * ket_bra_num_channels * num_spherical_harmonics +
                   j * num_spherical_harmonics + k] = zeta;
              CJ_1[j * ket_bra_num_channels * num_spherical_harmonics +
                   i * num_spherical_harmonics + k] = zeta;
              k++;
            } else {
              C3jBlm(a_angular_grid, a_angular_grid->getL(k),
                     a_angular_grid->getM(k), a_ket_bra_orbital->getType(),
                     a_ket_bra_orbital->getL(i), a_ket_bra_orbital->getM(i),
                     a_ket_bra_orbital->getType(), a_ket_bra_orbital->getL(j),
                     a_ket_bra_orbital->getM(j), zeta, beta);
              CJ_1[i * ket_bra_num_channels * num_spherical_harmonics +
                   j * num_spherical_harmonics + k] = zeta;
              CJ_1[j * ket_bra_num_channels * num_spherical_harmonics +
                   i * num_spherical_harmonics + k] = zeta;
              k++;
              CJ_1[i * ket_bra_num_channels * num_spherical_harmonics +
                   j * num_spherical_harmonics + k] = beta;
              CJ_1[j * ket_bra_num_channels * num_spherical_harmonics +
                   i * num_spherical_harmonics + k] = beta;
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

  size_t J_orbital_num_channels = a_J_orbital->getNumChannels();
  auto CJ_2 = std::make_unique<std::complex<double>[]>(J_orbital_num_channels *
                                                       J_orbital_num_channels *
                                                       num_spherical_harmonics);

  for (int i = 0; i < J_orbital_num_channels; ++i) {
    for (int j = 0; j < J_orbital_num_channels; ++j) {
      int k = 0;
      for (int lv = 0; lv < grid_Lmax; ++lv) {
        if (a_angular_grid->getL(k) >=
                abs(a_J_orbital->getL(i) - a_J_orbital->getL(j)) &&
            a_angular_grid->getL(k) <=
                a_J_orbital->getL(i) + a_J_orbital->getL(j)) {
          for (int mv = 0; mv <= lv; ++mv) {
            std::complex<double> zeta(0, 0), beta(0, 0);
            if (mv == 0) {
              C3jBlm(a_angular_grid, a_angular_grid->getL(k),
                     a_angular_grid->getM(k), a_J_orbital->getType(),
                     a_J_orbital->getL(i), a_J_orbital->getM(i),
                     a_J_orbital->getType(), a_J_orbital->getL(j),
                     a_J_orbital->getM(j), zeta, beta);
              CJ_2[i * J_orbital_num_channels * num_spherical_harmonics +
                   j * num_spherical_harmonics + k] = zeta;
              CJ_2[j * J_orbital_num_channels * num_spherical_harmonics +
                   i * num_spherical_harmonics + k] = zeta;
              k++;
            } else {
              C3jBlm(a_angular_grid, a_angular_grid->getL(k),
                     a_angular_grid->getM(k), a_J_orbital->getType(),
                     a_J_orbital->getL(i), a_J_orbital->getM(i),
                     a_J_orbital->getType(), a_J_orbital->getL(j),
                     a_J_orbital->getM(j), zeta, beta);
              CJ_2[i * J_orbital_num_channels * num_spherical_harmonics +
                   j * num_spherical_harmonics + k] = zeta;
              CJ_2[j * J_orbital_num_channels * num_spherical_harmonics +
                   i * num_spherical_harmonics + k] = zeta;
              k++;
              CJ_2[i * J_orbital_num_channels * num_spherical_harmonics +
                   j * num_spherical_harmonics + k] = beta;
              CJ_2[j * J_orbital_num_channels * num_spherical_harmonics +
                   i * num_spherical_harmonics + k] = beta;
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
      a_Nbas * a_Nbas * ket_bra_num_channels * ket_bra_num_channels);
  for (int i = 0; i < a_Nbas * ket_bra_num_channels; ++i) {
    int ket_r = i / ket_bra_num_channels;
    int ket_l = i % ket_bra_num_channels;
    for (int j = 0; j < a_Nbas * ket_bra_num_channels; ++j) {
      int bra_r = j / ket_bra_num_channels;
      int bra_l = j % ket_bra_num_channels;
      std::complex<double> tmp_J(0, 0);
      if (ket_r == bra_r) {
        for (int lk = 0; lk < a_Nbas * num_spherical_harmonics; ++lk) {
          int k = lk / num_spherical_harmonics;
          int L = lk % num_spherical_harmonics;
          if (a_angular_grid->getL(L) >= abs(a_ket_bra_orbital->getL(ket_l) -
                                             a_ket_bra_orbital->getL(bra_l)) &&
              a_angular_grid->getL(L) <= a_ket_bra_orbital->getL(ket_l) +
                                             a_ket_bra_orbital->getL(bra_l)) {
            std::complex<double> tmp_1 =
                4.0 * M_PI * pow(-1, a_angular_grid->getM(L)) *
                a_T->getTIXX(a_angular_grid->getL(L) * a_Nbas * a_Nbas +
                             k * a_Nbas + ket_r) /
                (2.0 * a_angular_grid->getL(L) + 1.0);
            printf("L %d k %d ket_r %d \n", L, k, ket_r);
            /* printf(" z2 %f  \n", */
            /*        a_T->getTIXX(a_angular_grid->getL(L) * a_Nbas * a_Nbas +
             */
            /* k * a_Nbas + ket_r)); */
            std::complex<double> tmp_2 =
                CJ_1[ket_l * ket_bra_num_channels * num_spherical_harmonics +
                     bra_l * num_spherical_harmonics + L];
            for (int l1 = 0; l1 < J_orbital_num_channels; ++l1) {
              for (int l2 = 0; l2 < J_orbital_num_channels; ++l2) {
                if (a_angular_grid->getL(L) >=
                        abs(a_J_orbital->getL(l1) - a_J_orbital->getL(l2)) &&
                    a_angular_grid->getL(L) <=
                        a_J_orbital->getL(l1) + a_J_orbital->getL(l2)) {
                  tmp_J += tmp_1 * tmp_2 *
                           CJ_2[l1 * J_orbital_num_channels *
                                    num_spherical_harmonics +
                                l2 * num_spherical_harmonics + L];
                  a_J_orbital->getPartialWaveRep(l1 * a_Nbas + k) *
                      a_J_orbital->getPartialWaveRep(l2 * a_Nbas + k);
                  /* printf(" z2 %f  z3 %f\n", tmp_1,
                   * tmp_2); */
                }
              }
            }
          }
        }
      }
      /* printf("tmp_J %f \n", tmp_J); */
      m_dvr_rep[i * a_Nbas * ket_bra_num_channels + j] = tmp_J;
    }
  }
}

Nuclearoperator::Nuclearoperator(
    const int &a_numprocs, const int &id, const size_t &a_Nbas,
    std::shared_ptr<AngularGrid> a_angular_grid,
    std::shared_ptr<MOPartialWaveRepresentation> a_ket_bra_orbital,
    std::shared_ptr<MOPartialWaveRepresentation> a_J_orbital,
    std::shared_ptr<Toperator> a_T) {

  size_t grid_Lmax = a_angular_grid->getLmax();
  size_t num_spherical_harmonics = (grid_Lmax - 1) * (grid_Lmax + 1) + 1;
  int test_smallerLmax = floor(grid_Lmax / 2);

  /* Building the angular coefficients */
  size_t ket_bra_num_channels = a_ket_bra_orbital->getNumChannels();
  auto CJ_1 = std::make_unique<std::complex<double>[]>(
      ket_bra_num_channels * ket_bra_num_channels * num_spherical_harmonics);

  for (int i = 0; i < ket_bra_num_channels; ++i) {
    for (int j = 0; j < ket_bra_num_channels; ++j) {
      int k = 0;
      for (int lv = 0; lv < grid_Lmax; ++lv) {
        if (a_angular_grid->getL(k) >=
                abs(a_ket_bra_orbital->getL(i) - a_ket_bra_orbital->getL(j)) &&
            a_angular_grid->getL(k) <=
                a_ket_bra_orbital->getL(i) + a_ket_bra_orbital->getL(j)) {
          for (int mv = 0; mv <= lv; ++mv) {
            std::complex<double> zeta(0, 0), beta(0, 0);
            if (mv == 0) {
              C3jBlm(a_angular_grid, a_angular_grid->getL(k),
                     a_angular_grid->getM(k), a_ket_bra_orbital->getType(),
                     a_ket_bra_orbital->getL(i), a_ket_bra_orbital->getM(i),
                     a_ket_bra_orbital->getType(), a_ket_bra_orbital->getL(j),
                     a_ket_bra_orbital->getM(j), zeta, beta);
              CJ_1[i * ket_bra_num_channels * num_spherical_harmonics +
                   j * num_spherical_harmonics + k] = zeta;
              CJ_1[j * ket_bra_num_channels * num_spherical_harmonics +
                   i * num_spherical_harmonics + k] = zeta;
              k++;
            } else {
              C3jBlm(a_angular_grid, a_angular_grid->getL(k),
                     a_angular_grid->getM(k), a_ket_bra_orbital->getType(),
                     a_ket_bra_orbital->getL(i), a_ket_bra_orbital->getM(i),
                     a_ket_bra_orbital->getType(), a_ket_bra_orbital->getL(j),
                     a_ket_bra_orbital->getM(j), zeta, beta);
              CJ_1[i * ket_bra_num_channels * num_spherical_harmonics +
                   j * num_spherical_harmonics + k] = zeta;
              CJ_1[j * ket_bra_num_channels * num_spherical_harmonics +
                   i * num_spherical_harmonics + k] = zeta;
              k++;
              CJ_1[i * ket_bra_num_channels * num_spherical_harmonics +
                   j * num_spherical_harmonics + k] = beta;
              CJ_1[j * ket_bra_num_channels * num_spherical_harmonics +
                   i * num_spherical_harmonics + k] = beta;
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

  size_t J_orbital_num_channels = a_J_orbital->getNumChannels();
  auto CJ_2 = std::make_unique<std::complex<double>[]>(J_orbital_num_channels *
                                                       J_orbital_num_channels *
                                                       num_spherical_harmonics);

  for (int i = 0; i < J_orbital_num_channels; ++i) {
    for (int j = 0; j < J_orbital_num_channels; ++j) {
      int k = 0;
      for (int lv = 0; lv < grid_Lmax; ++lv) {
        if (a_angular_grid->getL(k) >=
                abs(a_J_orbital->getL(i) - a_J_orbital->getL(j)) &&
            a_angular_grid->getL(k) <=
                a_J_orbital->getL(i) + a_J_orbital->getL(j)) {
          for (int mv = 0; mv <= lv; ++mv) {
            std::complex<double> zeta(0, 0), beta(0, 0);
            if (mv == 0) {
              C3jBlm(a_angular_grid, a_angular_grid->getL(k),
                     a_angular_grid->getM(k), a_J_orbital->getType(),
                     a_J_orbital->getL(i), a_J_orbital->getM(i),
                     a_J_orbital->getType(), a_J_orbital->getL(j),
                     a_J_orbital->getM(j), zeta, beta);
              CJ_2[i * J_orbital_num_channels * num_spherical_harmonics +
                   j * num_spherical_harmonics + k] = zeta;
              CJ_2[j * J_orbital_num_channels * num_spherical_harmonics +
                   i * num_spherical_harmonics + k] = zeta;
              k++;
            } else {
              C3jBlm(a_angular_grid, a_angular_grid->getL(k),
                     a_angular_grid->getM(k), a_J_orbital->getType(),
                     a_J_orbital->getL(i), a_J_orbital->getM(i),
                     a_J_orbital->getType(), a_J_orbital->getL(j),
                     a_J_orbital->getM(j), zeta, beta);
              CJ_2[i * J_orbital_num_channels * num_spherical_harmonics +
                   j * num_spherical_harmonics + k] = zeta;
              CJ_2[j * J_orbital_num_channels * num_spherical_harmonics +
                   i * num_spherical_harmonics + k] = zeta;
              k++;
              CJ_2[i * J_orbital_num_channels * num_spherical_harmonics +
                   j * num_spherical_harmonics + k] = beta;
              CJ_2[j * J_orbital_num_channels * num_spherical_harmonics +
                   i * num_spherical_harmonics + k] = beta;
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

  int local_n =
      (int)floor(a_Nbas * a_ket_bra_orbital->getNumChannels() / a_numprocs);

  m_dvr_rep = std::make_unique<std::complex<double>[]>(
      a_Nbas * a_Nbas * ket_bra_num_channels * ket_bra_num_channels);
  for (int i = 0; i < local_n; ++i) {
    int ket_l = floor(i + id * local_n / a_Nbas);
    int ket_r = i + id * local_n - a_Nbas * ket_l;
    for (int j = 0; j < a_Nbas * ket_bra_num_channels; ++j) {
      int bra_l = floor(j / a_Nbas);
      int bra_r = j - a_Nbas * bra_l;
      std::complex<double> tmp_J(1, 0);
      if (ket_r == bra_r) {
        printf("id %d ket_l %d ket_r %d bra_l %d bra_r %d \n", id, ket_l, ket_r,
               bra_l, bra_r);
        /* for (int lk = 0; lk < a_Nbas * num_spherical_harmonics; ++lk) { */
        /* int k = lk / num_spherical_harmonics; */
        /* int L = lk % num_spherical_harmonics; */
        for (int k = 0; k < a_Nbas; ++k) {
          for (int L = 0; L < num_spherical_harmonics; ++L) {
            /* printf("id %d L %d k %d ket_r %d \n", id, L, k, ket_r); */
            if (a_angular_grid->getL(L) >=
                    abs(a_ket_bra_orbital->getL(ket_l) -
                        a_ket_bra_orbital->getL(bra_l)) &&
                a_angular_grid->getL(L) <= a_ket_bra_orbital->getL(ket_l) +
                                               a_ket_bra_orbital->getL(bra_l)) {
              std::complex<double> tmp_1 =
                  4.0 * M_PI * pow(-1, a_angular_grid->getM(L)) *
                  a_T->getTIXX(a_angular_grid->getL(L) * a_Nbas * a_Nbas +
                               k * a_Nbas + ket_r) /
                  (2.0 * a_angular_grid->getL(L) + 1.0);
              std::complex<double> tmp_2 =
                  CJ_1[ket_l * ket_bra_num_channels * num_spherical_harmonics +
                       bra_l * num_spherical_harmonics + L];
              /* printf("z3 %f \n", tmp_2); */
              for (int l1 = 0; l1 < J_orbital_num_channels; ++l1) {
                for (int l2 = 0; l2 < J_orbital_num_channels; ++l2) {
                  if (a_angular_grid->getL(L) >=
                          abs(a_J_orbital->getL(l1) - a_J_orbital->getL(l2)) &&
                      a_angular_grid->getL(L) <=
                          a_J_orbital->getL(l1) + a_J_orbital->getL(l2)) {
                    tmp_J = tmp_1 * tmp_2 *
                            CJ_2[l1 * J_orbital_num_channels *
                                     num_spherical_harmonics +
                                 l2 * num_spherical_harmonics + L] *
                            a_J_orbital->getPartialWaveRep(l1 * a_Nbas + k) *
                            a_J_orbital->getPartialWaveRep(l2 * a_Nbas + k);
                  }
                }
              }
            }
          }
        }
      }
      m_dvr_rep[id * local_n * a_Nbas * ket_bra_num_channels +
                i * a_Nbas * ket_bra_num_channels + j] = tmp_J;
    }
  }
  /* MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &m_dvr_rep, */
  /*               local_n * a_Nbas * ket_bra_num_channels, MPI_DOUBLE_COMPLEX,
   */
  /*               MPI_COMM_WORLD); */
}

Nuclearoperator::Nuclearoperator(
    const int &a_numprocs, const std::array<int, 2> &a_proccessor_xy,
    const size_t &a_Nbas, std::shared_ptr<AngularGrid> a_angular_grid,
    std::shared_ptr<MOPartialWaveRepresentation> a_ket_bra_orbital,
    std::shared_ptr<MOPartialWaveRepresentation> a_J_orbital,
    std::shared_ptr<Toperator> a_T) {

  size_t grid_Lmax = a_angular_grid->getLmax();
  size_t num_spherical_harmonics = (grid_Lmax - 1) * (grid_Lmax + 1) + 1;
  int test_smallerLmax = floor(grid_Lmax / 2);

  /* Building the angular coefficients */
  size_t ket_bra_num_channels = a_ket_bra_orbital->getNumChannels();
  auto CJ_1 = std::make_unique<std::complex<double>[]>(
      ket_bra_num_channels * ket_bra_num_channels * num_spherical_harmonics);

  for (int i = 0; i < ket_bra_num_channels; ++i) {
    for (int j = 0; j < ket_bra_num_channels; ++j) {
      int k = 0;
      for (int lv = 0; lv < grid_Lmax; ++lv) {
        if (a_angular_grid->getL(k) >=
                abs(a_ket_bra_orbital->getL(i) - a_ket_bra_orbital->getL(j)) &&
            a_angular_grid->getL(k) <=
                a_ket_bra_orbital->getL(i) + a_ket_bra_orbital->getL(j)) {
          for (int mv = 0; mv <= lv; ++mv) {
            std::complex<double> zeta(0, 0), beta(0, 0);
            if (mv == 0) {
              C3jBlm(a_angular_grid, a_angular_grid->getL(k),
                     a_angular_grid->getM(k), a_ket_bra_orbital->getType(),
                     a_ket_bra_orbital->getL(i), a_ket_bra_orbital->getM(i),
                     a_ket_bra_orbital->getType(), a_ket_bra_orbital->getL(j),
                     a_ket_bra_orbital->getM(j), zeta, beta);
              CJ_1[i * ket_bra_num_channels * num_spherical_harmonics +
                   j * num_spherical_harmonics + k] = zeta;
              CJ_1[j * ket_bra_num_channels * num_spherical_harmonics +
                   i * num_spherical_harmonics + k] = zeta;
              k++;
            } else {
              C3jBlm(a_angular_grid, a_angular_grid->getL(k),
                     a_angular_grid->getM(k), a_ket_bra_orbital->getType(),
                     a_ket_bra_orbital->getL(i), a_ket_bra_orbital->getM(i),
                     a_ket_bra_orbital->getType(), a_ket_bra_orbital->getL(j),
                     a_ket_bra_orbital->getM(j), zeta, beta);
              CJ_1[i * ket_bra_num_channels * num_spherical_harmonics +
                   j * num_spherical_harmonics + k] = zeta;
              CJ_1[j * ket_bra_num_channels * num_spherical_harmonics +
                   i * num_spherical_harmonics + k] = zeta;
              k++;
              CJ_1[i * ket_bra_num_channels * num_spherical_harmonics +
                   j * num_spherical_harmonics + k] = beta;
              CJ_1[j * ket_bra_num_channels * num_spherical_harmonics +
                   i * num_spherical_harmonics + k] = beta;
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

  size_t J_orbital_num_channels = a_J_orbital->getNumChannels();
  auto CJ_2 = std::make_unique<std::complex<double>[]>(J_orbital_num_channels *
                                                       J_orbital_num_channels *
                                                       num_spherical_harmonics);

  for (int i = 0; i < J_orbital_num_channels; ++i) {
    for (int j = 0; j < J_orbital_num_channels; ++j) {
      int k = 0;
      for (int lv = 0; lv < grid_Lmax; ++lv) {
        if (a_angular_grid->getL(k) >=
                abs(a_J_orbital->getL(i) - a_J_orbital->getL(j)) &&
            a_angular_grid->getL(k) <=
                a_J_orbital->getL(i) + a_J_orbital->getL(j)) {
          for (int mv = 0; mv <= lv; ++mv) {
            std::complex<double> zeta(0, 0), beta(0, 0);
            if (mv == 0) {
              C3jBlm(a_angular_grid, a_angular_grid->getL(k),
                     a_angular_grid->getM(k), a_J_orbital->getType(),
                     a_J_orbital->getL(i), a_J_orbital->getM(i),
                     a_J_orbital->getType(), a_J_orbital->getL(j),
                     a_J_orbital->getM(j), zeta, beta);
              CJ_2[i * J_orbital_num_channels * num_spherical_harmonics +
                   j * num_spherical_harmonics + k] = zeta;
              CJ_2[j * J_orbital_num_channels * num_spherical_harmonics +
                   i * num_spherical_harmonics + k] = zeta;
              k++;
            } else {
              C3jBlm(a_angular_grid, a_angular_grid->getL(k),
                     a_angular_grid->getM(k), a_J_orbital->getType(),
                     a_J_orbital->getL(i), a_J_orbital->getM(i),
                     a_J_orbital->getType(), a_J_orbital->getL(j),
                     a_J_orbital->getM(j), zeta, beta);
              CJ_2[i * J_orbital_num_channels * num_spherical_harmonics +
                   j * num_spherical_harmonics + k] = zeta;
              CJ_2[j * J_orbital_num_channels * num_spherical_harmonics +
                   i * num_spherical_harmonics + k] = zeta;
              k++;
              CJ_2[i * J_orbital_num_channels * num_spherical_harmonics +
                   j * num_spherical_harmonics + k] = beta;
              CJ_2[j * J_orbital_num_channels * num_spherical_harmonics +
                   i * num_spherical_harmonics + k] = beta;
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
      a_Nbas * a_Nbas * ket_bra_num_channels * ket_bra_num_channels);
  for (int i = 0; i < a_Nbas * ket_bra_num_channels; ++i) {
    int ket_r = i / a_Nbas;
    int ket_l = i % a_Nbas;
    for (int j = 0; j < a_Nbas * ket_bra_num_channels; ++j) {
      int bra_r = j / a_Nbas;
      int bra_l = j % a_Nbas;
      std::complex<double> tmp_J(1, 0);
      if (ket_r == bra_r) {
        for (int lk = 0; lk < a_Nbas * num_spherical_harmonics; ++lk) {
          int k = lk / num_spherical_harmonics;
          int L = lk % num_spherical_harmonics;
          if (a_angular_grid->getL(L) >= abs(a_ket_bra_orbital->getL(ket_l) -
                                             a_ket_bra_orbital->getL(bra_l)) &&
              a_angular_grid->getL(L) <= a_ket_bra_orbital->getL(ket_l) +
                                             a_ket_bra_orbital->getL(bra_l)) {
            std::complex<double> tmp_1 =
                4.0 * M_PI * pow(-1, a_angular_grid->getM(L)) *
                a_T->getTIXX(a_angular_grid->getL(L) * a_Nbas * a_Nbas +
                             k * a_Nbas + ket_r) /
                (2.0 * a_angular_grid->getL(L) + 1.0);
            std::complex<double> tmp_2 =
                CJ_1[ket_l * ket_bra_num_channels * num_spherical_harmonics +
                     bra_l * num_spherical_harmonics + L];
            /* printf("z3 %f", tmp_2); */
            for (int l1 = 0; l1 < J_orbital_num_channels; ++l1) {
              for (int l2 = 0; l2 < J_orbital_num_channels; ++l2) {
                if (a_angular_grid->getL(L) >=
                        abs(a_J_orbital->getL(l1) - a_J_orbital->getL(l2)) &&
                    a_angular_grid->getL(L) <=
                        a_J_orbital->getL(l1) + a_J_orbital->getL(l2)) {
                  tmp_J = tmp_1 * tmp_2 *
                          CJ_2[l1 * J_orbital_num_channels *
                                   num_spherical_harmonics +
                               l2 * num_spherical_harmonics + L];
                  a_J_orbital->getPartialWaveRep(l1 * a_Nbas + k) *
                      a_J_orbital->getPartialWaveRep(l2 * a_Nbas + k);
                }
              }
            }
          }
        }
      }
      m_dvr_rep[i * a_Nbas * ket_bra_num_channels + j] = tmp_J;
    }
  }
}

Nuclearoperator::~Nuclearoperator() {}

std::complex<double> Nuclearoperator::getJ(int index) const {
  return m_dvr_rep[index];
}

void Nuclearoperator::C3jBlm(std::shared_ptr<AngularGrid> a_angular_grid,
                             const int &L, const int &M,
                             const std::string &type1, const int &l1,
                             const int &m1, const std::string &type2,
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
  a_zeta = 4.0 * M_PI * a_zeta;
  a_beta = 4.0 * M_PI * a_beta;
}

std::complex<double> Nuclearoperator::Ylm(const int l, const int m,
                                          const double theta,
                                          const double phi) {
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

double Nuclearoperator::Blm(const std::string type, const int l, const int m,
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
