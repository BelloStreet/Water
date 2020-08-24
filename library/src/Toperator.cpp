#include "../include/Toperator.hpp"
#include <iostream>
#include <mkl_lapacke.h>

Toperator::Toperator(std::shared_ptr<FEMDVR> a_femdvr_grid,
                     const int &a_lmax_times_2) {
  size_t nbas = a_femdvr_grid->getNbas();
  m_dvr_rep =
      std::make_unique<std::complex<double>[]>(nbas * nbas * a_lmax_times_2);
  for (int l = 0; l < a_lmax_times_2; ++l) {
    for (int i = 0; i < nbas; ++i) {
      for (int j = 0; j < nbas; ++j) {
        if (i == j) {
          m_dvr_rep[l * nbas * nbas + i * nbas + j] =
              a_femdvr_grid->getLaplacian(i * nbas + j) +
              (std::complex<double>)(l * (l + 1)) /
                  pow(a_femdvr_grid->getPoint(j), 2);
        } else {
          m_dvr_rep[l * nbas * nbas + i * nbas + j] =
              a_femdvr_grid->getLaplacian(i * nbas + j);
        }
      }
    }
  }
  m_inverse_dvr_rep =
      std::make_unique<std::complex<double>[]>(nbas * nbas * a_lmax_times_2);
  for (int l = 0; l < a_lmax_times_2; ++l) {
    auto tmp_laplacian = std::make_unique<std::complex<double>[]>(nbas * nbas);
    auto tmp_swap_value = std::make_unique<std::complex<double>[]>(nbas * nbas);
    for (int i = 0; i < nbas; ++i) {
      for (int j = 0; j < nbas; ++j) {
        tmp_laplacian[i * nbas + j] = m_dvr_rep[l * nbas * nbas + i * nbas + j];
        if (i == j) {
          tmp_swap_value[i * nbas + j] = 1.0;
        } else {
          tmp_swap_value[i * nbas + j] = 0.0;
        }
      }
    }
    auto tmp_ipiv = std::make_unique<int[]>(nbas);

    auto info = LAPACKE_zgesv(
        LAPACK_COL_MAJOR, nbas, nbas,
        reinterpret_cast<MKL_Complex16 *>(tmp_laplacian.get()), nbas,
        tmp_ipiv.get(), reinterpret_cast<MKL_Complex16 *>(tmp_swap_value.get()),
        nbas);
    for (int i = 0; i < nbas; ++i) {
      for (int j = 0; j < nbas; ++j) {
        m_inverse_dvr_rep[l * nbas * nbas + i * nbas + j] =
            tmp_swap_value[i * nbas + j];
      }
    }
  }

  std::complex<double> surface_term;
  int num_elements = a_femdvr_grid->getNElements();
  int radau_order = a_femdvr_grid->getRaduaOrder();
  int real_Nbas = a_femdvr_grid->getNRealbas();
  double alpha = a_femdvr_grid->getAlphaRad();
  std::complex<double> eit = a_femdvr_grid->getEit();
  std::complex<double> R0 = a_femdvr_grid->getR0();
  if (a_femdvr_grid->getRaduaOrder() != 0) {
    surface_term = 2.0 * static_cast<double>(radau_order) * eit / alpha + R0;
  } else {
    /* Bill says Frank is right! */
    /* surface_term = */
    /*     a_femdvr_grid->getRealBoundary(a_femdvr_grid->getNElements()); */
    /* TODO: Roger surface_term is different from the commented Frank surface
     * term */
    surface_term = a_femdvr_grid->getPoint(real_Nbas - 1);
  }

  for (int l = 0; l < a_lmax_times_2; ++l) {
    for (int i = 0; i < nbas; ++i) {
      std::complex<double> tmp_i_factor =
          1.0 /
          (a_femdvr_grid->getPoint(i) * pow(a_femdvr_grid->getWeight(i), 0.5));
      if ((radau_order != 0) && i > (real_Nbas)) {
        tmp_i_factor *=
            exp(-alpha * conj(eit) * (a_femdvr_grid->getPoint(i) - R0));
      }
      for (int j = 0; j < nbas; ++j) {
        std::complex<double> tmp_j_factor =
            1.0 / (a_femdvr_grid->getPoint(j) *
                   pow(a_femdvr_grid->getWeight(j), 0.5));
        if ((radau_order != 0) && j > (real_Nbas)) {
          tmp_j_factor *=
              exp(-alpha * conj(eit) * (a_femdvr_grid->getPoint(j) - R0));
        }
        std::complex<double> tmp_value =
            (2.0 * l + 1.0) *
            m_inverse_dvr_rep[l * nbas * nbas + i * nbas + j] * tmp_i_factor *
            tmp_j_factor;
        tmp_value =
            tmp_value +
            (pow(a_femdvr_grid->getPoint(i) * a_femdvr_grid->getPoint(j), l)) /
                (pow(surface_term, 2.0 * l + 1.0));
        m_inverse_dvr_rep[l * nbas * nbas + i * nbas + j] = tmp_value;
      }
    }
  }
}

Toperator::~Toperator() {}

std::complex<double> Toperator::getTXX(int index) const {
  return m_dvr_rep[index];
}

std::complex<double> Toperator::getTIXX(int index) const {
  return m_inverse_dvr_rep[index];
}
