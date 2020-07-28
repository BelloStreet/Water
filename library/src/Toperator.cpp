#include "../include/Toperator.hpp"
#include <iostream>
#include <mkl_lapacke.h>

Toperator::Toperator(std::shared_ptr<FEMDVR> a_femdvr_grid,
                     std::shared_ptr<AngularGrid> a_angular_grid) {

  m_dvr_rep = std::make_unique<std::complex<double>[]>(
      a_femdvr_grid->getNbas() * a_femdvr_grid->getNbas() *
      a_angular_grid->getLmax());
  for (int l = 0; l < a_angular_grid->getLmax(); ++l) {
    for (int i = 0; i < a_femdvr_grid->getNbas(); ++i) {
      for (int j = 0; j < a_femdvr_grid->getNbas(); ++j) {
        if (i == j) {
          m_dvr_rep[l * a_femdvr_grid->getNbas() * a_femdvr_grid->getNbas() +
                    i * a_femdvr_grid->getNbas() + j] =
              a_femdvr_grid->getLaplacian(i * a_femdvr_grid->getNbas() + j) +
              (std::complex<double>)(l * (l + 1)) /
                  pow(a_femdvr_grid->getPoint(i), 2);
        } else {
          m_dvr_rep[l * a_femdvr_grid->getNbas() * a_femdvr_grid->getNbas() +
                    i * a_femdvr_grid->getNbas() + j] =
              a_femdvr_grid->getLaplacian(i * a_femdvr_grid->getNbas() + j);
        }
      }
    }
  }
  m_inverse_dvr_rep = std::make_unique<std::complex<double>[]>(
      a_femdvr_grid->getNbas() * a_femdvr_grid->getNbas() *
      a_angular_grid->getLmax());
  for (int l = 0; l < a_angular_grid->getLmax(); ++l) {
    auto tmp_laplacian = std::make_unique<std::complex<double>[]>(
        a_femdvr_grid->getNbas() * a_femdvr_grid->getNbas());
    auto tmp_swap_value = std::make_unique<std::complex<double>[]>(
        a_femdvr_grid->getNbas() * a_femdvr_grid->getNbas());
    for (int i = 0; i < a_femdvr_grid->getNbas(); ++i) {
      for (int j = 0; j < a_femdvr_grid->getNbas(); ++j) {
        tmp_laplacian[i * a_femdvr_grid->getNbas() + j] =
            m_dvr_rep[l * a_femdvr_grid->getNbas() * a_femdvr_grid->getNbas() +
                      i * a_femdvr_grid->getNbas() + j];
        if (i == j) {
          tmp_swap_value[i * a_femdvr_grid->getNbas() + j] = 1.0;
        } else {
          tmp_swap_value[i * a_femdvr_grid->getNbas() + j] = 0.0;
        }
      }
    }
    auto tmp_ipiv = std::make_unique<int[]>(a_femdvr_grid->getNbas());

    auto info = LAPACKE_zgesv(
        LAPACK_COL_MAJOR, a_femdvr_grid->getNbas(), a_femdvr_grid->getNbas(),
        reinterpret_cast<MKL_Complex16 *>(tmp_laplacian.get()),
        a_femdvr_grid->getNbas(), tmp_ipiv.get(),
        reinterpret_cast<MKL_Complex16 *>(tmp_swap_value.get()),
        a_femdvr_grid->getNbas());
    for (int i = 0; i < a_femdvr_grid->getNbas(); ++i) {
      for (int j = 0; j < a_femdvr_grid->getNbas(); ++j) {
        m_inverse_dvr_rep[l * a_femdvr_grid->getNbas() *
                              a_femdvr_grid->getNbas() +
                          i * a_femdvr_grid->getNbas() + j] =
            tmp_swap_value[i * a_femdvr_grid->getNbas() + j];
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
