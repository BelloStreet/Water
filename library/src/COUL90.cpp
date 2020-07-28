// cSpell:disable

#include "../include/COUL90.hpp"

#include <gsl/gsl_sf_coulomb.h>

#include <complex>
#include <iostream>
#include <vector>

#include "../include/FEMDVR.hpp"

COUL90::COUL90(std::shared_ptr<FEMDVR> a_grid, int &a_charge,
               const double a_k) {
  m_k = a_k;
  m_eta = -a_charge / a_k;
  m_grid = a_grid;
  double lmin = 0.0;
  int lmax = lmin + a_k;
  auto fc = std::make_unique<double[]>(lmax);
  auto fcp = std::make_unique<double[]>(lmax);
  auto gc = std::make_unique<double[]>(lmax);
  auto gcp = std::make_unique<double[]>(lmax);
  auto exp_F = std::make_unique<double[]>(lmax);
  auto exp_G = std::make_unique<double[]>(lmax);
  m_Psik = std::make_unique<std::complex<double>[]>(a_grid->getNbas() * lmax);
  m_Psikprime =
      std::make_unique<std::complex<double>[]>(a_grid->getNbas() * lmax);
  m_dvrRep = std::make_unique<std::complex<double>[]>(a_grid->getNbas() * lmax);
  int GSL_OVERFLOW;
  for (int i = 1; i < a_grid->getNbas(); ++i) {
    double rval = a_k * a_grid->getPoint(i).real();
    GSL_OVERFLOW = gsl_sf_coulomb_wave_FGp_array(
        lmin, a_k, m_eta, rval, fc.get(), fcp.get(), gc.get(), gcp.get(),
        exp_F.get(), exp_G.get());
    for (int lval = 0; lval < lmax; ++lval) {
      m_dvrRep[i * lval + lval] =
          fc[lval] * sqrt(a_grid->getWeight(i).real()) / m_k;
      m_Psik[i * lval + lval] = fc[lval];
      m_Psikprime[i * lval + lval] = fcp[lval] * m_k;
    }
  }
}

COUL90::~COUL90() {}

inline std::complex<double> COUL90::dvrRep(int index) const {
  return m_dvrRep[index];
}

inline std::complex<double> COUL90::getPsik(int index) const {
  return m_Psik[index];
}

inline std::complex<double> COUL90::getPsikprime(int index) const {
  return m_Psikprime[index];
}
