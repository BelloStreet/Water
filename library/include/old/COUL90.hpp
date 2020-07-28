// cSpell:disable
#ifndef _COUL90_HPP_
#define _COUL90_HPP_

#include <complex>
#include <vector>

#include "FEMDVR.hpp"

class COUL90 {
public:
  /// A wrapper for the COUL90 Fortran routines that computes hypergeometric
  /// functions
  COUL90(std::shared_ptr<FEMDVR> a_grid, int &a_charge, const double a_k);

  /// Destructor
  ~COUL90();

  /// Getter for the representation of a coulomb wave on the DVR grid
  std::complex<double> *dvrRep() const;

  /// Getter for the coulomb wave
  std::complex<double> *getPsik() const;

  /// Getter for the coulomb wave derivative
  std::complex<double> *getPsikprime() const;

private:
  double m_k, m_eta;
  std::unique_ptr<std::complex<double>[]> m_Psik, m_Psikprime, m_dvrRep;
  std::shared_ptr<FEMDVR> m_grid;
};
#endif
