// cSpell:disable
#ifndef _TOPERATOR_HPP_
#define _TOPERATOR_HPP_

#include <complex>

#include "FEMDVR.hpp"

class Toperator {
public:
  /// Toperator: $T = \Del^2 + \frac{l(l+1)}{2x^2}$
  Toperator(std::shared_ptr<FEMDVR> a_femdvr_grid, const int &a_lmax_times_2);

  /// Destructor
  ~Toperator();

  /// Getter for the T operator in the DVR representation.
  std::complex<double> getTXX(int index) const;

  /// Getter for the inverse of the T operator in the DVR representation.
  /// Used in the poison solution for $\frac{1}{|r_1-r_2|}.
  std::complex<double> getTIXX(int index) const;

private:
  std::unique_ptr<std::complex<double>[]> m_dvr_rep, m_inverse_dvr_rep;
};
#endif
