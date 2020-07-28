// cSpell:disable
#ifndef _JOPERATOR_HPP_
#define _JOPERATOR_HPP_

#include <complex>

#include "FEMDVR.hpp"
#include "Global_DataTypes.hpp"
#include "MOPartialWaveRepresentation.hpp"

class Joperator {
public:
  /// General Hamiltonian
  /// There should be a constructor for a bound state Hamiltonian.
  Joperator(std::shared_ptr<FEMDVR> a_femdvr_grid,
            std::shared_ptr<AngularGrid> a_angular_grid);

  /// Destructor
  ~Joperator();

  /// Getter for the T operator in the DVR representation.
  std::complex<double> getTXX(int index) const;

  /// Getter for the inverse of the T operator in the DVR representation.
  std::complex<double> getTIXX(int index) const;

private:
  std::unique_ptr<std::complex<double>[]> m_dvr_rep, m_inverse_dvr_rep;
};
#endif
