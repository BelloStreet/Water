#ifndef _MOPARTIALWAVEREPRESENTATION_HPP_
#define _MOPARTIALWAVEREPRESENTATION_HPP_

#include "FEMDVR.hpp"
#include <complex>
#include <vector>

#include "FEMDVR.hpp"
#include "Global_DataTypes.hpp"

class MOPartialWaveRepresentation {
public:
  /// Constructor with default lmax=6 and lebedev quadrature = 2354, pass in
  /// different values to change. Have to change lebedev function to match new
  /// order

  MOPartialWaveRepresentation(std::shared_ptr<FEMDVR> a_grid,
                              const molecule_t &a_molecule,
                              const std::string &a_irrep,
                              const unsigned int &a_lmax = 19,
                              const int &a_angular_order = 2354);

  /// Destructor
  ~MOPartialWaveRepresentation();

  /// Getter for the Partial Wave Representation for given molecular orbital.
  std::complex<double> *getPartialWaveRep() const;

  /// Getter for the number of channels generated from lmax
  unsigned int getNumChannels() const;

  // Good practice is to make helper functions private
private:
  unsigned int m_lmax, m_num_channels;
  int m_angular_order;
  std::unique_ptr<std::complex<double>[]> m_partial_wave_rep_orbital;
  double Blm(const std::string, const int, const int, const double,
             const double);
};
#endif
