#ifndef _MOPARTIALWAVEREPRESENTATION_HPP_
#define _MOPARTIALWAVEREPRESENTATION_HPP_

#include <complex>
#include <vector>

#include "AngularGrid.hpp"
#include "FEMDVR.hpp"
#include "Global_DataTypes.hpp"

class MOPartialWaveRepresentation
{
public:
  /// Constructor with default lmax=6 and lebedev quadrature = 2354, pass in
  /// different values to change. Have to change lebedev function to match new
  /// lebedev quadrature order

  MOPartialWaveRepresentation(std::shared_ptr<FEMDVR>      a_femdvr_grid,
                              std::shared_ptr<AngularGrid> a_angular_grid,
                              const molecule_t &           a_molecule,
                              const unsigned int &         a_which_MO);

  /// Destructor
  ~MOPartialWaveRepresentation();

  /// Getter for the Partial Wave Representation for given molecular orbital.
  std::complex<double>
  getPartialWaveRep(int index) const;

  /// Getter for the number of channels generated from lmax.
  size_t
  getNumChannels() const;

  /// Getter for the type of Blm, either sine or cosine like.
  std::string
  getType() const;

  /// Getter for the L quantum number generated for this irrep.
  int
  getL(int index) const;

  /// Getter for the M quantum number generated for this irrep.
  int
  getM(int index) const;

  // Good practice is to make helper functions private
private:
  int                    m_lmax;
  size_t                 m_num_channels;
  std::string            m_type;
  std::unique_ptr<int[]> m_quantum_number_l, m_quantum_number_m;
  std::unique_ptr<std::complex<double>[]> m_partial_wave_rep_orbital;
  double
  Blm(const std::string, const int, const int, const double, const double);
};
#endif
