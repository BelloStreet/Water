#ifndef _ANGULAREGRID_HPP_
#define _ANGULAREGRID_HPP_

#include <memory>

class AngularGrid {
public:
  /// Wrapper for the lebedev quadrature that creates the angular grid.
  /// The order pased in must match the lebedev routined called inside!
  AngularGrid(const unsigned int &a_lmax = 19,
              const unsigned int &a_angular_order = 2354);

  /// Destructor
  ~AngularGrid();

  /// Getter for angular weight used for later spherical integraton.
  double getAngularWeight(int index) const;

  /// Getter for the X points in the angular integration.
  double getX(int index) const;

  /// Getter for the Y points in the angular integration.
  double getY(int index) const;

  /// Getter for the Z points in the angular integration.
  double getZ(int index) const;

  /// Getter for lmax
  unsigned int getLmax() const;

  /// Getter for lmax
  unsigned int getNumChannels() const;

  /// Getter for the L quantum number generated for the angular grid's max
  /// angular momentum
  int getL(int index) const;

  /// Getter for the M quantum number generated for the angular grid's max
  /// angular momentum
  int getM(int index) const;

  unsigned int getAngularOrder() const;

  // Good practice is to make helper functions private
private:
  unsigned int m_lmax, m_angular_order, m_num_channels;
  std::unique_ptr<int[]> m_quantum_number_l, m_quantum_number_m;
  std::unique_ptr<double[]> m_angular_weight, m_angular_x, m_angular_y,
      m_angular_z;
};
#endif
