#ifndef _FEMDVR_HPP_
#define _FEMDVR_HPP_
#include <complex>
#include <memory>
#include <vector>

#include "Quadrature_Lobatto.hpp"
#include "Quadrature_Radau.hpp"

class FEMDVR {
public:
  /// Interface class for FEM-DVR grid

  /** Constructor arguments: lobatto and radau order, number of elements.
  Passing the quadratures in as a unique_ptr to naturally take ownership.
  This requires the caller to use move semantics because you cannot copy a
  unique_ptr
  */

  /// Constructor for real bound grid.
  FEMDVR(std::unique_ptr<Quadrature_Lobatto> a_lobattoQuad,
         const unsigned int &a_Nelem);

  /// Constructor for complex scattering grid with Radau tail.
  FEMDVR(std::unique_ptr<Quadrature_Lobatto> a_lobattoQuad,
         std::unique_ptr<Quadrature_Radau> a_radauQuad,
         const unsigned int &a_Nelem);

  /// Destructor.
  ~FEMDVR();

  /// Getter for points.
  std::complex<double> getPoint(int index) const;

  /// Getter for weights.
  std::complex<double> getWeight(int index) const;

  /// Getter for laplacian.
  std::complex<double> getLaplacian(int index) const;

  /// Getter for the real boundaries.
  /* TODO::Must change the function prototype once correct for boundary logic */
  std::complex<double> getRealBoundary(int index) const;

  /// Get the Lobotto order used
  size_t getLobattoOrder() const;

  /// Get the Radau order used. This is set to zero for real grid.
  size_t getRaduaOrder() const;

  /// Get the number of elements
  size_t getNElements() const;

  /// Get number of total basis functions, including complex contour if present.
  size_t getNbas() const;

  /// Get number of real basis functions
  size_t getNRealbas() const;

  /// Get R0: the point where the grid is complex scaled.
  /// This point is set to zero if the grid is completely real.
  std::complex<double> getR0() const;

  /// Get the scaling factor alpha used in the radau grid mappings.
  double getAlphaRad() const;

  /// Get the angle used to rotate by into the complex plane.
  double getTheta() const;

  /// Get the complex rotated factor.
  std::complex<double> getEit() const;

  void print() const;

private:
  unsigned int m_Lobatto_order, m_Radau_order, m_Nelem, m_R0_index;
  size_t m_NRealbas, m_Nbas;
  double m_alphaRad, m_theta;
  std::complex<double> m_R0, m_eit;
  /* TODO:Either get rid of bounds variables or make them pointers */
  std::vector<std::complex<double>> m_realbounds, m_complexbounds;
  std::unique_ptr<std::complex<double>[]> m_points, m_weights, m_laplacian;
};
#endif
