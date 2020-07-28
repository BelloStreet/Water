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

  /// Get number of basis functions
  size_t &getNbas();
  /* const unsigned int &getNbas(); */

  void print() const;

private:
  unsigned int m_Nelem;
  size_t m_Nbas;
  /* TODO:Either get rid of bounds variables or make them pointers */
  std::vector<std::complex<double>> m_realbounds, m_complexbounds;
  std::unique_ptr<std::complex<double>[]> m_points, m_weights, m_laplacian;
};
#endif
