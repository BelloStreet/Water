// cSpell:disable
#ifndef _QUADRATURE_LOBATTO_HPP_
#define _QUADRATURE_LOBATTO_HPP_

#include <memory>

#include "Quadrature.hpp"

class Quadrature_Lobatto : public Quadrature {
public:
  /// Lobatto for real part of the grid.
  Quadrature_Lobatto(int &a_dvr);

  virtual ~Quadrature_Lobatto();

  virtual std::complex<double> getPoint(int index) const;

  virtual std::complex<double> getWeight(int index) const;

  virtual void print() const;

private:
  int m_kind, m_kpts, m_dvr;
  double m_alpha, m_beta;
  std::unique_ptr<std::complex<double>[]> m_points, m_weights;
};
#endif
