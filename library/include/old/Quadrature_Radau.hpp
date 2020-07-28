// cSpell:disable
#ifndef _QUADRATURE_RADAU_HPP_
#define _QUADRATURE_RADAU_HPP_

#include <memory>

#include "Quadrature.hpp"

class Quadrature_Radau : public Quadrature {
public:
  /// Radau for complex tail.
  Quadrature_Radau(int &a_dvr);

  virtual ~Quadrature_Radau();

  virtual std::complex<double> *getPoints() const;

  virtual std::complex<double> *getWeights() const;

  virtual void print() const;

private:
  int m_kind, m_kpts;
  double m_alpha, m_beta;
  std::unique_ptr<std::complex<double>[]> m_points, m_weights;
};
#endif
