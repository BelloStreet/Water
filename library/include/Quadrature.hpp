// cSpell:disable
#ifndef _QUADRATURE_HPP_
#define _QUADRATURE_HPP_

#include <complex>

extern "C" {
void gaussq_(int *, int *, double *, double *, int *, double *, double *,
             double *, double *);
}

class Quadrature {
public:
  /// Quadrature base class that wraps fortran routine gaussq and defines the
  /// interface.
  Quadrature(int &a_dvr) { m_dvr = a_dvr; }

  virtual ~Quadrature(){};

  virtual std::complex<double> getPoint(int index) const = 0;

  virtual std::complex<double> getWeight(int index) const = 0;

  virtual void print() const = 0;

  // Access function
  const int &getDVROrder() { return m_dvr; }

protected:
  int m_dvr;
};
#endif
