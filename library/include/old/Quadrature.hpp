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

  virtual std::complex<double> *getPoints() const = 0;

  virtual std::complex<double> *getWeights() const = 0;

  virtual void print() const = 0;

  // Access function
  const unsigned int &getDVROrder() { return m_dvr; }

protected:
  unsigned int m_dvr;
};
#endif
