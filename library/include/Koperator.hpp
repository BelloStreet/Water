#ifndef _KOPERATOR_HPP_
#define _KOPERATOR_HPP_

#include <complex>

#include "MOPartialWaveRepresentation.hpp"
#include "Toperator.hpp"

class Koperator {
public:
  /// The exchange operator takes in the angular grid, a ket and bra orbitals
  /// with the T operator that provides $T^{-1}$ for the $\frac{1}{|r_1-r_2|}$
  /// term.
  Koperator(std::shared_ptr<AngularGrid> a_angular_grid,
            std::shared_ptr<MOPartialWaveRepresentation> a_ket_orbital,
            std::shared_ptr<MOPartialWaveRepresentation> a_bra_orbital,
            std::shared_ptr<Toperator> a_T);

  /// Destructor
  ~Koperator();

  /// Getter for the T operator in the DVR representation.
  std::complex<double> getK(int index) const;

private:
  std::unique_ptr<std::complex<double>[]> m_dvr_rep;
  void C3jBlm(std::shared_ptr<AngularGrid>, const int &, const int &,
              const std::string &, const int &, const int &,
              const std::string &, const int &, const int &,
              std::complex<double> &, std::complex<double> &);
  std::complex<double> Ylm(const int, const int, const double, const double);
  double Blm(const std::string, const int, const int, const double,
             const double);
};
#endif
