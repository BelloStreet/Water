#ifndef _JOPERATOR_HPP_
#define _JOPERATOR_HPP_

#include <complex>

#include "MOPartialWaveRepresentation.hpp"
#include "Toperator.hpp"

class Joperator
{
public:
  /// The coulomb operator takes in the angular grid, a ket and bra orbitals
  /// with the T operator that provides $T^{-1}$ for the $\frac{1}{|r_1-r_2|}$
  /// term.
  Joperator(const size_t &                               a_Nbas,
            std::shared_ptr<AngularGrid>                 a_angular_grid,
            std::shared_ptr<MOPartialWaveRepresentation> a_ket_orbital,
            std::shared_ptr<MOPartialWaveRepresentation> a_bra_orbital,
            std::shared_ptr<Toperator>                   a_T);

  Joperator(const int &                                  a_numprocs,
            const int &                                  id,
            const size_t &                               a_Nbas,
            std::shared_ptr<AngularGrid>                 a_angular_grid,
            std::shared_ptr<MOPartialWaveRepresentation> a_ket_orbital,
            std::shared_ptr<MOPartialWaveRepresentation> a_bra_orbital,
            std::shared_ptr<Toperator>                   a_T);

  Joperator(const int &                                  a_numprocs,
            const std::array<int, 2> &                   a_proccessor_xy,
            const size_t &                               a_Nbas,
            std::shared_ptr<AngularGrid>                 a_angular_grid,
            std::shared_ptr<MOPartialWaveRepresentation> a_ket_orbital,
            std::shared_ptr<MOPartialWaveRepresentation> a_bra_orbital,
            std::shared_ptr<Toperator>                   a_T);

  /// Destructor
  ~Joperator();

  /// Getter for the T operator in the DVR representation.
  std::complex<double>
  getJ(int index) const;

private:
  std::unique_ptr<std::complex<double>[]> m_dvr_rep, m_CJ1, m_CJ2;
  void
  C3jBlm(std::shared_ptr<AngularGrid>,
         const int &,
         const int &,
         const std::string &,
         const int &,
         const int &,
         const std::string &,
         const int &,
         const int &,
         std::complex<double> &,
         std::complex<double> &);
  std::complex<double>
  Ylm(const int, const int, const double, const double);
  double
  Blm(const std::string, const int, const int, const double, const double);
};
#endif
