// cSpell:disable
#ifndef _HAMILTONIAN_HPP_
#define _HAMILTONIAN_HPP_

#include <complex>

#include "FEMDVR.hpp"
#include "Global_DataTypes.hpp"
#include "MOPartialWaveRepresentation.hpp"

class Hamiltonian {
public:
  /// General Hamiltonian
  /// There should be a constructor for a bound state Hamiltonian.
  Hamiltonian(std::shared_ptr<FEMDVR> a_grid,
              std::shared_ptr<MOPartialWaveRepresentation> a_MO_orbitals,
              molecule_t a_molecule);

  /// Should be a scatter Hamiltonian constructor.
  /* Hamiltonian(std::shared_ptr<FEMDVR> a_grid, molecule_t a_molecule, */
  /* ); */

  /// Destructor
  ~Hamiltonian();

  /// Getter for the Hamiltonian in the DVR representation.
  const std::complex<double> getHamiltionianDVRRep(int index) const;

  /// Getter for the J operator DVR representation.
  const std::complex<double> getJ(int index) const;

  /// Getter for the K operator DVR representation.
  const std::complex<double> getK(int index) const;

private:
  std::unique_ptr<std::complex<double>[]> m_Hamiltonian, m_J_operator,
      m_K_operator;
  std::shared_ptr<FEMDVR> m_grid;
};
#endif
