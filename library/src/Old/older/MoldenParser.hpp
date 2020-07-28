#ifndef _MOLDENPARSER_HPP_
#define _MOLDENPARSER_HPP_
#include "Global_DataTypes.hpp"
#include <complex>
#include <memory>
#include <vector>

class MoldenParser {
public:
  /// Class to make parsing Molden files more transparent.

  MoldenParser(const std::string &filename);

  /// Destructor.
  ~MoldenParser();

  /// Getter for atoms' information.

  molecule_t getMoleculeInfo() const;

  /// Getter for Gaussian orbital information.

  orbital_t getGaussiansInfo() const;

  /// Getter for molecular orbital information.

  std::vector<MO_t> getMOsInfo() const;

private:
  molecule_t m_molecule;
  orbital_t m_gaussians;
  std::vector<MO_t> m_MOs;
};
#endif
