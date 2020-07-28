#include "Global_DataTypes.hpp"
#include "MoldenParser.hpp"

#include <iostream>
#include <stdio.h>

int main(int argc, char **argv) {

  /* MoldenParser molden_file(argv[1]); */
  molecule_t molecule;
  moldenParser(argv[1], molecule);

  /* molecule_t molecule = molden_file.getMoleculeInfo(); */
  std::cout << "Point Group " << molecule.pointgroup << "\n";
  for (auto i : molecule.atoms) {
    std::cout << "Symbol " << i.symbol << "\n";
    std::cout << "Charge " << i.charge << "\n";
    std::cout << "Number of S orbitals "
              << i.AOorbital.s_orbitals[0].coeff_s.size() << "\n";
    for (auto j : i.AOorbital.s_orbitals[0].coeff_s) {
      std::cout << j << "\n";
    }
  }
  for (auto i : molecule.molecular_orbitals) {
    /* std::cout << i.symmetry << "\n"; */
  }

  return 0;
}
