#include "AngularGrid.hpp"
#include "FEMDVR.hpp"
#include "Global_DataTypes.hpp"
#include "MOPartialWaveRepresentation.hpp"
#include "MoldenParser.hpp"
#include "Quadrature_Lobatto.hpp"
#include "Quadrature_Radau.hpp"

#include <iostream>
#include <stdio.h>

int main(int argc, char **argv) {

  molecule_t molecule;
  moldenParser(argv[1], molecule);

  int ndvr = 10;
  int nel = 13;

  auto LobattoQuad = std::make_unique<Quadrature_Lobatto>(ndvr);
  auto realGrid = std::make_shared<FEMDVR>(std::move(LobattoQuad), nel);

  auto angularGrid = std::make_shared<AngularGrid>();
  std::cout << angularGrid->getNumChannels() << "\n";
  /* std::cout << angularGrid->getLmax() << "\n"; */

  int which_MO = 0;
  MOPartialWaveRepresentation MOPWRep(realGrid, angularGrid, molecule,
                                      which_MO);

  std::complex<double> check_norm = 0.0;
  for (int i = 0; i < realGrid->getNbas() * MOPWRep.getNumChannels(); i++) {
    check_norm +=
        MOPWRep.getPartialWaveRep(i) * conj(MOPWRep.getPartialWaveRep(i));
  }
  std::cout << check_norm << "\n";

  return 0;
}
