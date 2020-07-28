#include "AngularGrid.hpp"
#include "FEMDVR.hpp"
#include "Global_DataTypes.hpp"
#include "Toperator.hpp"

#include <iostream>
#include <stdio.h>

int main(int argc, char **argv) {

  int ndvr = 5;
  int nel = 3, lmax = 3;

  auto LobattoQuad = std::make_unique<Quadrature_Lobatto>(ndvr);
  auto realGrid = std::make_shared<FEMDVR>(std::move(LobattoQuad), nel);
  auto angularGrid = std::make_shared<AngularGrid>(lmax);

  Toperator laplacian(realGrid, angularGrid);

  for (int l = 0; l < angularGrid->getLmax(); ++l) {
    for (int i = 0; i < realGrid->getNbas(); ++i) {
      for (int j = 0; j < realGrid->getNbas(); ++j) {
        std::cout << laplacian.getTIXX(l * realGrid->getNbas() *
                                           realGrid->getNbas() +
                                       i * realGrid->getNbas() + j)
                  << "\n";
      }
    }
  }
  /* for (int l = 0; l < angularGrid->getLmax(); ++l) { */
  /* for (int i = 0; i < realGrid->getNbas(); ++i) { */
  /* std::cout << laplacian.getTXX(l * realGrid->getNbas() + */
  /*                                  i * realGrid->getNbas() + i) */
  /*           << "\n"; */
  /* } */
  /* } */

  return 0;
}
