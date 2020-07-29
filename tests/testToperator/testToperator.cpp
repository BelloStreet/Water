#include "AngularGrid.hpp"
#include "FEMDVR.hpp"
#include "Global_DataTypes.hpp"
#include "Toperator.hpp"

#include <iostream>
#include <stdio.h>

int main(int argc, char **argv) {

  int ndvr = 5;
  int nel = 3, lmax = 0, lmax_times_2 = 3;

  auto LobattoQuad = std::make_unique<Quadrature_Lobatto>(ndvr);
  auto RadauQuad = std::make_unique<Quadrature_Radau>(ndvr);
  auto scatterGrid = std::make_shared<FEMDVR>(std::move(LobattoQuad),
                                              std::move(RadauQuad), nel);

  Toperator T(scatterGrid, lmax_times_2);

  for (int l = 0; l < lmax_times_2; ++l) {
    for (int i = 0; i < scatterGrid->getNbas(); ++i) {
      for (int j = 0; j < scatterGrid->getNbas(); ++j) {
        std::cout << T.getTIXX(l * scatterGrid->getNbas() *
                                   scatterGrid->getNbas() +
                               i * scatterGrid->getNbas() + j)
                  << "\n";
      }
    }
  }
  /* for (int l = 0; l < lmax_times_2; ++l) { */
  /*   for (int i = 0; i < scatterGrid->getNbas(); ++i) { */
  /*     std::cout << T.getTXX(l * realGrid->getNbas() + */
  /*     i * realGrid->getNbas() + i) */
  /*     << "\n"; */
  /*   } */
  /* } */

  return 0;
}
