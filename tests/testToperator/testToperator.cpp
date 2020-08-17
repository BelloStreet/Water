/* #include "AngularGrid.hpp" */
/* #include "FEMDVR.hpp" */
/* #include "Global_DataTypes.hpp" */
#include "Toperator.hpp"

#include <iostream>
#include <stdio.h>

int main(int argc, char **argv) {

  int ndvr = 5;
  int nel = 3, lmax = 0, lmax_times_2 = 3;

  auto rLobattoQuad = std::make_unique<Quadrature_Lobatto>(ndvr);
  auto sLobattoQuad = std::make_unique<Quadrature_Lobatto>(ndvr);
  auto RadauQuad = std::make_unique<Quadrature_Radau>(ndvr);
  auto realGrid = std::make_shared<FEMDVR>(std::move(rLobattoQuad), nel);
  auto scatterGrid = std::make_shared<FEMDVR>(std::move(sLobattoQuad),
                                              std::move(RadauQuad), nel);

  Toperator T(realGrid, lmax_times_2);
  /* Toperator T(scatterGrid, lmax_times_2); */

  printf(" lmax %d nbas %zu \n", lmax_times_2, realGrid->getNbas());
  for (int l = 0; l < lmax_times_2; ++l) {
    for (int i = 0; i < realGrid->getNbas(); ++i) {
      for (int j = 0; j < realGrid->getNbas(); ++j) {
        std::cout << T.getTIXX(l * realGrid->getNbas() * realGrid->getNbas() +
                               i * realGrid->getNbas() + j)
                  << "\n";
      }
    }
  }
  /* printf("lamax %d nbas %zu \n", lmax_times_2, realGrid->getNbas()); */
  /* for (int l = 0; l < lmax_times_2; ++l) { */
  /*   for (int i = 0; i < realGrid->getNbas(); ++i) { */
  /*     std::cout << T.getTXX(l * realGrid->getNbas() * realGrid->getNbas() +
   */
  /*                           i * realGrid->getNbas() + i) */
  /*               << "\n"; */
  /*   } */
  /* } */

  return 0;
}
