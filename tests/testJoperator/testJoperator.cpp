#include "Joperator.hpp"
#include "MoldenParser.hpp"
/* #include "mxx/comm.hpp" */
/* #include "mxx/distribution.hpp" */

#include <iostream>
#include <mpi.h>
#include <stdio.h>

int main(int argc, char **argv) {

  int ndvr = 5;
  int nel = 13, lmax = 0, lmax_times_2 = 1;
  int numprocs, id;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &id);

  molecule_t molecule;
  moldenParser(argv[1], molecule);

  auto LobattoQuad = std::make_unique<Quadrature_Lobatto>(ndvr);
  auto realGrid = std::make_shared<FEMDVR>(std::move(LobattoQuad), nel);
  auto angularGrid = std::make_shared<AngularGrid>();

  int first_a1 = 0, b2 = 2;
  auto first_A1 = std::make_shared<MOPartialWaveRepresentation>(
      realGrid, angularGrid, molecule, first_a1);
  auto B2 = std::make_shared<MOPartialWaveRepresentation>(realGrid, angularGrid,
                                                          molecule, b2);
  auto T = std::make_shared<Toperator>(realGrid, 19);

  Joperator J(realGrid->getNbas(), angularGrid, B2, first_A1, T);
  /* Joperator J(numprocs, id, realGrid->getNbas(), angularGrid, B2, first_A1,
   * T); */

  return 0;
}
