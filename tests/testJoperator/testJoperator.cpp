/* #include "../../extern/mxx/include/mxx/bcast.hpp" */
#include "AngularGrid.hpp"
#include "FEMDVR.hpp"
#include "Global_DataTypes.hpp"
#include "Toperator.hpp"
/* #include "mxx/comm.hpp" */
/* #include "mxx/distribution.hpp" */

#include <iostream>
#include <mpi.h>
#include <stdio.h>

int main(int argc, char **argv) {

  int ndvr = 5;
  int nel = 3, lmax = 0, lmax_times_2 = 3;
  int numprocs, id;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &id);
  std::cout << "numper of processes " << numprocs << "  id " << id << "\n";

  /* auto LobattoQuad = std::make_unique<Quadrature_Lobatto>(ndvr); */
  /* auto RadauQuad = std::make_unique<Quadrature_Radau>(ndvr); */
  /* auto scatterGrid = std::make_shared<FEMDVR>(std::move(LobattoQuad), */
  /*                                             std::move(RadauQuad), nel); */

  /* Toperator T(scatterGrid, lmax_times_2); */

  return 0;
}
