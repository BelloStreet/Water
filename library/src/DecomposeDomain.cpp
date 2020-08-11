#include "../include/DecomposeDomain.hpp"
#include <math.h>

/**
  Function that decomposes the MPI processors' domain into a 2D grid with the
  equations $coords_y = FLOOR(rank/numprocs); coords_x = rank -
  coords_y*numprocs$ and returns the x, y pair.  Local min and max of an array
  can be then found with the equations $nx_local=nx/numprocs; x_cell_min_local =
  nx_local*coords_x; x_cell_max_local = nx_local*(coords_x+1)-1$ (same equations
  to get the y-dimension).
  */
std::array<int, 2> decomposeDomain(const int &a_numprocs, const int &a_rank) {
  int coords_y = floor(a_rank / a_numprocs);
  int coords_x = a_rank - coords_y * a_numprocs;
  std::array<int, 2> processor_coords = {{coords_x, coords_y}};
  return processor_coords;
}
