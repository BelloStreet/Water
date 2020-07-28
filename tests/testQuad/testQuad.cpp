#include <chrono>
#include <iostream>

#include "Quadrature_Lobatto.hpp"
#include "Quadrature_Radau.hpp"

std::chrono::time_point<std::chrono::steady_clock> start_time;
std::chrono::time_point<std::chrono::steady_clock> end_time;
static std::chrono::duration<double> time0;

int main(int argc, char **argv) {
  int ndvr = 5, mdvr = 5;

  start_time = std::chrono::steady_clock::now();
  Quadrature_Lobatto gridLobatto(ndvr);
  Quadrature_Radau gridRadau(ndvr);
  end_time = std::chrono::steady_clock::now();
  time0 += end_time - start_time;
  double seconds = time0.count();
  std::cout << "Radau finished in microseconds " << seconds * 1000000 << "\n";
  std::cout << "\n\n";
  /* auto x = gridLobatto.getPoints(); */
  /* for (int i = 0; i < ndvr; ++i) { */
  /* std::cout << "Points " << x[i] << "\n"; */
  /* } */
  /* auto Xlob = gridLobatto.getPoints(); */
  /* for (int i = 0; i < ndvr; ++i) std::cout << Xlob[i] << "\n"; */

  gridLobatto.print();
  /* std::cout<<std::endl; */
  /* gridRadau.print(); */

  return 0;
}

