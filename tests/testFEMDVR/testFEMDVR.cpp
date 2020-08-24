#include "FEMDVR.hpp"
#include <iostream>
#include <memory>

int main(int argc, char **argv) {
  int ndvr = 5, mdvr = 6;
  unsigned int nel = 3;

  auto real_lobattoQuad = std::make_unique<Quadrature_Lobatto>(ndvr);
  FEMDVR realGrid(std::move(real_lobattoQuad), nel);
  realGrid.print();

  auto real_size = realGrid.getNbas();
  for (int i = 0; i < real_size; ++i) {
    for (int j = 0; j < real_size; ++j) {
      std::cout << realGrid.getPoint(i) * realGrid.getPoint(j) << "  i" << i
                << "  j" << j << "\n";
    }
  }

  std::cout << "\n";
  std::cout << "\n";
  std::cout << "\n";
  std::cout << "\n";

  auto scatter_lobattoQuad = std::make_unique<Quadrature_Lobatto>(ndvr);
  auto radauQuad = std::make_unique<Quadrature_Radau>(mdvr);
  FEMDVR scatterGrid(std::move(scatter_lobattoQuad), std::move(radauQuad), nel);
  scatterGrid.print();

  auto size = scatterGrid.getNbas();
  for (int i = 0; i < size; ++i) {
    for (int j = 0; j < size; ++j) {
      std::cout << scatterGrid.getLaplacian(i * size + j) << "  i" << i << "  j"
                << j << "\n";
    }
  }

  return 0;
}
