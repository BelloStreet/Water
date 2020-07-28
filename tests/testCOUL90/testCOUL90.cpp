// cSpell:disable
#include <gnuplot-iostream.h>

#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "COUL90.hpp"
#include "FEMDVR.hpp"

int main(int argc, char **argv) {
  // Exception pratice
  std::string arg = argv[1];
  int plot;
  try {
    std::size_t pos;
    plot = std::stoi(arg, &pos);
    if (plot < 0 || plot > 1) {
      throw std::runtime_error(
          "Input must be 1 for plotting or 0 for turning plotting off.");
    }
  } catch (std::exception &exp) {
    std::cerr << "Error: " << exp.what() << '\n';
    return 0;
  }

  int ndvr = 10, mdvr = 10;
  unsigned int nel = 4;
  double theta = 30.0, alphaRad = 0.3;

  auto lobattoQuad = std::make_shared<Quadrature_Lobatto>(ndvr);
  auto radauQuad = std::make_shared<Quadrature_Radau>(mdvr);
  auto grid =
      std::make_shared<FEMDVR>(lobattoQuad, radauQuad, nel, theta, alphaRad);

  int charge = 2;
  double k = 1;
  COUL90 coulomb(grid, charge, k);
  auto points = grid->getPoints();
  auto wave = coulomb.dvrRep();

  /* The gnuplot interface requires the use of container classes */
  std::vector<double> x, y;
  x.reserve(grid->getNbas());
  y.reserve(grid->getNbas());
  for (int i = 0; i < grid->getNbas(); ++i) {
    x.push_back((points + i)->real());
    y.push_back((wave + i)->real());
  }

  if (plot) {
    Gnuplot gp;
    // '-' means read from stdin.  The send1d() function sends data to
    // gnuplot's stdin.
    /* gp << "set yrange [-1:1]\n"; */
    gp << "plot '-' with lines \n";
    gp.send1d(boost::make_tuple(x, y));
  }

  return 0;
}
