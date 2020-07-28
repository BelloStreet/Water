#ifndef _GLOBAL_DATATYPES_HPP_
#define _GLOBAL_DATATYPES_HPP_
#include <string>
#include <vector>

struct Sorbitalparameters_t {
  std::vector<double> expon_s;
  std::vector<double> coeff_s;
};
struct Porbitalparameters_t {
  std::vector<double> expon_p;
  std::vector<double> coeff_p;
};
struct Dorbitalparameters_t {
  std::vector<double> expon_d;
  std::vector<double> coeff_d;
};
struct Forbitalparameters_t {
  std::vector<double> expon_f;
  std::vector<double> coeff_f;
};

struct AOorbital_t {
  std::vector<Sorbitalparameters_t> s_orbitals;
  std::vector<Porbitalparameters_t> p_orbitals;
  std::vector<Dorbitalparameters_t> d_orbitals;
  std::vector<Forbitalparameters_t> f_orbitals;
};

struct atom_t {
  AOorbital_t AOorbital;
  std::string symbol;
  double x;
  double y;
  double z;
  int charge;
};

struct MO_t {
  std::vector<double> value;
  double energy;
  double occupation;
  std::string symmetry;
  std::string spin;
};

struct molecule_t {
  std::vector<atom_t> atoms;
  std::vector<MO_t> molecular_orbitals;
  std::string pointgroup;
};

#endif
