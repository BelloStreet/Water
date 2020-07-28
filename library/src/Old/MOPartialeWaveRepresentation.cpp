
#include <complex>
#include <gsl/gsl_sf_legendre.h>
#include <iostream>
#include <vector>

#include "../include/MOPartialWaveRepresentation.hpp"
#include "../include/sphere_lebedev_rule.hpp"

MOPartialWaveRepresentation::MOPartialWaveRepresentation(
    std::shared_ptr<FEMDVR> a_grid, const molecule_t &a_molecule,
    const std::string &a_irrep, const unsigned int &a_lmax,
    const int &a_angular_order) {
  m_lmax = a_lmax;
  m_angular_order = a_angular_order;

  // Generating lm pairs
  int l = (int)(ceil(m_lmax / 2));
  std::vector<int> l0_list;
  std::vector<int> m0_list;
  int m0, l0, which_MO;
  std::string type;
  if (a_molecule.pointgroup != "NAN") {
    if (a_molecule.pointgroup == "C2v") {
      if (a_irrep == "A1") {
        m0 = 0;
        type = "cosine";
        which_MO = 0;
      } else if (a_irrep == "B1") {
        m0 = 2;
        type = "cosine";
        which_MO = 1;
      } else if (a_irrep == "B2") {
        m0 = 1;
        type = "sine";
        which_MO = 2;
      } else if (a_irrep == "A2") {
        m0 = 1;
        type = "sine";
        which_MO = 3;
      }
      for (int i = 0; i < l; ++i) {
        for (int m = m0; m <= i; m = m + 2) {
          l0_list.push_back(i);
          m0_list.push_back(m);
        }
      }
    } else if (a_molecule.pointgroup == "D2h") {
      if (a_irrep == "Ag") {
        m0 = 0;
        l0 = 0;
        type = "cosine";
        which_MO = 0;
      } else if (a_irrep == "B3u") {
        m0 = 1;
        l0 = 1;
        type = "cosine";
        which_MO = 1;
      } else if (a_irrep == "B2u") {
        m0 = 1;
        l0 = 1;
        type = "sine";
        which_MO = 2;
      } else if (a_irrep == "B1g") {
        m0 = 0;
        l0 = 0;
        type = "cosine";
        which_MO = 3;
      } else if (a_irrep == "B1u") {
        m0 = 0;
        l0 = 1;
        type = "cosine";
        which_MO = 4;
      } else if (a_irrep == "B2g") {
        m0 = 1;
        l0 = 0;
        type = "sine";
        which_MO = 5;
      } else if (a_irrep == "B3g") {
        m0 = 1;
        l0 = 0;
        type = "sine";
        which_MO = 6;
      } else if (a_irrep == "Au") {
        m0 = 0;
        l0 = 1;
        type = "sine";
        which_MO = 7;
      }
      for (int i = l0; i < l; i = i + 2) {
        for (int m = m0; m <= i; m = m + 2) {
          l0_list.push_back(i);
          m0_list.push_back(m);
        }
      }
    } else if (a_molecule.pointgroup == "Cs") {
      if (a_irrep == "A'") {
      } else if (a_irrep == "A''") {
      }
    }
  } else {
    if (a_irrep == "A1") {
      m0 = 0;
      type = "cosine";
      which_MO = 0;
    } else if (a_irrep == "B1") {
      m0 = 2;
      type = "cosine";
      which_MO = 1;
    } else if (a_irrep == "B2") {
      m0 = 1;
      type = "sine";
      which_MO = 2;
    } else if (a_irrep == "A2") {
      m0 = 1;
      type = "sine";
      which_MO = 3;
    }
    for (int i = 0; i < l; ++i) {
      for (int m = m0; m <= i; m = m + 2) {
        l0_list.push_back(i);
        m0_list.push_back(m);
      }
    }
  }

  /* Lebedev Quadrature */
  auto x_ang = std::make_unique<double[]>(a_angular_order);
  auto y_ang = std::make_unique<double[]>(a_angular_order);
  auto z_ang = std::make_unique<double[]>(a_angular_order);
  auto w_ang = std::make_unique<double[]>(a_angular_order);
  ld2354(x_ang.get(), y_ang.get(), z_ang.get(), w_ang.get());
  // the number in the function's name has to be the order passed!

  auto X = std::make_unique<double[]>(a_angular_order * a_grid->getNbas());
  auto Y = std::make_unique<double[]>(a_angular_order * a_grid->getNbas());
  auto Z = std::make_unique<double[]>(a_angular_order * a_grid->getNbas());
  for (int i = 0; i < a_grid->getNbas(); ++i) {
    for (int j = 0; j < a_angular_order; ++j) {
      X[i * a_angular_order + j] = x_ang[j] * real(a_grid->getPoints()[i]);
      Y[i * a_angular_order + j] = y_ang[j] * real(a_grid->getPoints()[i]);
      Z[i * a_angular_order + j] = z_ang[j] * real(a_grid->getPoints()[i]);
    }
  }

  auto weight = std::make_unique<double[]>(4);
  for (int i = 0; i < 4; ++i) {
    weight[i] = 1.0;
    for (int j = 2 * i - 1; j > 0; j = j - 2) {
      weight[i] *= sqrt(1.0 / j);
    }
  }
  auto radial_basis =
      std::make_unique<double[]>(a_angular_order * a_grid->getNbas());
  for (int xgrid = 0; xgrid < a_angular_order * a_grid->getNbas(); ++xgrid) {
    double value = 0.0;
    int mo_index = 0;
    for (auto atom : a_molecule.atoms) {
      double x_diff = X[xgrid] - atom.x;
      double y_diff = Y[xgrid] - atom.y;
      double z_diff = Z[xgrid] - atom.z;
      for (auto s_orbital : atom.AOorbital.s_orbitals) {
        double S_xyz = 0.0;
        for (int i = 0; i < s_orbital.coeff_s.size(); ++i) {
          S_xyz += s_orbital.coeff_s[i] *
                   exp(-s_orbital.expon_s[i] *
                       (x_diff * x_diff + y_diff * y_diff + z_diff * z_diff));
        }
        value +=
            S_xyz * a_molecule.molecular_orbitals[which_MO].value[mo_index];
        mo_index++;
      }
      for (auto p_orbital : atom.AOorbital.p_orbitals) {
        double S_xyz = 0.0;
        for (int i = 0; i < p_orbital.coeff_p.size(); ++i) {
          S_xyz += p_orbital.coeff_p[i] *
                   exp(-p_orbital.expon_p[i] *
                       (x_diff * x_diff + y_diff * y_diff + z_diff * z_diff));
        }
        value += S_xyz * x_diff *
                 a_molecule.molecular_orbitals[which_MO].value[mo_index];
        mo_index++;
        value += S_xyz * y_diff *
                 a_molecule.molecular_orbitals[which_MO].value[mo_index];
        mo_index++;
        value += S_xyz * z_diff *
                 a_molecule.molecular_orbitals[which_MO].value[mo_index];
        mo_index++;
      }
      for (auto d_orbital : atom.AOorbital.d_orbitals) {
        double S_xyz = 0.0;
        for (int i = 0; i < d_orbital.coeff_d.size(); ++i) {
          S_xyz += d_orbital.coeff_d[i] *
                   exp(-d_orbital.expon_d[i] *
                       (x_diff * x_diff + y_diff * y_diff + z_diff * z_diff));
        }
        value += S_xyz * pow(x_diff, 2) *
                 a_molecule.molecular_orbitals[which_MO].value[mo_index] *
                 weight[2];
        mo_index++;
        value += S_xyz * pow(y_diff, 2) *
                 a_molecule.molecular_orbitals[which_MO].value[mo_index] *
                 weight[2];
        mo_index++;
        value += S_xyz * pow(z_diff, 2) *
                 a_molecule.molecular_orbitals[which_MO].value[mo_index] *
                 weight[2];
        mo_index++;
        value += S_xyz * x_diff * y_diff *
                 a_molecule.molecular_orbitals[which_MO].value[mo_index];
        mo_index++;
        value += S_xyz * x_diff * z_diff *
                 a_molecule.molecular_orbitals[which_MO].value[mo_index];
        mo_index++;
        value += S_xyz * y_diff * z_diff *
                 a_molecule.molecular_orbitals[which_MO].value[mo_index];
        mo_index++;
      }
      for (auto f_orbital : atom.AOorbital.f_orbitals) {
        double S_xyz = 0.0;
        for (int i = 0; i < f_orbital.coeff_f.size(); ++i) {
          S_xyz += f_orbital.coeff_f[i] *
                   exp(-f_orbital.expon_f[i] *
                       (x_diff * x_diff + y_diff * y_diff + z_diff * z_diff));
        }
        value += S_xyz * pow(x_diff, 3) *
                 a_molecule.molecular_orbitals[which_MO].value[mo_index] *
                 weight[3];
        mo_index++;
        value += S_xyz * pow(y_diff, 3) *
                 a_molecule.molecular_orbitals[which_MO].value[mo_index] *
                 weight[3];
        mo_index++;
        value += S_xyz * pow(z_diff, 3) *
                 a_molecule.molecular_orbitals[which_MO].value[mo_index] *
                 weight[3];
        mo_index++;
        value += S_xyz * pow(y_diff, 2) *
                 a_molecule.molecular_orbitals[which_MO].value[mo_index] *
                 weight[2];
        mo_index++;
        value += S_xyz * pow(x_diff, 2) *
                 a_molecule.molecular_orbitals[which_MO].value[mo_index] *
                 weight[2];
        mo_index++;
        value += S_xyz * pow(x_diff, 2) *
                 a_molecule.molecular_orbitals[which_MO].value[mo_index] *
                 weight[2];
        mo_index++;
        value += S_xyz * pow(z_diff, 2) *
                 a_molecule.molecular_orbitals[which_MO].value[mo_index] *
                 weight[2];
        mo_index++;
        value += S_xyz * pow(z_diff, 2) *
                 a_molecule.molecular_orbitals[which_MO].value[mo_index] *
                 weight[2];
        mo_index++;
        value += S_xyz * pow(y_diff, 2) *
                 a_molecule.molecular_orbitals[which_MO].value[mo_index] *
                 weight[2];
        mo_index++;
        value += S_xyz * x_diff * y_diff * z_diff *
                 a_molecule.molecular_orbitals[which_MO].value[mo_index];
        mo_index++;
      }
    }
    radial_basis[xgrid] = value;
  }

  m_num_channels = l0_list.size();
  m_partial_wave_rep_orbital = std::make_unique<std::complex<double>[]>(
      a_grid->getNbas() * m_num_channels);
  for (int channel = 0; channel < m_num_channels; ++channel) {
    for (int radial_point = 0; radial_point < a_grid->getNbas();
         ++radial_point) {
      std::complex<double> beta = 0.0;
      for (int angular_point = 0; angular_point < a_angular_order;
           ++angular_point) {
        double theta = acos(z_ang[angular_point]);
        double phi = atan2(y_ang[angular_point], x_ang[angular_point]);
        beta += sqrt(real(a_grid->getWeights()[radial_point])) *
                w_ang[angular_point] *
                radial_basis[angular_point + a_angular_order * radial_point] *
                Blm(type, l0_list[channel], m0_list[channel], theta, phi) *
                real(a_grid->getPoints()[radial_point]);
      }
      m_partial_wave_rep_orbital[channel * a_grid->getNbas() + radial_point] =
          4.0 * M_PI * beta;
    }
  }
}

MOPartialWaveRepresentation::~MOPartialWaveRepresentation() {}

std::complex<double> *MOPartialWaveRepresentation::getPartialWaveRep() const {
  return m_partial_wave_rep_orbital.get();
}

unsigned int MOPartialWaveRepresentation::getNumChannels() const {
  return m_num_channels;
}

double MOPartialWaveRepresentation::Blm(const std::string type, const int l,
                                        const int m, const double theta,
                                        const double phi) {
  double value = 0.0, y = 0.0;
  int m1 = abs(m);
  if (m != 0) {
    y = sqrt(2.0) * gsl_sf_legendre_sphPlm(l, m1, cos(theta));
  } else {
    y = gsl_sf_legendre_sphPlm(l, m1, cos(theta));
  }
  if (type == "cosine") {
    value = y * cos(m1 * phi);
  } else if (type == "sine") {
    value = y * sin(m1 * phi);
  }
  return value;
}
