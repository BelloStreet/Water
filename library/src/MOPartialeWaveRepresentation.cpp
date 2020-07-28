
#include <complex>
#include <gsl/gsl_sf_legendre.h>
#include <iostream>
#include <vector>

#include "../include/MOPartialWaveRepresentation.hpp"
#include "../include/sphere_lebedev_rule.hpp"

MOPartialWaveRepresentation::MOPartialWaveRepresentation(
    std::shared_ptr<FEMDVR> a_femdvr_grid,
    std::shared_ptr<AngularGrid> a_angular_grid, const molecule_t &a_molecule,
    const unsigned int &a_which_MO) {

  // Generating lm pairs for this irrep
  int l = (int)(ceil(a_angular_grid->getLmax() / 2));
  std::vector<int> l0_vec;
  std::vector<int> m0_vec;
  int m0, l0;
  std::string type;
  if (a_molecule.pointgroup != "NAN") {
    if (a_molecule.pointgroup == "C2v") {
      if (a_molecule.molecular_orbitals[a_which_MO].symmetry == "A1") {
        m0 = 0;
        type = "cosine";
      } else if (a_molecule.molecular_orbitals[a_which_MO].symmetry == "B1") {
        m0 = 1;
        type = "cosine";
      } else if (a_molecule.molecular_orbitals[a_which_MO].symmetry == "B2") {
        m0 = 1;
        type = "sine";
      } else if (a_molecule.molecular_orbitals[a_which_MO].symmetry == "A2") {
        m0 = 2;
        type = "sine";
      }
      for (int i = 0; i < l; ++i) {
        for (int m = m0; m <= i; m = m + 2) {
          l0_vec.push_back(i);
          m0_vec.push_back(m);
        }
      }
    } else if (a_molecule.pointgroup == "D2h") {
      if (a_molecule.molecular_orbitals[a_which_MO].symmetry == "Ag") {
        m0 = 0;
        l0 = 0;
        type = "cosine";
      } else if (a_molecule.molecular_orbitals[a_which_MO].symmetry == "B3u") {
        m0 = 1;
        l0 = 1;
        type = "cosine";
      } else if (a_molecule.molecular_orbitals[a_which_MO].symmetry == "B2u") {
        m0 = 1;
        l0 = 1;
        type = "sine";
      } else if (a_molecule.molecular_orbitals[a_which_MO].symmetry == "B1g") {
        m0 = 0;
        l0 = 0;
        type = "sine";
      } else if (a_molecule.molecular_orbitals[a_which_MO].symmetry == "B1u") {
        m0 = 0;
        l0 = 1;
        type = "cosine";
      } else if (a_molecule.molecular_orbitals[a_which_MO].symmetry == "B2g") {
        m0 = 1;
        l0 = 0;
        type = "cosine";
      } else if (a_molecule.molecular_orbitals[a_which_MO].symmetry == "B3g") {
        m0 = 1;
        l0 = 0;
        type = "sine";
      } else if (a_molecule.molecular_orbitals[a_which_MO].symmetry == "Au") {
        m0 = 0;
        l0 = 1;
        type = "sine";
      }
      for (int i = l0; i < l; i = i + 2) {
        for (int m = m0; m <= i; m = m + 2) {
          l0_vec.push_back(i);
          m0_vec.push_back(m);
        }
      }
    } else if (a_molecule.pointgroup == "Cs") {
      if (a_molecule.molecular_orbitals[a_which_MO].symmetry == "A'") {
      } else if (a_molecule.molecular_orbitals[a_which_MO].symmetry == "A''") {
      }
    }
  } else {
    if (a_molecule.molecular_orbitals[a_which_MO].symmetry == "A1") {
      m0 = 0;
      type = "cosine";
    } else if (a_molecule.molecular_orbitals[a_which_MO].symmetry == "B1") {
      m0 = 1;
      type = "cosine";
    } else if (a_molecule.molecular_orbitals[a_which_MO].symmetry == "B2") {
      m0 = 1;
      type = "sine";
    } else if (a_molecule.molecular_orbitals[a_which_MO].symmetry == "A2") {
      m0 = 2;
      type = "sine";
    }
    for (int i = 0; i < l; ++i) {
      for (int m = m0; m <= i; m = m + 2) {
        l0_vec.push_back(i);
        m0_vec.push_back(m);
      }
    }
  }
  m_num_channels = l0_vec.size();
  m_quantum_number_l = std::make_unique<int[]>(m_num_channels);
  m_quantum_number_m = std::make_unique<int[]>(m_num_channels);
  for (int i = 0; i < m_num_channels; ++i) {
    m_quantum_number_l[i] = l0_vec[i];
    m_quantum_number_m[i] = m0_vec[i];
  }

  auto X = std::make_unique<double[]>(a_angular_grid->getAngularOrder() *
                                      a_femdvr_grid->getNbas());
  auto Y = std::make_unique<double[]>(a_angular_grid->getAngularOrder() *
                                      a_femdvr_grid->getNbas());
  auto Z = std::make_unique<double[]>(a_angular_grid->getAngularOrder() *
                                      a_femdvr_grid->getNbas());
  for (int i = 0; i < a_femdvr_grid->getNbas(); ++i) {
    for (int j = 0; j < a_angular_grid->getAngularOrder(); ++j) {
      X[i * a_angular_grid->getAngularOrder() + j] =
          a_angular_grid->getX(j) * a_femdvr_grid->getPoint(i).real();
      Y[i * a_angular_grid->getAngularOrder() + j] =
          a_angular_grid->getY(j) * a_femdvr_grid->getPoint(i).real();
      Z[i * a_angular_grid->getAngularOrder() + j] =
          a_angular_grid->getZ(j) * a_femdvr_grid->getPoint(i).real();
    }
  }

  auto weight = std::make_unique<double[]>(4);
  for (int i = 0; i < 4; ++i) {
    weight[i] = 1.0;
    for (int j = 2 * i - 1; j > 0; j = j - 2) {
      weight[i] *= sqrt(1.0 / j);
    }
  }
  auto radial_basis = std::make_unique<double[]>(
      a_angular_grid->getAngularOrder() * a_femdvr_grid->getNbas());
  for (int xgrid = 0;
       xgrid < a_angular_grid->getAngularOrder() * a_femdvr_grid->getNbas();
       ++xgrid) {
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
            S_xyz * a_molecule.molecular_orbitals[a_which_MO].value[mo_index];
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
                 a_molecule.molecular_orbitals[a_which_MO].value[mo_index];
        mo_index++;
        value += S_xyz * y_diff *
                 a_molecule.molecular_orbitals[a_which_MO].value[mo_index];
        mo_index++;
        value += S_xyz * z_diff *
                 a_molecule.molecular_orbitals[a_which_MO].value[mo_index];
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
                 a_molecule.molecular_orbitals[a_which_MO].value[mo_index] *
                 weight[2];
        mo_index++;
        value += S_xyz * pow(y_diff, 2) *
                 a_molecule.molecular_orbitals[a_which_MO].value[mo_index] *
                 weight[2];
        mo_index++;
        value += S_xyz * pow(z_diff, 2) *
                 a_molecule.molecular_orbitals[a_which_MO].value[mo_index] *
                 weight[2];
        mo_index++;
        value += S_xyz * x_diff * y_diff *
                 a_molecule.molecular_orbitals[a_which_MO].value[mo_index];
        mo_index++;
        value += S_xyz * x_diff * z_diff *
                 a_molecule.molecular_orbitals[a_which_MO].value[mo_index];
        mo_index++;
        value += S_xyz * y_diff * z_diff *
                 a_molecule.molecular_orbitals[a_which_MO].value[mo_index];
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
                 a_molecule.molecular_orbitals[a_which_MO].value[mo_index] *
                 weight[3];
        mo_index++;
        value += S_xyz * pow(y_diff, 3) *
                 a_molecule.molecular_orbitals[a_which_MO].value[mo_index] *
                 weight[3];
        mo_index++;
        value += S_xyz * pow(z_diff, 3) *
                 a_molecule.molecular_orbitals[a_which_MO].value[mo_index] *
                 weight[3];
        mo_index++;
        value += S_xyz * pow(y_diff, 2) *
                 a_molecule.molecular_orbitals[a_which_MO].value[mo_index] *
                 weight[2];
        mo_index++;
        value += S_xyz * pow(x_diff, 2) *
                 a_molecule.molecular_orbitals[a_which_MO].value[mo_index] *
                 weight[2];
        mo_index++;
        value += S_xyz * pow(x_diff, 2) *
                 a_molecule.molecular_orbitals[a_which_MO].value[mo_index] *
                 weight[2];
        mo_index++;
        value += S_xyz * pow(z_diff, 2) *
                 a_molecule.molecular_orbitals[a_which_MO].value[mo_index] *
                 weight[2];
        mo_index++;
        value += S_xyz * pow(z_diff, 2) *
                 a_molecule.molecular_orbitals[a_which_MO].value[mo_index] *
                 weight[2];
        mo_index++;
        value += S_xyz * pow(y_diff, 2) *
                 a_molecule.molecular_orbitals[a_which_MO].value[mo_index] *
                 weight[2];
        mo_index++;
        value += S_xyz * x_diff * y_diff * z_diff *
                 a_molecule.molecular_orbitals[a_which_MO].value[mo_index];
        mo_index++;
      }
    }
    radial_basis[xgrid] = value;
  }

  m_blm = std::make_unique<double[]>(m_num_channels);
  m_partial_wave_rep_orbital = std::make_unique<std::complex<double>[]>(
      a_femdvr_grid->getNbas() * m_num_channels);
  for (int channel = 0; channel < m_num_channels; ++channel) {
    for (int radial_point = 0; radial_point < a_femdvr_grid->getNbas();
         ++radial_point) {
      std::complex<double> beta = 0.0;
      for (int angular_point = 0;
           angular_point < a_angular_grid->getAngularOrder(); ++angular_point) {
        double theta = acos(a_angular_grid->getZ(angular_point));
        double phi = atan2(a_angular_grid->getY(angular_point),
                           a_angular_grid->getX(angular_point));
        m_blm[channel] =
            Blm(type, l0_vec[channel], m0_vec[channel], theta, phi);
        beta += sqrt(a_femdvr_grid->getWeight(radial_point).real()) *
                a_angular_grid->getAngularWeight(angular_point) *
                radial_basis[angular_point +
                             a_angular_grid->getAngularOrder() * radial_point] *
                Blm(type, l0_vec[channel], m0_vec[channel], theta, phi) *
                a_femdvr_grid->getPoint(radial_point).real();
      }
      m_partial_wave_rep_orbital[channel * a_femdvr_grid->getNbas() +
                                 radial_point] = 4.0 * M_PI * beta;
    }
  }
}

MOPartialWaveRepresentation::~MOPartialWaveRepresentation() {}

std::complex<double>
MOPartialWaveRepresentation::getPartialWaveRep(int index) const {
  return m_partial_wave_rep_orbital[index];
}

double MOPartialWaveRepresentation::getBLM(int index) const {
  return m_blm[index];
}

unsigned int MOPartialWaveRepresentation::getNumChannels() const {
  return m_num_channels;
}

int MOPartialWaveRepresentation::getQuantumNumberL(int index) const {
  return m_quantum_number_l[index];
}

int MOPartialWaveRepresentation::getQuantumNumberM(int index) const {
  return m_quantum_number_m[index];
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
