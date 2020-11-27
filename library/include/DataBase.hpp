#pragma once

#include "Joperator.hpp"
#include "Koperator.hpp"
#include "MOPartialWaveRepresentation.hpp"

#include <stdexcept>
#include <vector>

struct OrbitalDataBase {
  OrbitalDataBase(const FEMDVR &a_radial_grid,
                  const AngularGrid &a_angular_grid,
                  const uint32_t &a_number_orbitals)
      : m_nbas(a_radial_grid.getNbas()) {
    m_orbital_coefs.reserve(a_number_orbitals * a_radial_grid.getNbas() *
                            a_angular_grid.getNumChannels());
    m_orbital_coefs_index.reserve(a_number_orbitals * a_radial_grid.getNbas() *
                                  a_angular_grid.getNumChannels());
    m_orbital_l_quantum_number.reserve(a_number_orbitals *
                                       a_angular_grid.getNumChannels());
    m_orbital_m_quantum_number.reserve(a_number_orbitals *
                                       a_angular_grid.getNumChannels());
    m_quantum_number_index.reserve(a_number_orbitals *
                                   a_angular_grid.getNumChannels());
    m_orbital_type.reserve(a_number_orbitals);
  }

  void StoreData(const uint32_t &a_ith_orbital,
                 const MOPartialWaveRepresentation &a_MO) {

    m_number_of_lm_pairs_per_orbital.push_back(a_MO.getNumChannels());
    m_orbital_type.push_back(a_MO.getType());

    m_orbital_coefs.reserve(a_ith_orbital * m_nbas * a_MO.getNumChannels());
    m_orbital_l_quantum_number.reserve(a_ith_orbital * a_MO.getNumChannels());
    m_orbital_m_quantum_number.reserve(a_ith_orbital * a_MO.getNumChannels());

    for (uint32_t i = 0; i < a_MO.getNumChannels(); ++i) {
      m_orbital_l_quantum_number.push_back(a_MO.getL(i));
      m_orbital_m_quantum_number.push_back(a_MO.getM(i));
      m_quantum_number_index.push_back(a_ith_orbital * a_MO.getNumChannels() +
                                       i);
    }
    for (uint32_t i = 0; i < m_nbas * a_MO.getNumChannels(); ++i) {
      m_orbital_coefs.push_back(a_MO.getPartialWaveRep(i));
      m_orbital_coefs_index.push_back(
          a_ith_orbital * m_nbas * m_number_of_lm_pairs_per_orbital.back() + i);
    }
  }

  std::complex<double> accessCoefficient(const uint32_t &ith_orbital,
                                         const uint32_t &i) const {
    /* Assert(i < orbital_data.m_orbital_coefs.size(), */
    /*        ExcMessage("Index is out of coef bounds!")); */
    uint32_t coef_index = m_orbital_coefs_index
        [ith_orbital * m_nbas * m_number_of_lm_pairs_per_orbital[ith_orbital] +
         i];
    return m_orbital_coefs[coef_index];
  }

  std::string accessType(const uint32_t &ith_orbital) const {
    /* Assert(i < orbital_data.m_orbital_coefs.size(), */
    /*        ExcMessage("Index is out of coef bounds!")); */
    return m_orbital_type[ith_orbital];
  }

  uint32_t accessQuantumNumber(const uint32_t &ith_orbital, const uint32_t &i,
                               const std::string which_quantum_number) const {
    /* Assert(i < orbital_data.m_orbital_coefs.size(), */
    /*        ExcMessage("Index is out of coef bounds!")); */
    uint32_t index = m_quantum_number_index
        [ith_orbital * m_number_of_lm_pairs_per_orbital[ith_orbital] + i];
    if (which_quantum_number == "l")
      return m_orbital_l_quantum_number[index];
    else if (which_quantum_number == "m")
      return m_orbital_m_quantum_number[index];
    else
      throw std::invalid_argument(
          "received a character other than 'l' or 'm' for third argument");
  }

  /* private: */
  uint32_t m_nbas;
  std::vector<uint32_t> m_number_of_lm_pairs_per_orbital;
  std::vector<std::complex<double>> m_orbital_coefs;
  std::vector<int> m_orbital_coefs_index, m_orbital_l_quantum_number,
      m_orbital_m_quantum_number, m_quantum_number_index;
  std::vector<std::string> m_orbital_type;
};

struct OperatorDataBase {
  OperatorDataBase(const FEMDVR &a_radial_grid,
                   const AngularGrid &a_angular_grid,
                   const uint32_t &a_number_orbitals)
      : m_nbas(a_radial_grid.getNbas()) {
    m_J.reserve(a_number_orbitals * a_radial_grid.getNbas() *
                a_radial_grid.getNbas() * a_angular_grid.getNumChannels() *
                a_angular_grid.getNumChannels());
    m_K.reserve(a_number_orbitals * a_radial_grid.getNbas() *
                a_angular_grid.getNumChannels() * a_radial_grid.getNbas() *
                a_angular_grid.getNumChannels());
    m_operator_index.reserve(a_number_orbitals * a_radial_grid.getNbas() *
                             a_angular_grid.getNumChannels() *
                             a_radial_grid.getNbas() *
                             a_angular_grid.getNumChannels());
    m_number_channels.reserve(a_number_orbitals);
  }

  void StoreData(const uint32_t &a_ith_orbital, const Joperator &a_J,
                 const Koperator &a_K) {

    m_number_channels.push_back(a_J.getNumChannels());
    m_J.reserve(a_J.JSize());
    m_K.reserve(a_K.KSize());
    m_operator_index.reserve(a_ith_orbital * a_J.JSize());

    for (uint32_t i = 0; i < a_J.JSize(); ++i) {
      m_J.push_back(a_J.getJ(i));
      m_K.push_back(a_K.getK(i));
      m_operator_index.push_back(a_ith_orbital * a_J.JSize() + i);
    }
  }

  std::complex<double> accessJ(const uint32_t &ith_orbital,
                               const uint32_t &i) const {
    /* Assert(i < orbital_data.m_orbital_coefs.size(), */
    /*        ExcMessage("Index is out of coef bounds!")); */
    uint32_t J_index = m_operator_index[ith_orbital * m_nbas * m_nbas *
                                            m_number_channels[ith_orbital] *
                                            m_number_channels[ith_orbital] +
                                        i];
    return m_J[J_index];
  }

  std::complex<double> accessK(const uint32_t &ith_orbital,
                               const uint32_t &i) const {
    /* Assert(i < orbital_data.m_orbital_coefs.size(), */
    /*        ExcMessage("Index is out of coef bounds!")); */
    uint32_t K_index = m_operator_index[ith_orbital * m_nbas * m_nbas *
                                            m_number_channels[ith_orbital] *
                                            m_number_channels[ith_orbital] +
                                        i];
    return m_K[K_index];
  }

  uint32_t accessNumChannels(const uint32_t &ith_orbital) const {
    /* Assert(i < orbital_data.m_orbital_coefs.size(), */
    /*        ExcMessage("Index is out of coef bounds!")); */
    return m_number_channels[ith_orbital];
  }

  /* private: */
  uint32_t m_nbas;
  std::vector<std::complex<double>> m_J, m_K;
  std::vector<int> m_operator_index;
  std::vector<uint32_t> m_number_channels;
};
