#include "../include/Quadrature_Radau.hpp"

#include <complex>
#include <iostream>

#include "../include/Quadrature.hpp"

Quadrature_Radau::Quadrature_Radau(int &a_dvr) : Quadrature(a_dvr) {
  m_dvr = a_dvr;
  m_kind = 6, m_kpts = 1, m_alpha = 0, m_beta = 0;
  auto x = std::make_unique<double[]>(a_dvr);
  auto w = std::make_unique<double[]>(a_dvr);
  auto src = std::make_unique<double[]>(a_dvr);
  auto endpts = std::make_unique<double[]>(2);
  endpts[0] = 0.0;
  endpts[1] = 1.0;
  gaussq_(&m_kind, &a_dvr, &m_alpha, &m_beta, &m_kpts, endpts.get(), src.get(),
          x.get(), w.get());
  m_points = std::make_unique<std::complex<double>[]>(a_dvr);
  m_weights = std::make_unique<std::complex<double>[]>(a_dvr);
  for (int i = 0; i < a_dvr; ++i) {
    m_points[i] = x[i];
    m_weights[i] = w[i];
  }
}

Quadrature_Radau::~Quadrature_Radau(){};

inline std::complex<double> Quadrature_Radau::getPoint(int index) const {
  return m_points[index];
}

inline std::complex<double> Quadrature_Radau::getWeight(int index) const {
  return m_weights[index];
}

void Quadrature_Radau::print() const {
  std::cout << "Radau Quadrature" << std::endl;
  for (int i = 0; i < m_dvr; ++i) {
    std::cout << "Points " << m_points[i] << " "
              << "Weights " << m_weights[i] << std::endl;
  }
  std::cout << "\n";
}
