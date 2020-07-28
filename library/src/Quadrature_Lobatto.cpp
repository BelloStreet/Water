#include "../include/Quadrature_Lobatto.hpp"

#include <iostream>

#include "../include/Quadrature.hpp"

Quadrature_Lobatto::Quadrature_Lobatto(int &a_dvr) : Quadrature(a_dvr) {
  m_kind = 1, m_kpts = 2, m_alpha = 0, m_beta = 0;
  m_dvr = a_dvr;
  auto x = std::make_unique<double[]>(a_dvr);
  auto w = std::make_unique<double[]>(a_dvr);
  auto src = std::make_unique<double[]>(a_dvr);
  auto endpts = std::make_unique<double[]>(2);
  endpts[0] = -1.0;
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

Quadrature_Lobatto::~Quadrature_Lobatto(){};

inline std::complex<double> Quadrature_Lobatto::getPoint(int index) const {
  return m_points[index];
}

inline std::complex<double> Quadrature_Lobatto::getWeight(int index) const {
  return m_weights[index];
}

void Quadrature_Lobatto::print() const {
  std::cout << "Lobatto Quadrature" << std::endl;
  for (int i = 0; i < m_dvr; ++i) {
    std::cout << "Points " << m_points[i] << " "
              << "Weights " << m_weights[i] << std::endl;
  }
}
