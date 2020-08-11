
#include "../include/AngularGrid.hpp"
#include "../include/sphere_lebedev_rule.hpp"
#include <iostream>
#include <math.h>

AngularGrid::AngularGrid(const unsigned int &a_lmax,
                         const unsigned int &a_angular_order) {
  /* TODO: Using smaller lmax for testing like Roger */
  /* m_lmax = a_lmax; */
  m_lmax = floor(a_lmax / 2);
  m_angular_order = a_angular_order;

  /* Lebedev Quadrature */
  m_angular_x = std::make_unique<double[]>(m_angular_order);
  m_angular_y = std::make_unique<double[]>(m_angular_order);
  m_angular_z = std::make_unique<double[]>(m_angular_order);
  m_angular_weight = std::make_unique<double[]>(m_angular_order);
  ld2354(m_angular_x.get(), m_angular_y.get(), m_angular_z.get(),
         m_angular_weight.get());
  // the number in the function's name has to be the order passed!

  /* Generating lm pairs for grid lmax passed in */
  m_num_channels = (m_lmax - 1) * (m_lmax + 1) + 1;
  m_quantum_number_l = std::make_unique<int[]>(m_num_channels);
  m_quantum_number_m = std::make_unique<int[]>(m_num_channels);
  int k = 0;
  for (int i = 0; i < m_lmax; ++i) {
    for (int j = 0; j <= i; ++j) {
      if (j == 0) {
        m_quantum_number_l[k] = i;
        m_quantum_number_m[k] = j;
        k++;
      } else {
        m_quantum_number_l[k] = i;
        m_quantum_number_m[k] = j;
        k++;
        m_quantum_number_l[k] = i;
        m_quantum_number_m[k] = -j;
        k++;
      }
    }
  }
}

AngularGrid::~AngularGrid() {}

double AngularGrid::getAngularWeight(int index) const {
  return m_angular_weight[index];
}

double AngularGrid::getX(int index) const { return m_angular_x[index]; }

double AngularGrid::getY(int index) const { return m_angular_y[index]; }

double AngularGrid::getZ(int index) const { return m_angular_z[index]; }

unsigned int AngularGrid::getLmax() const { return m_lmax; }

unsigned int AngularGrid::getNumChannels() const { return m_num_channels; }

int AngularGrid::getL(int index) const { return m_quantum_number_l[index]; }

int AngularGrid::getM(int index) const { return m_quantum_number_m[index]; }

unsigned int AngularGrid::getAngularOrder() const { return m_angular_order; }
