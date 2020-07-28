// cSpell:disable
#include "../include/FEMDVR.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

/* TODO:Should use the geometry and Z charge of the system in the interface to
 * optimize the grid for that particular system */

FEMDVR::FEMDVR(std::unique_ptr<Quadrature_Lobatto> a_lobattoQuad,
               const unsigned int &a_Nelem) {
  m_Nelem = a_Nelem;
  m_Nbas = a_Nelem * (a_lobattoQuad->getDVROrder() - 1) - 1;
  std::complex<double> I(1, 0);

  /* TODO:Once the system's characteristics is passed in, then can use that to
   * adapt the grid */

  /* TODO: Redudancy for the boundaries but these are here only for testing */

  /* Testing with Roger have to pass 13 elements for this test*/
  /* m_realbounds.reserve(14); */
  /* m_realbounds.push_back(0.0); */
  /* m_realbounds.push_back(0.04); */
  /* m_realbounds.push_back(0.1); */
  /* m_realbounds.push_back(0.4); */
  /* m_realbounds.push_back(0.75); */
  /* m_realbounds.push_back(1.0); */
  /* m_realbounds.push_back(1.3); */
  /* m_realbounds.push_back(1.7782787998); */
  /* m_realbounds.push_back(2.0); */
  /* m_realbounds.push_back(4.0); */
  /* m_realbounds.push_back(8.0); */
  /* m_realbounds.push_back(12.0); */
  /* m_realbounds.push_back(16.0); */
  /* m_realbounds.push_back(20.0); */
  // Default finite element spacing
  m_realbounds.reserve(a_Nelem + 1);
  m_complexbounds.reserve(a_Nelem + 1);
  for (int i = 0; i < a_Nelem + 1; ++i) {
    m_realbounds.push_back(i * 10);
    /* m_realbounds.push_back(i * 0.05); */
  }
  /* TODO:Get rid of this... */
  m_complexbounds = m_realbounds;

  m_points = std::make_unique<std::complex<double>[]>(m_Nbas);
  m_weights = std::make_unique<std::complex<double>[]>(m_Nbas);

  auto lobbatto_order = a_lobattoQuad->getDVROrder();
  for (int element = 0; element < m_Nelem - 1; ++element) {
    auto lobtempX = std::make_unique<std::complex<double>[]>(lobbatto_order);
    auto lobtempW = std::make_unique<std::complex<double>[]>(lobbatto_order);
    for (int dvr = 0; dvr < lobbatto_order; ++dvr) {
      lobtempX[dvr] = m_realbounds[element] +
                      0.5 * (a_lobattoQuad->getPoint(dvr) + I) *
                          (m_realbounds[element + 1] - m_realbounds[element]);
      lobtempW[dvr] = 0.5 * a_lobattoQuad->getWeight(dvr) *
                      (m_realbounds[element + 1] - m_realbounds[element]);
    }
    if (element != 0) {
      m_weights[element * lobbatto_order - element - 1] += lobtempW[0];
    }
    for (int dvr = 1; dvr < lobbatto_order; ++dvr) {
      m_points[element * lobbatto_order + dvr - element - 1] = lobtempX[dvr];
      m_weights[element * lobbatto_order + dvr - element - 1] = lobtempW[dvr];
    }
  }
  auto lobtempX_last = std::make_unique<std::complex<double>[]>(lobbatto_order);
  auto lobtempW_last = std::make_unique<std::complex<double>[]>(lobbatto_order);
  for (int dvr = 0; dvr < lobbatto_order; ++dvr) {
    lobtempX_last[dvr] =
        m_realbounds[m_Nelem - 1] +
        0.5 * (a_lobattoQuad->getPoint(dvr) + I) *
            (m_realbounds[m_Nelem] - m_realbounds[m_Nelem - 1]);
    lobtempW_last[dvr] = 0.5 * a_lobattoQuad->getWeight(dvr) *
                         (m_realbounds[m_Nelem] - m_realbounds[m_Nelem - 1]);
  }
  m_weights[(m_Nelem - 1) * lobbatto_order - m_Nelem] += lobtempW_last[0];
  for (int dvr = 1; dvr < lobbatto_order - 1; ++dvr) {
    m_points[(m_Nelem - 1) * lobbatto_order + dvr - m_Nelem] =
        lobtempX_last[dvr];
    m_weights[(m_Nelem - 1) * lobbatto_order + dvr - m_Nelem] =
        lobtempW_last[dvr];
  }

  auto lobbatto_derivatives =
      std::make_unique<std::complex<double>[]>(lobbatto_order * lobbatto_order);
  for (int i = 0; i < lobbatto_order; ++i) {
    lobbatto_derivatives[i * lobbatto_order + i] = 0.0;
    for (int k = 0; k < lobbatto_order; ++k) {
      if (i != k)
        lobbatto_derivatives[i * lobbatto_order + i] +=
            1.0 / (a_lobattoQuad->getPoint(i) - a_lobattoQuad->getPoint(k));
    }
    for (int j = 0; j < lobbatto_order; ++j) {
      if (i != j) {
        std::complex<double> tmp =
            1.0 / (a_lobattoQuad->getPoint(j) - a_lobattoQuad->getPoint(i));
        for (int k = 0; k < lobbatto_order; ++k) {
          if ((k != i) && (k != j))
            tmp *= (a_lobattoQuad->getPoint(i) - a_lobattoQuad->getPoint(k)) /
                   (a_lobattoQuad->getPoint(j) - a_lobattoQuad->getPoint(k));
        }
        lobbatto_derivatives[i * lobbatto_order + j] = tmp;
      }
    }
  }
  m_laplacian = std::make_unique<std::complex<double>[]>(m_Nbas * m_Nbas);
  for (int element = 0; element < a_Nelem; ++element) {
    int start = 0;
    int end = lobbatto_order;
    if (element == 0)
      start = 1; // omit first DVR point
    if (element == a_Nelem - 1)
      end = lobbatto_order - 1; // omit last DVR point
    int offset = (lobbatto_order - 1) * element;
    for (int i = start; i < end; ++i) {
      int ith_index = (i - 1) + offset;
      for (int j = start; j < end; ++j) {
        int jth_index = (j - 1) + offset;
        std::complex<double> tmp = 0.0;
        for (int m = 0; m < lobbatto_order; ++m) {
          tmp += lobbatto_derivatives[m * lobbatto_order + j] *
                 lobbatto_derivatives[m * lobbatto_order + i] *
                 a_lobattoQuad->getWeight(m);
        }
        m_laplacian[ith_index * m_Nbas + jth_index] +=
            tmp * 2.0 /
            ((m_realbounds[element + 1] - m_realbounds[element]) *
             sqrt(m_weights[ith_index] * m_weights[jth_index]));
      }
    }
  }

} // end of constructor

FEMDVR::FEMDVR(std::unique_ptr<Quadrature_Lobatto> a_lobattoQuad,
               std::unique_ptr<Quadrature_Radau> a_radauQuad,
               const unsigned int &a_Nelem) {
  m_Nelem = a_Nelem;
  double theta = 30.0;
  double alphaRad = 0.3;
  double NRealbas = (a_Nelem - 1) * (a_lobattoQuad->getDVROrder() - 1) - 1;
  m_Nbas = (a_Nelem - 1) * (a_lobattoQuad->getDVROrder() - 1) +
           a_radauQuad->getDVROrder() - 1;
  std::complex<double> I(1, 0);

  /* TODO:Once the system's characteristics is passed in, then can use that to
   * adapt the grid */

  // Default finite element spacing
  m_realbounds.reserve(a_Nelem + 1);
  m_complexbounds.reserve(a_Nelem + 1);
  for (int i = 0; i < a_Nelem + 1; ++i) {
    m_realbounds.push_back(i * 10);
  }

  m_points = std::make_unique<std::complex<double>[]>(m_Nbas);
  m_weights = std::make_unique<std::complex<double>[]>(m_Nbas);

  auto lobbatto_order = a_lobattoQuad->getDVROrder();
  for (int element = 0; element < m_Nelem - 1; ++element) {
    auto lobtempX = std::make_unique<std::complex<double>[]>(lobbatto_order);
    auto lobtempW = std::make_unique<std::complex<double>[]>(lobbatto_order);
    for (int dvr = 0; dvr < lobbatto_order; ++dvr) {
      lobtempX[dvr] = m_realbounds[element] +
                      0.5 * (a_lobattoQuad->getPoint(dvr) + I) *
                          (m_realbounds[element + 1] - m_realbounds[element]);
      lobtempW[dvr] = 0.5 * a_lobattoQuad->getWeight(dvr) *
                      (m_realbounds[element + 1] - m_realbounds[element]);
    }
    if (element != 0) {
      m_weights[element * lobbatto_order - element - 1] += lobtempW[0];
    }
    for (int dvr = 1; dvr < lobbatto_order; ++dvr) {
      m_points[element * lobbatto_order + dvr - element - 1] = lobtempX[dvr];
      m_weights[element * lobbatto_order + dvr - element - 1] = lobtempW[dvr];
    }
  }

  m_complexbounds = m_realbounds;
  std::complex<double> R0 = m_realbounds[a_Nelem - 1];
  std::complex<double> eit = std::polar(1.0, M_PI * theta / 180.0);
  m_complexbounds[a_Nelem] = R0 + eit * (m_realbounds[a_Nelem] - R0);
  auto radau_order = a_radauQuad->getDVROrder();
  auto radtempX = std::make_unique<std::complex<double>[]>(radau_order);
  auto radtempW = std::make_unique<std::complex<double>[]>(radau_order);
  std::complex<double> radscale = eit / (2.0 * alphaRad);
  for (int i = 0; i < radau_order; ++i) {
    radtempX[i] = radscale * a_radauQuad->getPoint(i) + R0;
    radtempW[i] = radscale * a_radauQuad->getWeight(i);
  }
  m_weights[(m_Nelem - 1) * (radau_order - fabs(radau_order - lobbatto_order)) -
            m_Nelem] += radtempW[0];
  for (int dvr = 1; dvr < radau_order; ++dvr) {
    m_points[(m_Nelem - 1) *
                 (radau_order - fabs(radau_order - lobbatto_order)) +
             dvr - m_Nelem] = radtempX[dvr];
    m_weights[(m_Nelem - 1) *
                  (radau_order - fabs(radau_order - lobbatto_order)) +
              dvr - m_Nelem] = radtempW[dvr];
  }

  m_complexbounds[m_Nelem] = radtempX[radau_order - 1];
  m_realbounds[m_Nelem] = (m_complexbounds[m_Nelem] - R0) * conj(eit) + R0;

  // Laplacian
  auto lobbatto_derivatives =
      std::make_unique<std::complex<double>[]>(lobbatto_order * lobbatto_order);
  for (int i = 0; i < lobbatto_order; ++i) {
    lobbatto_derivatives[i * lobbatto_order + i] = 0.0;
    for (int k = 0; k < lobbatto_order; ++k) {
      if (i != k)
        lobbatto_derivatives[i * lobbatto_order + i] +=
            1.0 / (a_lobattoQuad->getPoint(i) - a_lobattoQuad->getPoint(k));
    }
    for (int j = 0; j < lobbatto_order; ++j) {
      if (i != j) {
        std::complex<double> tmp =
            1.0 / (a_lobattoQuad->getPoint(j) - a_lobattoQuad->getPoint(i));
        for (int k = 0; k < lobbatto_order; ++k) {
          if ((k != i) && (k != j))
            tmp *= (a_lobattoQuad->getPoint(i) - a_lobattoQuad->getPoint(k)) /
                   (a_lobattoQuad->getPoint(j) - a_lobattoQuad->getPoint(k));
        }
        lobbatto_derivatives[i * lobbatto_order + j] = tmp;
      }
    }
  }

  m_laplacian = std::make_unique<std::complex<double>[]>(m_Nbas * m_Nbas);
  for (int element = 0; element < a_Nelem - 1; ++element) {
    int start = 0;
    int end = lobbatto_order;
    if (element == 0)
      start = 1; // omit first DVR point
    int offset = (lobbatto_order - 1) * element;
    for (int i = start; i < end; ++i) {
      int ith_index = (i - 1) + offset;
      for (int j = start; j < end; ++j) {
        int jth_index = (j - 1) + offset;
        std::complex<double> tmp = 0.0;
        for (int m = 0; m < lobbatto_order; ++m) {
          tmp += lobbatto_derivatives[m * end + i] *
                 lobbatto_derivatives[m * end + j] *
                 a_lobattoQuad->getWeight(m);
        }
        m_laplacian[ith_index * m_Nbas + jth_index] +=
            tmp * 2.0 /
            ((m_complexbounds[element + 1] - m_complexbounds[element]) *
             sqrt(m_weights[ith_index] * m_weights[jth_index]));
      }
    }
  }
  // finding the derivatives on the complex contour so using rotated radau
  // points and weights
  auto radau_derivatives =
      std::make_unique<std::complex<double>[]>(radau_order * lobbatto_order);
  for (int i = 0; i < radau_order; ++i) {
    radau_derivatives[i * radau_order + i] = 0.0;
    for (int k = 0; k < radau_order; ++k) {
      if (i != k)
        radau_derivatives[i * radau_order + i] +=
            1.0 / (a_radauQuad->getPoint(i) - a_radauQuad->getPoint(k));
    }
    for (int j = 0; j < radau_order; ++j) {
      if (i != j) {
        std::complex<double> tmp =
            1.0 / (a_radauQuad->getPoint(j) - a_radauQuad->getPoint(i));
        for (int k = 0; k < radau_order; ++k) {
          if ((k != i) && (k != j))
            tmp *= (a_radauQuad->getPoint(i) - a_radauQuad->getPoint(k)) /
                   (a_radauQuad->getPoint(j) - a_radauQuad->getPoint(k));
        }
        radau_derivatives[i * radau_order + j] = tmp;
      }
    }
  }

  std::complex<double> eitm2 = 1.0 / pow(eit, 2);
  int offset = (lobbatto_order - 1) * (m_Nelem - 1);
  for (int i = 0; i < radau_order; ++i) {
    int ith_index = (i - 1) + offset;
    for (int j = 0; j < radau_order; ++j) {
      int jth_index = (j - 1) + offset;
      std::complex<double> tmp = 0.0;
      for (int m = 0; m < radau_order; ++m) {
        tmp += radau_derivatives[m * radau_order + i] *
               radau_derivatives[m * radau_order + j] * radtempW[m];
      }
      if (i == j) {
        tmp = tmp -
              alphaRad * conj(eit) *
                  (radtempW[i] * radau_derivatives[i * radau_order + j] +
                   radtempW[j] * radau_derivatives[j * radau_order + i]) +
              (pow(alphaRad, 2)) * eitm2 * radtempW[i];
        m_laplacian[ith_index * m_Nbas + jth_index] +=
            tmp / sqrt(m_weights[ith_index] * m_weights[jth_index]);
      } else {
        tmp = tmp - alphaRad * conj(eit) *
                        (radtempW[i] * radau_derivatives[i * radau_order + j] +
                         radtempW[j] * radau_derivatives[j * radau_order + i]);
        m_laplacian[ith_index * m_Nbas + jth_index] +=
            tmp / sqrt(m_weights[ith_index] * m_weights[jth_index]);
      }
    }
  }

} // end of constructor

/// Destructor.
FEMDVR::~FEMDVR() {}

/// Getter for points.
std::complex<double> FEMDVR::getPoint(int index) const {
  return m_points[index];
}

/// Getter for weights.
std::complex<double> FEMDVR::getWeight(int index) const {
  return m_weights[index];
}

/// Getter for weights.
std::complex<double> FEMDVR::getLaplacian(int index) const {
  return m_laplacian[index];
}

/// Get number of basis functions
size_t &FEMDVR::getNbas() { return m_Nbas; }

void FEMDVR::print() const {
  std::cout << "\n";
  std::cout << "*** FEM-DVR Grid***" << std::endl;
  std::cout << "\n";
  std::cout << "Number of Finite-Elements = " << m_Nelem << "; "
            << "Number of Basis Functions = " << m_Nbas << "\n";
  std::cout << "\n";
  for (int i = 0; i < m_Nbas; ++i) {
    std::cout << "Point " << m_points[i] << std::setw(20) << "Weights "
              << m_weights[i] << "\n";
  }
}
