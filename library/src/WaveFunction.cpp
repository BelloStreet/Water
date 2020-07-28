// cSpell:disable
/* #include "../include/WaveFunction.hpp" */

/* WaveFunction::WaveFunction(){}; */

/* WaveFunction::~WaveFunction(){}; */

/// Gives the DVR representation of a function at the point x on the grid.
/* std::complex<double> WaveFunction::dvrRep( */
/*     const double &a_x, const std::vector<std::complex<double>> &a_func) { */
/*   assert((a_x > m_complexbounds[0].real()) && */
/*          (a_x <= m_complexbounds[m_Nelem].real())); */

/*   int FEM_loc = 0; */
/*   for (int element = 0; element < m_Nelem; ++element) { */
/*     if (a_x < element) { */
/*       ++FEM_loc; */
/*     } else */
/*       break; */
/*   } */
/*   int index, start, end; */
/*   std::complex<double> fac; */
/*   std::vector<std::complex<double>> tempPoints, tempWeights; */
/*   if (FEM_loc == 0) { */
/*     index = 0; */
/*     start = 1; */
/*     end = m_lobattoQuad->getDVROrder(); */
/*     fac = 1.0; */
/*     tempPoints.push_back(0.0); */
/*     tempWeights.push_back(0.0); */
/*     std::copy(&m_points[0] + 1, &m_points[0] + m_lobattoQuad->getDVROrder() -
 * 1, */
/*               std::back_inserter(tempPoints)); */
/*     std::copy(&m_weights[0] + 1, */
/*               &m_weights[0] + m_lobattoQuad->getDVROrder() - 1, */
/*               std::back_inserter(tempWeights)); */
/*   } else if (FEM_loc == m_Nelem - 1) { */
/*     index = m_Nbas - m_radauQuad->getDVROrder(); */
/*     start = 0; */
/*     end = m_radauQuad->getDVROrder(); */
/*     fac = std::exp(-m_alphaRad * conj(m_eit) * (a_x - m_R0)); */
/*     std::copy(&m_points[0] + m_Nbas - m_radauQuad->getDVROrder(), */
/*               &m_points[0] + m_Nbas, std::back_inserter(tempPoints)); */
/*     std::copy(&m_weights[0] + m_Nbas - m_radauQuad->getDVROrder(), */
/*               &m_weights[0] + m_Nbas, std::back_inserter(tempWeights)); */
/*   } else { */
/*     index = FEM_loc * (m_lobattoQuad->getDVROrder() - 1); */
/*     start = 0; */
/*     end = m_lobattoQuad->getDVROrder(); */
/*     fac = 1.0; */
/*     std::copy(&m_points[0] + (FEM_loc) * (m_lobattoQuad->getDVROrder() - 2),
 */
/*               &m_points[0] + (FEM_loc + 1) * (m_lobattoQuad->getDVROrder() -
 * 1), */
/*               std::back_inserter(tempPoints)); */
/*     std::copy( */
/*         &m_weights[0] + (FEM_loc) * (m_lobattoQuad->getDVROrder() - 2), */
/*         &m_weights[0] + (FEM_loc + 1) * (m_lobattoQuad->getDVROrder() - 1),
 */
/*         std::back_inserter(tempWeights)); */
/*   } */
/*   std::complex<double> accum = 0.0; */
/*   for (int j = start; j < end; ++j) { */
/*     std::complex<double> prod = 1.0; */
/*     for (int k = start; k < end; ++k) { */
/*       if (k != j) { */
/*         prod *= (a_x - tempPoints[k]) / (tempPoints[j] - tempPoints[k]); */
/*       } */
/*     } */
/*     accum += prod * a_func[index] / sqrt(tempWeights[index]); */
/*     ++index; */
/*   } */
/*   return accum * fac; */

/* }  // dvrRep */
