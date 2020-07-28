#include <benchmark/benchmark.h>
#include <iostream>
#include <vector>

static void lobattoBench(benchmark::State &s) {
  std::vector<double> m_points, m_weights;
std::vector<std::vector<double> > m_derivatives;
  int m_kind = 1, m_kpts = 2, m_dvr = 5;
  double m_alpha = 0.0, m_beta= 0.0;
  double *x = new double[m_dvr];
  double *w = new double[m_dvr];
  double *src = new double[m_dvr];
  double *endpts = new double[2];
  endpts[0] = -1.0;
  endpts[1] = 1.0;
  gaussq_(&m_kind, &m_dvr, &m_alpha, &m_beta, &m_kpts, endpts, src, x, w);
  for (int i = 0; i < m_dvr; ++i) {
    m_points.push_back(x[i]);
    m_weights.push_back(w[i]);
  }
  m_derivatives.resize(m_dvr);
  for (int i = 0; i < m_dvr; ++i) {
    m_derivatives[i].resize(m_dvr);
    m_derivatives[i][i] = 0.0;
    for (int k = 0; k < m_dvr; ++k) {
      if (i != k)
        m_derivatives[i][i] = m_derivatives[i][i] + 1.0 / (x[i] - x[k]);
    }
    for (int j = 0; j < m_dvr; ++j) {
      if (i != j) {
        double temp = 1.0 / (x[j] - x[i]);
        for (int k = 0; k < m_dvr; ++k) {
          if ((k != i) && (k != j))
            temp = temp * (x[i] - x[k]) / (x[j] - x[k]);
        }
        m_derivatives[i][j] = temp;
      }
    }
  }
}

// Register the benchmark
// DenseRange allows us to generate a set of inputs
// ReportAggregatesOnly allows us to limit our output
BENCHMARK(lobattoBench)->DenseRange(13, 26)->ReportAggregatesOnly(true);

// This is basically our main function
BENCHMARK_MAIN();
