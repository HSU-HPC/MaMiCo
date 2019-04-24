#include "totalVariation.h"

double totalVariation(const double* Q, int order, int numberOfVariables,
                      int numberOfParameters,
                      const tarch::la::Vector<DIMENSIONS, double>& dx,
                      bool correctForVolume) {
  const auto basisSize = order + 1;
  const auto numberOfData = numberOfVariables + numberOfParameters;
#if DIMENSIONS == 2
  const auto idx_lQi = kernels::idx3(basisSize, basisSize, numberOfData);
#else
  const auto idx_lQi =
      kernels::idx4(basisSize, basisSize, basisSize, numberOfData);
#endif

  const auto& quadratureWeights = kernels::gaussLegendreWeights[order];
  double tv = 0.0;
#if DIMENSIONS == 2
  // x direction (independent from the y derivatives)
  for (int k = 0; k < basisSize; k++) {  // k == y
    // Matrix operation
    for (int l = 0; l < basisSize; l++) {  // l == x
      const auto w = quadratureWeights[k] * quadratureWeights[l];
      for (int m = 0; m < numberOfVariables; m++) {
        double curGrad = 0.0;
        for (int n = 0; n < basisSize; n++) {  // n == matmul x
          const auto t = Q[idx_lQi(k, n, m)] * kernels::dudx[order][l][n];
          curGrad += t;
        }
        tv += w * std::abs(curGrad);
      }
    }
  }

  // y direction (independent from the x derivatives)
  for (int k = 0; k < basisSize; k++) {
    // Matrix operation
    for (int l = 0; l < basisSize; l++) {  // l == y
      const auto w = quadratureWeights[k] * quadratureWeights[l];
      for (int m = 0; m < numberOfVariables; m++) {
        double curGrad = 0.0;
        for (int n = 0; n < basisSize; n++) {  // n = matmul y
          const auto t = Q[idx_lQi(n, k, m)] *
                         kernels::dudx[order][l][n]; /* l,n: transpose */
          curGrad += t;
        }
        tv += w * std::abs(curGrad);
      }
    }
  }
#else
  // x direction (independent from the y and z derivatives)
  for (int j = 0; j < basisSize; j++) {    // z
    for (int k = 0; k < basisSize; k++) {  // y
      // Matrix operation
      for (int l = 0; l < basisSize; l++) {  // x
        const auto w =
            quadratureWeights[j] * quadratureWeights[k] * quadratureWeights[l];
        for (int m = 0; m < numberOfVariables; m++) {
          double curGrad = 0.0;
          for (int n = 0; n < basisSize; n++) {
            const auto t = Q[idx_lQi(j, k, n, m)] *
                           kernels::dudx[order][l][n];
            curGrad += t;
          }
          tv += w * std::abs(curGrad);
        }
      }
    }
  }

  // y direction (independent from the x and z derivatives)
  for (int j = 0; j < basisSize; j++) {    // z
    for (int k = 0; k < basisSize; k++) {  // x
      // Matrix operation
      for (int l = 0; l < basisSize; l++) {  // y
        const auto w =
            quadratureWeights[j] * quadratureWeights[k] * quadratureWeights[l];
        for (int m = 0; m < numberOfVariables; m++) {
          double curGrad = 0.0;
          for (int n = 0; n < basisSize; n++) {
            const auto t = Q[idx_lQi(j, n, k, m)] *
                           kernels::dudx[order][l][n];
            curGrad += t;
          }
          tv += w * std::abs(curGrad);
        }
      }
    }
  }

  // z direction (independent from the x and y derivatives)
  for (int j = 0; j < basisSize; j++) {    // y
    for (int k = 0; k < basisSize; k++) {  // x
      // Matrix operation
      for (int l = 0; l < basisSize; l++) {  // z
        const auto w =
            quadratureWeights[j] * quadratureWeights[k] * quadratureWeights[l];
        for (int m = 0; m < numberOfVariables; m++) {
          double curGrad = 0.0;
          for (int n = 0; n < basisSize; n++) {
            const auto t = Q[idx_lQi(n, j, k, m)] *
                           kernels::dudx[order][l][n];
            curGrad += t;
          }
          tv += w * std::abs(curGrad);
        }
      }
    }
  }
#endif

  if (correctForVolume) {
#if DIMENSIONS == 2
    const auto volume = dx[0] * dx[1];
#else
    const auto volume = dx[0] * dx[1] * dx[2];
#endif
    return volume * tv;
  }

  return tv;
}
