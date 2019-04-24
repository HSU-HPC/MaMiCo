#ifndef COMPRESSIBLENAVIERSTOKES_TOTALVARIATION_H
#define COMPRESSIBLENAVIERSTOKES_TOTALVARIATION_H
#include <vector>

#include "peano/utils/Dimensions.h"
#include "tarch/la/Vector.h"

#include "kernels/GaussLegendreQuadrature.h"
#include "kernels/KernelUtils.h"
#include "kernels/aderdg/generic/Kernels.h"
using computeFunc = double (*)(const double* const);

double totalVariation(const double* Q, int order, int numberOfVariables,
                      int numberOfParameters,
                      const tarch::la::Vector<DIMENSIONS, double>& dx,
                      bool correctForVolume);
/*
double totalVariation(const double* Q, int order, int numberOfVariables, int
numberOfParameters, const tarch::la::Vector<DIMENSIONS,double>& dx, bool
correctForVolume, computeFunc mapObservable);
*/

template <typename F>
double totalVariation(const double* Q, int order, int numberOfVariables,
                      int numberOfParameters,
                      const tarch::la::Vector<DIMENSIONS, double>& dx,
                      bool correctForVolume, F mapObservable) {
  const auto basisSize = order + 1;
  const auto numberOfData = numberOfVariables + numberOfParameters;

  auto observable = std::vector<double>(std::pow(basisSize, DIMENSIONS));
#if DIMENSIONS == 2
  const auto idxQ = kernels::idx3(basisSize, basisSize, numberOfData);
  const auto idxObservable = kernels::idx3(basisSize, basisSize, 1);
  for (int k = 0; k < basisSize; k++) {
    for (int l = 0; l < basisSize; l++) {
      observable[idxObservable(k, l, 0)] = mapObservable(Q + idxQ(k, l, 0));
    }
  }
#else
  const auto idxQ =
      kernels::idx4(basisSize, basisSize, basisSize, numberOfData);
  const auto idxObservable = kernels::idx4(basisSize, basisSize, basisSize, 1);
  for (int j = 0; j < basisSize; j++) {
    for (int k = 0; k < basisSize; k++) {
      for (int l = 0; l < basisSize; l++) {
        observable[idxObservable(j, k, l, 0)] =
            mapObservable(Q + idxQ(j, k, l, 0));
      }
    }
  }
#endif

  return totalVariation(observable.data(), order, 1, 0, dx, correctForVolume);
}

inline double forwardDiff(double l, double c, double h) {
  assert(std::isfinite( (1./h) * (l - c)));
  return (1./h) * (l - c);
}

inline double centralDiff(double l, double r, double h) {
  assert(std::isfinite( (1./(2*h)) * (r - l)));
  return (1./(2 * h)) * (r - l);
}

inline double backwardDiff(double c, double r, double h) {
  assert(std::isfinite( (1./h) * (c - r)));
  return (1./h) * (c - r);
}

inline double stableDiff(double l,
		  double c,
		  double r,
		  int idxC,  
		  double h,
		  size_t ghostLayerWidth,
		  size_t patchSize) {
  // TODO(Lukas) Check for off by one errors
  // TODO(Lukas) Use central differences for center?
  // TODO(Lukas) Check if order of arguments is correct.
  
  // idxC must be signed to avoid underflow!
  const auto idxL = idxC - 1;
  const auto idxR = idxC + 1;
  
  //std::cout << idxC << std::endl;
  if (idxL <= static_cast<int>(ghostLayerWidth)) {
    return backwardDiff(c, r, h);
  } else if (idxR >= static_cast<int>(patchSize + ghostLayerWidth)) {
    return forwardDiff(l, c, h);
  } else {
    return centralDiff(l, r, h);
  }
}

#endif  // COMPRESSIBLENAVIERSTOKES_TOTALVARIATION_H
