/**
 * This file is part of the ExaHyPE project.
 * Copyright (c) 2016  http://exahype.eu
 * All rights reserved.
 *
 * The project has received funding from the European Union's Horizon 
 * 2020 research and innovation programme under grant agreement
 * No 671698. For copyrights and licensing, please consult the webpage.
 *
 * Released under the BSD 3 Open Source License.
 * For the full license text, see LICENSE.txt
 **/
#include <array>
#include <algorithm>
#include <cmath>
#include <limits>

#include "tarch/la/Vector.h" 
#include "kernels/KernelUtils.h"

// Debugging
#include <iostream>

template <typename SolverType>
double stableDiffusiveTimeStepSize(
    SolverType& solver,
    const double* const luh,
    const tarch::la::Vector<DIMENSIONS, double>& dx) {
  constexpr int numberOfVariables  = SolverType::NumberOfVariables;
  constexpr int numberOfParameters = SolverType::NumberOfParameters;
  constexpr int numberOfData       = numberOfVariables+numberOfParameters;
  constexpr int order              = SolverType::Order;
  constexpr int basisSize          = order+1;
  constexpr double cflFactor       = SolverType::CFL;

  // Obtained from von Neumann analysis.
  // PNPM[N] <= 1/(2N+1)
  constexpr double PNPM[10]  = {1.0,   0.33,  0.17, 0.1,  0.069,
                                0.045, 0.038, 0.03, 0.02, 0.015};

  auto dt = std::numeric_limits<double>::max();
#if DIMENSIONS == 2
	kernels::idx3 idx_luh(basisSize, basisSize, numberOfData);
#elif DIMENSIONS == 3
	kernels::idx4 idx_luh(basisSize, basisSize, basisSize, numberOfData);
#endif

	auto minDx = std::numeric_limits<double>::max();
	for (int i = 0; i < DIMENSIONS; i++) {
		minDx = std::min(minDx, dx[i]);
	}

	auto maxHyperbolicEigenvalue = std::numeric_limits<double>::min();
	auto maxDiffusiveEigenvalue = std::numeric_limits<double>::min();

  // Iterate over dofs.
  for (int i = 0; i < basisSize; i++) {
    for (int j = 0; j < basisSize; j++) {
#if DIMENSIONS == 3
      for (int k = 0; k < basisSize; k++) {
#endif
	for (int dim = 0; dim < DIMENSIONS; dim++) {
	  // First compute max eigenvalues of hyperbolic and diffusive part.
	  auto hyperbolicEigenvalues = std::array<double,numberOfVariables>();
	  auto diffusiveEigenvalues = std::array<double,numberOfVariables>();
#if DIMENSIONS == 2
        const auto idx = idx_luh(i,j,0);
#elif DIMENSIONS == 3
        const auto idx = idx_luh(i,j,k,0);
#endif
	  solver.eigenvalues(luh + idx, dim, hyperbolicEigenvalues.data());
	  solver.viscousEigenvalues(luh + idx, dim, diffusiveEigenvalues.data());

	  for (const auto eigen : hyperbolicEigenvalues) {
	    maxHyperbolicEigenvalue = std::max(maxHyperbolicEigenvalue, std::abs(eigen));
	  }
	  for (const auto eigen : diffusiveEigenvalues) {
	    maxDiffusiveEigenvalue = std::max(maxDiffusiveEigenvalue, std::abs(eigen));
	  }

	} // dim
#if DIMENSIONS == 3
	  } // k
#endif
    } // j
   } // i
  //dt = (cflFactor * minDx * PNPM[order]) / DIMENSIONS * 1./ (
  //			maxHyperbolicEigenvalue + maxDiffusiveEigenvalue * 2 * (2 * order + 1) / minDx);
  dt = (cflFactor * minDx * PNPM[order]) / DIMENSIONS * 1./ (
			maxHyperbolicEigenvalue + maxDiffusiveEigenvalue * 2 * (1/PNPM[order]) / minDx);

  return dt;
}
