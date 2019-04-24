#ifndef __DIFFUSIVE_RIEMANN_SOLVER_HEADER__
#define __DIFFUSIVE_RIEMANN_SOLVER_HEADER__

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

// included in ../../Kernels.h

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>

#include "kernels/GaussLegendreQuadrature.h"
#include "kernels/KernelUtils.h"

/**
 * We implement a very simple Rusanov scheme with scalar dissipation
 * (smax*Id).
 *
 * We need to consider material parameters
 * in QL and QR.
 * We don't need to consider material parameters
 * in FL,FR.
 */
template <bool useNCP, typename SolverType>
void riemannSolverNonlinear(
    SolverType& solver, double* FL, double* FR, const double* const QL,
    const double* const QR,
    const tarch::la::Vector<DIMENSIONS, double>& characteristicLength,
    const double dt, const int direction) {
  constexpr int numberOfVariables = SolverType::NumberOfVariables;
  constexpr int numberOfParameters = SolverType::NumberOfParameters;
  constexpr int numberOfData = numberOfVariables + numberOfParameters;
  constexpr int order = SolverType::Order;
  constexpr int basisSize = order + 1;

  using namespace kernels;
  idx2 idx_FLR(basisSize, numberOfVariables);
  idx2 idx_QLR(basisSize, numberOfData);

  // Compute maximal eigenvalues for entire element.
  double maxHyperbolicEigenvalue = std::numeric_limits<double>::min();
  double maxDiffusiveEigenvalue = std::numeric_limits<double>::min();
  for (int j = 0; j < basisSize; j++) {
    // do not need to store material parameters
    double LL[numberOfVariables] = {0.0};
    double LR[numberOfVariables] = {0.0};
    solver.eigenvalues(&QL[idx_QLR(j, 0)], direction, LL);
    solver.eigenvalues(&QR[idx_QLR(j, 0)], direction, LR);

    // skip parameters
    // hyperbolic eigenvalues
    std::transform(LL, LL + numberOfVariables, LL, std::abs<double>);
    std::transform(LR, LR + numberOfVariables, LR, std::abs<double>);
    const double smax_L = *std::max_element(LL, LL + numberOfVariables);
    const double smax_R = *std::max_element(LR, LR + numberOfVariables);
    maxHyperbolicEigenvalue =
        std::max(maxHyperbolicEigenvalue, std::max(smax_L, smax_R));

    // diffusive eigenvalues
    solver.viscousEigenvalues(&QL[idx_QLR(j, 0)], direction, LL);
    solver.viscousEigenvalues(&QR[idx_QLR(j, 0)], direction, LR);
    std::transform(LL, LL + numberOfVariables, LL, std::abs<double>);
    std::transform(LR, LR + numberOfVariables, LR, std::abs<double>);
    const double smaxDiffusive_L =
        *std::max_element(LL, LL + numberOfVariables);
    const double smaxDiffusive_R =
        *std::max_element(LR, LR + numberOfVariables);
    maxDiffusiveEigenvalue = std::max(
        maxDiffusiveEigenvalue, std::max(smaxDiffusive_L, smaxDiffusive_R));
  }

  // Compute penalty term for entire element.
  const auto pi = std::acos(-1); // not constexpr for intel compilers
  const double factor = 2 * (2 * order + 1) / (characteristicLength[direction]
          * std::sqrt(0.5 * pi));
  const double penalty =
      maxHyperbolicEigenvalue + factor * maxDiffusiveEigenvalue;

  // compute fluxes (and fluctuations for non-conservative PDEs)
  double Qavg[numberOfData];
  idx2 idx_gradQ(DIMENSIONS, numberOfVariables);
  double gradQ[DIMENSIONS][numberOfVariables] = {0.0};
  double ncp[numberOfVariables] = {0.0};
  {
    for (int j = 0; j < basisSize; j++) {
      if (useNCP) {  // we don't use matrixB but the NCP call here.
        for (int l = 0; l < numberOfVariables; l++) {
          gradQ[direction][l] = QR[idx_QLR(j, l)] - QL[idx_QLR(j, l)];
          Qavg[l] = 0.5 * (QR[idx_QLR(j, l)] + QL[idx_QLR(j, l)]);
        }

        solver.nonConservativeProduct(Qavg, gradQ[0], ncp);
      }

      // skip parameters
      for (int k = 0; k < numberOfVariables; k++) {
        FL[idx_FLR(j, k)] =
            0.5 * (FR[idx_FLR(j, k)] + FL[idx_FLR(j, k)]) -
            0.5 * penalty * (QR[idx_QLR(j, k)] - QL[idx_QLR(j, k)]);
        assertion6(std::isfinite(FL[idx_FLR(j, k)]), j, k, factor, penalty,
                   maxHyperbolicEigenvalue, maxDiffusiveEigenvalue);
        if (useNCP) {
          FR[idx_FLR(j, k)] = FL[idx_FLR(j, k)] - 0.5 * ncp[k];
          FL[idx_FLR(j, k)] = FL[idx_FLR(j, k)] + 0.5 * ncp[k];
        } else {
          FR[idx_FLR(j, k)] = FL[idx_FLR(j, k)];
        }
      }

    }
  }
}

#endif
