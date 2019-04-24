#ifndef COMPRESSIBLENAVIERSTOKES_CRITERION_H
#define COMPRESSIBLENAVIERSTOKES_CRITERION_H

#include <string>
#include <vector>

#include "kernels/limiter/generic/Limiter.h"

#include "AMRSettings.h"
#include "PDE.h"
#include "totalVariation.h"
#include "VarianceHelper.h"

namespace NavierStokes {

/** NEW **/
template <typename GlobalObservables>
void resetGlobalObservables(GlobalObservables& obs)  {
  obs.gobs(0) = -1.0;
  obs.gobs(1) = -1.0;
  obs.gobs(2) = 0.0;
}

template <
  typename GlobalObservables,
  typename ReadOnlyGlobalObservables
>
void mergeGlobalObservables(
    GlobalObservables&         obs,
    ReadOnlyGlobalObservables& other)  {
    if (obs.size() == 0) {
    return;
  }

  assertion2(obs.size() == other.size(),
             obs.size(), other.size());

  auto *reducedGlobalObservables = obs.data();
  auto const *curGlobalObservables = other.data();
  
  const auto mean0 = reducedGlobalObservables[0];
  const auto mean1 = curGlobalObservables[0];
  const auto var0 = reducedGlobalObservables[1];
  const auto var1 = curGlobalObservables[1];
  const auto count0 = reducedGlobalObservables[2];
  const auto count1 = curGlobalObservables[2];

  auto mergedMean = 0.0;
  auto mergedVariance = 0.0;

  std::tie(mergedMean, mergedVariance) =
      mergeVariance(mean0, mean1, var0, var1, count0, count1);
  reducedGlobalObservables[0] = mergedMean;
  reducedGlobalObservables[1] = mergedVariance;
  reducedGlobalObservables[2] = count0 + count1;
    
  }

template <typename GlobalObservables>
void mapGlobalObservablesDG(
    GlobalObservables& globalObservables,
    const double* const Q,
    const tarch::la::Vector<DIMENSIONS,double>& dx,
    const std::string &scenarioName,
    const PDE &ns,
    const AMRSettings &amrSettings,
    int Order,
    int NumberOfVariables,
    int NumberOfParameters)  {
  
  if (globalObservables.size() == 0) {
    return;
  }

  auto indicator = amrSettings.indicator;
  auto useTV = amrSettings.useTotalVariation;

  indicator = IndicatorVariable::potentialTemperature;
  useTV = true;
  assert(indicator == IndicatorVariable::potentialTemperature);
  assert(useTV);

  auto computeIndicator = [&](const double *const Q) {
    const auto vars = ReadOnlyVariables{Q};
    if (indicator == IndicatorVariable::potentialTemperature) {
      const auto pressure =
          ns.evaluatePressure(vars.E(), vars.rho(), vars.j(), ns.getZ(Q),
                  ns.getHeight(Q));
      const auto temperature = ns.evaluateTemperature(vars.rho(), pressure);
      return ns.evaluatePotentialTemperature(temperature, pressure);
    }
    if (indicator == IndicatorVariable::Z) {
      return ns.getZ(Q)/vars.rho();
    }
    auto backgroundRho = 0.0;
    auto backgroundPressure = 0.0;
    std::tie(backgroundRho, backgroundPressure) = ns.getBackgroundState(Q);


    if (indicator == IndicatorVariable::rho) {
      return vars.rho() - backgroundRho;
    } else {
      // Pressure
      return ns.evaluatePressure(vars.E(), vars.rho(), vars.j(),
              ns.getZ(Q), ns.getHeight(Q)) -
             backgroundPressure;
    }
  };

  auto *observables = globalObservables.data();
  // TODO(Lukas) Implement global observables for 3D!

  assertion(useTV);

  if (useTV) {
    const auto tv =
        totalVariation(Q, Order, NumberOfVariables, NumberOfParameters, dx,
                       amrSettings.correctForVolume, computeIndicator);
    observables[0] = tv;
    observables[1] = 0;
    observables[2] = 1;
  } else {
    throw -1;
  }
}

template <typename GlobalObservables>
void mapGlobalObservablesFV(
    GlobalObservables&                          globalObservables,
    const double* const                         luh,
    const tarch::la::Vector<DIMENSIONS,double>& cellSize)  {
  // TODO(Lukas): Please implement
}

template <typename GlobalObservables, int PatchSize, int GhostLayerWidth,
	  int NumberOfVariables, int NumberOfParameters>
void mapGlobalObservablesFV(
					   GlobalObservables& globalObservables,
					   const double* const Q,
					   const tarch::la::Vector<DIMENSIONS, double>& dx,
					   const std::string& scenarioName,
					   const PDE& ns,
					   const AMRSettings& amrSettings) {
  // Idea: Map FV-Data to DG-Representation, and compute global variables there.
  // This way, we can use the same implementation of total variance for both
  // solvers. Additionally, we have the same TV before and after limiting. It is
  // really slow though, so it might be a good idea to approximate TV
  // differently.
  // TODO(Lukas) Skip mapping if not using TV, e.g. non-gradient based crit.
  if (globalObservables.size() == 0) {
    return;
  }

  assert((PatchSize - 1) % 2 == 0);
  constexpr int Order = (PatchSize - 1) / 2;
  constexpr int NumberOfData = NumberOfVariables + NumberOfParameters;

  // TODO(Lukas) Stack allocate this?
  auto QDG =
      std::vector<double>(std::pow(Order + 1, DIMENSIONS) * NumberOfData);
  kernels::limiter::generic::c::projectOnDGSpace<Order + 1,
       NumberOfData, GhostLayerWidth>(Q, QDG.data());

  NavierStokes::mapGlobalObservablesDG(globalObservables,
				       QDG.data(),
				       dx,
				       scenarioName,
				       ns,
				       amrSettings,
				       Order,
				       NumberOfVariables,
				       NumberOfParameters);
}

}  // namespace NavierStokes

#endif  // COMPRESSIBLENAVIERSTOKES_CRITERION_H
