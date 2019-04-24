#include "NavierStokesSolver_FV.h"

#include "NavierStokesSolver_FV_Variables.h"
#include "NavierStokesSolver_ADERDG_Variables.h"

#include "PDE.h"
#include "AMR/Criterion.h"

#include "Scenarios/Atmosphere.h"
#include "SetupHelper.h"

tarch::logging::Log NavierStokes::NavierStokesSolver_FV::_log( "NavierStokes::NavierStokesSolver_FV" );

void NavierStokes::NavierStokesSolver_FV::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
  auto parsedConfig = parseConfig(cmdlineargs, constants, NumberOfVariables, NumberOfParameters, NumberOfGlobalObservables);
  ns = std::move(parsedConfig.ns);
  scenarioName = std::move(parsedConfig.scenarioName);
  scenario = std::move(parsedConfig.scenario);
  amrSettings = std::move(parsedConfig.amrSettings);
}

void NavierStokes::NavierStokesSolver_FV::adjustSolution(const double* const x,const double t,const double dt, double* Q) {
  if (tarch::la::equals(t, 0.0)) {
    ns.setHeight(Q, x[DIMENSIONS-1]);
    ns.setBackgroundState(Q, 0.0, 0.0);

    AbstractNavierStokesSolver_ADERDG::Variables vars(Q);
    scenario->initialValues(x, ns, vars);
    for (int i = 0; i < vars.variables(); ++i) {
      assertion2(std::isfinite(Q[i]), i, Q[i]);
    }
  }

  const auto vars = ReadOnlyVariables{Q};
  const auto pressure = ns.evaluatePressure(vars.E(), vars.rho(), vars.j(), ns.getZ(Q));
  assertion5(pressure >= 0, pressure, vars.E(), vars.rho(), vars.j(), ns.getZ(Q));
}

void NavierStokes::NavierStokesSolver_FV::eigenvalues(const double* const Q, const int dIndex, double* lambda) {
  ns.evaluateEigenvalues(Q, dIndex, lambda);
}

void NavierStokes::NavierStokesSolver_FV::boundaryValues(
    const double* const x,
    const double t,const double dt,
    const int faceIndex,
    const int normalNonZero,
    const double* const stateIn,
    double* stateOut) {
  // No slip, 2D
  ReadOnlyVariables varsIn(stateIn);
  Variables varsOut(stateOut);

  // Not supported!
  assertion(scenario->getBoundaryType(faceIndex) != BoundaryType::analytical);

  // All bcs here are walls!
  assertion(scenario->getBoundaryType(faceIndex) == BoundaryType::wall ||
                 scenario->getBoundaryType(faceIndex) == BoundaryType::hydrostaticWall ||
                 scenario->getBoundaryType(faceIndex) == BoundaryType::movingWall ||
                 scenario->getBoundaryType(faceIndex) == BoundaryType::freeSlipWall);

  ns.setHeight(stateOut, x[DIMENSIONS-1]);
  ns.setBackgroundState(stateOut, 0.0, 0.0);
  // Need to reconstruct the background state in this case.
  // Extrapolating does not work!
  assert(!ns.useBackgroundState || scenario->getBoundaryType(faceIndex) ==
      BoundaryType::hydrostaticWall);

  // Rho/E extrapolated, velocity mirrored.
  std::copy_n(stateIn, NumberOfVariables, stateOut);

  if (scenario->getBoundaryType(faceIndex) == BoundaryType::hydrostaticWall ||
      scenario->getBoundaryType(faceIndex) == BoundaryType::freeSlipWall) {
    // Normal velocity zero after Riemann.
    varsOut.j(normalNonZero) = -varsIn.j(normalNonZero);
  } else {
    // No-slip
    // All velocities zero after Riemann.
    varsOut.j(0) = -varsIn.j(0);
    varsOut.j(1) = -varsIn.j(1);
#if DIMENSIONS == 3
    varsOut.j(2) = -varsIn.j(2);
#endif
  }

  if (scenario->getBoundaryType(faceIndex) == BoundaryType::movingWall) {
    const auto wallSpeed = 1.0;
    varsOut.j(0) = 2 * wallSpeed - varsIn.j(0);
  }

  if (scenario->getBoundaryType(faceIndex) == BoundaryType::hydrostaticWall) {
    // Note: This boundary condition is incorrect for the viscous case, as we do not
    // reconstruct the temperature diffusion coming from the wall!
    const auto posZ = x[DIMENSIONS-1];

    // We also need to reconstruct the temperature at the border.
    // This corresponds to a heated wall.
    const auto pressure = computeHydrostaticPressure(ns, scenario->getGravity(),
                                                     posZ, scenario->getBackgroundPotentialTemperature());
    const auto T = potentialTToT(ns, pressure, scenario->getBackgroundPotentialTemperature());
    const auto rho = pressure / (ns.gasConstant * T);
    // TODO(Lukas) What should we do in case of an advection-scenarios?

    // TODO(Lukas) Is background state necessary here?
    auto E = -1;
    if (ns.useGravity) {
      E = ns.evaluateEnergy(rho, pressure, varsOut.j(), ns.getZ(stateIn), x[DIMENSIONS - 1]);
    } else {
      E = ns.evaluateEnergy(rho, pressure, varsOut.j(), ns.getZ(stateIn));
    }
    varsOut.E() = varsIn.E();//E;
    varsOut.rho() = varsIn.rho();// rho;
    ns.setBackgroundState(stateOut, varsOut.rho(), ns.evaluatePressure(
            varsOut.E(), varsOut.rho(), varsOut.j(), 0.0, ns.getHeight(stateOut)
            ));
  }
}

void NavierStokes::NavierStokesSolver_FV::viscousFlux(const double* const Q,const double* const gradQ, double** F) {
  ns.evaluateFlux(Q, gradQ, F, true);
}

void NavierStokes::NavierStokesSolver_FV::viscousEigenvalues(const double* const Q, const int dIndex, double* lambda) {
  ns.evaluateDiffusiveEigenvalues(Q, dIndex, lambda);
}

//You can either implement this method or modify fusedSource
void NavierStokes::NavierStokesSolver_FV::algebraicSource(const tarch::la::Vector<DIMENSIONS, double>& x, double t, const double *const Q, double *S) {
  // TODO: Actually use coordinates!
  scenario->source(x, t, ns, Q, S);
}

void NavierStokes::NavierStokesSolver_FV::resetGlobalObservables(GlobalObservables& globalObservables) const  {
  NavierStokes::resetGlobalObservables(globalObservables);
}

void NavierStokes::NavierStokesSolver_FV::mapGlobalObservables(
    GlobalObservables&                          globalObservables,
    const double* const                         luh,
    const tarch::la::Vector<DIMENSIONS,double>& cellSize)  const {
  
  // Ignore efficiency for now.
  // TODO: Are derivatives correct?
  // TODO: Is ghost layer handled correctly?
  // TODO: Don't hardcore indicator variable.
  assert(DIMENSIONS == 2);
  constexpr auto numberOfData = NumberOfVariables + NumberOfParameters;
  kernels::idx3 idx(PatchSize+2*GhostLayerWidth,PatchSize+2*GhostLayerWidth,numberOfData);
  kernels::idx2 idx_obs(PatchSize+2*GhostLayerWidth,PatchSize+2*GhostLayerWidth);

  kernels::idx3 idx_slope(PatchSize+2*GhostLayerWidth,
                          PatchSize+2*GhostLayerWidth,
                          DIMENSIONS);
  const auto subcellSize = (1./PatchSize) * cellSize; 
  assert(std::isfinite(subcellSize[0]));
  assert(std::isfinite(subcellSize[1]));

  // Compute slope by finite differences
  constexpr auto variablesPerPatch = (PatchSize+2*GhostLayerWidth)*(PatchSize+2*GhostLayerWidth)*numberOfData;
  constexpr int patchBegin = GhostLayerWidth; // patchBegin cell is inside domain
  constexpr int patchEnd = patchBegin+PatchSize; // patchEnd cell is outside domain
  double slope[variablesPerPatch*DIMENSIONS] = {0.0};
  auto observables = std::vector<double>((PatchSize+2*GhostLayerWidth)*(PatchSize+2*GhostLayerWidth));

  const size_t numberOfIndicators = 1;
  auto computeIndicator = [&](const double *const Q) {
    const auto vars = ReadOnlyVariables{Q};
    const auto pressure =
    ns.evaluatePressure(vars.E(), vars.rho(), vars.j(), ns.getZ(Q),
                        ns.getHeight(Q));
    const auto temperature = ns.evaluateTemperature(vars.rho(), pressure);
    return ns.evaluatePotentialTemperature(temperature, pressure);
  };

  // Compute indicator variables
  for (int j = patchBegin; j < patchEnd; j++) {
    for (int k = patchBegin; k < patchEnd; k++) {
      observables[idx_obs(j,k)] = computeIndicator(luh + idx(j,k,0));
    }
  }

  
  // Compute slopes    
  // slopex
  for (int j = patchBegin; j < patchEnd; j++) { // y
    for (int k = patchBegin; k < patchEnd; k++) { // x
      const auto left = observables[idx_obs(j, k-1)];
      const auto center = observables[idx_obs(j, k+0)];
      const auto right = observables[idx_obs(j, k+1)];

      slope[idx_slope(j, k, 0)] =
        stableDiff(left, center, right,
                   j,
                   subcellSize[0],
                   GhostLayerWidth,
                   PatchSize);
    }
  }
  // slopey
  for (int j = patchBegin; j < patchEnd; j++) { // y
    for (int k = patchBegin; k < patchEnd; k++) { // x
      const auto left = observables[idx_obs(j-1, k)];
      const auto center = observables[idx_obs(j, k)];
      const auto right = observables[idx_obs(j+1, k)];

      slope[idx_slope(j, k, 0)] =
        stableDiff(left, center, right,
                   k,
                   subcellSize[1],
                   GhostLayerWidth,
                   PatchSize);

    }
  }

  double meanReduced = -1.0;
  double varReduced = -1.0;
  double n = 0.0;

  for (int i = patchBegin - 1; i < patchEnd; ++i) {
    for (int j = patchBegin-1; j < patchEnd; ++j) {
      for (int d = 0; d < DIMENSIONS; ++d) {
        const auto curTv = slope[idx_slope(i,j,d)];
        std::tie(meanReduced, varReduced) = mergeVariance(meanReduced,
							  curTv,
							  varReduced,
							  0.0,
							  n,
							  1);
        n += 1;
      }
    }
  }



  auto *rawGlobalObservables = globalObservables.data();
  rawGlobalObservables[0] = meanReduced;
  rawGlobalObservables[1] = varReduced;
  rawGlobalObservables[2] = n;
}

void NavierStokes::NavierStokesSolver_FV::mergeGlobalObservables(
    GlobalObservables&         globalObservables,
    ReadOnlyGlobalObservables& otherObservables) const  {
  NavierStokes::mergeGlobalObservables(globalObservables,otherObservables);
}
