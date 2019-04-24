#include <array>

#include "Atmosphere.h"
#include "DensityCurrent.h"

void NavierStokes::DensityCurrent::initialValues(const double* const x,
                                                 const PDE& ns,
                                                 Variables& vars) {
  // For details see:
  // TODO(Lukas) Add reference to correct paper!
  const auto posX = x[0];
  const auto posZ = (DIMENSIONS == 2) ? x[1] : x[2];

  const double backgroundPotentialT = getBackgroundPotentialTemperature();  // [K], background potential temperature

  double potentialT = backgroundPotentialT;

  const auto centerX = 0;
  const auto centerZ = 3000;
  const auto sizeX = 4000;
  const auto sizeZ = 2000;

  const auto distX = (posX - centerX) / sizeX;
  const auto distZ = (posZ - centerZ) / sizeZ;
  const auto r = std::sqrt(distX * distX + distZ * distZ);
  const auto r_c = 1.0;

  if (r <= r_c) {
    const auto pertubationSize = -15;
    const auto pi = std::acos(-1);
    potentialT += pertubationSize / 2 * (1 + std::cos(pi * r));
  }

  // Air is initially at rest.
#if DIMENSIONS == 2
  vars.j(0, 0);
#elif DIMENSIONS == 3
  vars.j(0, 0, 0);
#endif

  const auto pressure = computeHydrostaticPressure(ns, getGravity(), posZ, backgroundPotentialT);
  const auto temperature = potentialTToT(ns, pressure, potentialT);
  vars.rho() = pressure / (ns.gasConstant * temperature);
  vars.E() = ns.evaluateEnergy(vars.rho(), pressure, vars.j(), 0.0, ns.getHeight(vars.data()));

  if (ns.useBackgroundState) {
    // TODO(Lukas) Refactor?
    // Then compute background state (without pot.T. pertubation
    const auto backgroundTemperature = potentialTToT(ns, pressure, potentialT);
    const auto backgroundRho = pressure / (ns.gasConstant * backgroundTemperature);
    ns.setBackgroundState(vars.data(), backgroundRho, pressure);
  }
}

void NavierStokes::DensityCurrent::source(
    const tarch::la::Vector<DIMENSIONS, double>& x, double t, const PDE& ns,
    const double* const Q, double* S) {
  Scenario::source(x, t, ns, Q, S);
  double rhoPertubation = Q[0];
  if (ns.useBackgroundState) {
    auto backgroundRho = 0.0;
    auto backgroundPressure = 0.0;
    std::tie(backgroundRho, backgroundPressure) = ns.getBackgroundState(Q);

    rhoPertubation -= backgroundRho;
  }
  S[DIMENSIONS] = -1 * rhoPertubation * getGravity();

  // Only use this source term if the gravitational force is not already
  // included in the pressure.
  if (!ns.useGravity) {
    S[DIMENSIONS + 1] = -1 * Q[2] * getGravity();
  }
}

double NavierStokes::DensityCurrent::getGamma() const { return gamma; }

double NavierStokes::DensityCurrent::getPr() const { return Pr; }

double NavierStokes::DensityCurrent::getC_v() const { return c_v; }

double NavierStokes::DensityCurrent::getC_p() const { return c_p; }

double NavierStokes::DensityCurrent::getGasConstant() const {
  return gasConstant;
}

double NavierStokes::DensityCurrent::getReferencePressure() const {
  return referencePressure;
}

double NavierStokes::DensityCurrent::getGravity() const {
  return 9.81;
}

double NavierStokes::DensityCurrent::getBackgroundPotentialTemperature() const {
    return 300;
}

NavierStokes::BoundaryType NavierStokes::DensityCurrent::getBoundaryType(int faceId) {
  // For boundaries in Z direction, we need to reconstruct
  // the temperature flux coming from the boundary.
#if DIMENSIONS == 2
  const bool isZFace = faceId == 2 || faceId == 3;
#else
  const bool isZFace = faceId == 4 || faceId == 5;
#endif
  if (isZFace) {
    return BoundaryType::hydrostaticWall;
  }
  return BoundaryType::hydrostaticWall;
  //return BoundaryType::wall;

}
