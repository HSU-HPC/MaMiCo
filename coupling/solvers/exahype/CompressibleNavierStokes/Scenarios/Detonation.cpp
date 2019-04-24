#include "Detonation.h"

void NavierStokes::Detonation::initialValues(const double* const x,
                                             const PDE& ns, Variables& vars) {
  assert(ns.useAdvection);

  const auto centerX = 0.0;
  const auto centerY = 0.0;

  // B: Burned, U: Unburned
  const auto rhoB = 1.4;
  const auto rhoU = 0.887565;
  const auto pB = 1.0;
  const auto pU = 0.191709;
  const auto ZB = 0.0;
  const auto ZU = 1.0;

  const auto alpha = std::atan2(x[1] - centerX, x[0] - centerY); // Angle in polar coordinates
  const auto uB = 0.0;
  const auto uU = -0.577350 * std::cos(alpha);
  const auto vB = 0.0;
  const auto vU = -0.577350 * std::sin(alpha);

  const auto distX = x[0];
  const auto distY = x[1];
  const auto distToCenter = std::sqrt(
          distX * distX + distY * distY
          );
  const auto radiusBurned = 0.3;
  const auto isBurned = distToCenter < radiusBurned;

  double pressure = -1;
  if (isBurned) {
     vars.rho() = rhoB;
     pressure = pB;
#if DIMENSIONS == 2
     vars.j(uB * rhoB, vB * rhoB);
#else
    vars.j(uB * rhoB, vB * rhoB, 0.0);
#endif
     ns.setZ(vars.data(), ZB * rhoB);
  } else {
     vars.rho() = rhoU;
     pressure = pU;
#if DIMENSIONS == 2
     vars.j(uU * rhoU, vU * rhoU);
#else
     vars.j(uU * rhoU, vU * rhoU, 0.0);
#endif
     ns.setZ(vars.data(), ZU * rhoU);
  }

  vars.E() = ns.evaluateEnergy(vars.rho(), pressure, vars.j(), ns.getZ(vars.data()));

  assert(ns.getZ(vars.data()) == 0 || ns.getZ(vars.data()) == rhoU);
}
double NavierStokes::Detonation::getMolecularDiffusionCoeff() const {
  return 0.0;
}

double NavierStokes::Detonation::getQ0() const {
  return 1.0;
}
void NavierStokes::Detonation::source(
    const tarch::la::Vector<DIMENSIONS, double>& x, double t, const PDE& ns,
    const double* const Q, double* S) {
  Scenario::source(x, t, ns, Q, S);
  const auto vars = ReadOnlyVariables{Q};

  const auto pressure = ns.evaluatePressure(vars.E(), vars.rho(), vars.j(), ns.getZ(Q));
  const auto T = ns.evaluateTemperature(vars.rho(), pressure);

  const auto reactionTimeScale = 0.1; // TODO(Lukas) Try stiff ts?
  const auto ignitionT = 0.26;
  //const auto ignitionT = 0.4;
  if (T > ignitionT) {
    // Critical temperature, start decay.
    ns.setZ(S, -1.0/reactionTimeScale * ns.getZ(Q));
  } else {
    ns.setZ(S, 0.0);
  }
}

double NavierStokes::Detonation::getGamma() const { return 1.4; }

double NavierStokes::Detonation::getPr() const {
  return 0.75;
}

double NavierStokes::Detonation::getC_v() const { return 2.5; }

double NavierStokes::Detonation::getC_p() const { return 3.5; }

double NavierStokes::Detonation::getGasConstant() const {
  return 1.0;
}

bool NavierStokes::Detonation::getUseAdvection() const {
  return true;
}

NavierStokes::BoundaryType NavierStokes::Detonation::getBoundaryType(int faceId) {
  return BoundaryType::wall;
}