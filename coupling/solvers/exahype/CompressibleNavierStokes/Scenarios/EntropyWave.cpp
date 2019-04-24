#include "EntropyWave.h"

void NavierStokes::EntropyWave::analyticalSolution(const double* const x,
                                                   double t,
                                                   const NavierStokes::PDE& ns,
                                                   Variables& vars,
                                                   double* gradState) {
  // Solution from chapter 7.13.2 in "I do like CFD, VOL.1" by Katate Masatsuka.
  assert(ns.referenceViscosity == 0.0);
  // Set gradient to zero as this is only needed
  // to validate simulations with zero viscosity!
  const auto gradSize = vars.variables() * DIMENSIONS;
  std::fill_n(gradState, gradSize, 0.0);

  const auto width = 0.3;
  const auto rhoOffset = 0.5;
  const auto pressure = 1.0;
#if DIMENSIONS == 2
  const auto v = tarch::la::Vector<DIMENSIONS, double>(0.5, 0.5);
  const auto vSum = v[0] + v[1];
  const auto positionSum = x[0] + x[1];
#else
  const auto v = tarch::la::Vector<DIMENSIONS, double>(0.5, 0.5, 0.5);
  const auto vSum = v[0] + v[1] + v[2];
  const auto positionSum = x[0] + x[1] + x[2];
#endif

  const auto pi = std::acos(-1);

  vars.rho() = rhoOffset + width * std::sin(pi * (positionSum - vSum * t));
  vars.j(vars.rho() * v);
  vars.E() = ns.evaluateEnergy(vars.rho(), pressure, vars.j());
}

NavierStokes::BoundaryType NavierStokes::EntropyWave::getBoundaryType(
    int faceId) {
  return BoundaryType::analytical;
}
