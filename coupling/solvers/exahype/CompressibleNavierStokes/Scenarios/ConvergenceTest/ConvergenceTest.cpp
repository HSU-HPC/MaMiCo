#include "ConvergenceTest.h"
#include "ManufacturedSolution.h"

void NavierStokes::ConvergenceTest::analyticalSolution(const double* const x,
                                                       double t, const PDE& ns,
                                                       Variables& vars,
                                                       double* gradState) {
  evaluateQ(ns.gamma, ns.referenceViscosity, t, x[0], x[1], vars.data());
  evaluateGradQ(ns.gamma, ns.referenceViscosity, t, x[0], x[1], gradState);
}

void NavierStokes::ConvergenceTest::source(
    const tarch::la::Vector<DIMENSIONS, double>& x, double t, const PDE& ns,
    const double* const Q, double* S) {
  constexpr auto NumberOfVariables = DIMENSIONS + 2;  // TODO(Lukas) generalise?
  std::fill_n(S, NumberOfVariables, 0.0);
  evaluateSource(ns.gasConstant, ns.gamma,
                 ns.evaluateHeatConductionCoeff(ns.referenceViscosity),
                 ns.referenceViscosity, t, x[0], x[1], S);
}

NavierStokes::BoundaryType NavierStokes::ConvergenceTest::getBoundaryType(
    int faceId) {
  return BoundaryType::analytical;
}
