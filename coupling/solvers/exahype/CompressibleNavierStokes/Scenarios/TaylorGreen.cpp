#include "TaylorGreen.h"
#include "tarch/Assertions.h"

void NavierStokes::TaylorGreen::analyticalSolution(const double* const x,
                                                   const double t,
                                                   const PDE& ns,
                                                   Variables& vars,
                                                   double* gradState) {
  assertionEqualsMsg(DIMENSIONS, 2, "Taylor-Green is only supported for 2D!");

  kernels::idx2 idxGradQ(DIMENSIONS, vars.SizeVariables);

  const double Ft = std::exp(-2 * ns.referenceViscosity * t);

  vars.rho() = 1.0;
  vars.j(0) = 1 * std::sin(x[0]) * std::cos(x[1]) * Ft;
  vars.j(1) = -1 * std::cos(x[0]) * std::sin(x[1]) * Ft;

  const double p_0 = 100. / ns.gamma;
  const double pressure = p_0 + ((vars.rho() * Ft * Ft) / 4) *
                                    (std::cos(2 * x[0]) + std::cos(2 * x[1]));

  vars.E() = ns.evaluateEnergy(vars.rho(), pressure, vars.j());

  // Assuming rho is constant.
  // j(0)
  gradState[idxGradQ(0, 1)] = 1 * std::cos(x[0]) * std::cos(x[1]) * Ft;
  gradState[idxGradQ(1, 1)] = -1 * std::sin(x[0]) * std::sin(x[1]) * Ft;

  // j(1)
  gradState[idxGradQ(0, 2)] = 1 * std::sin(x[0]) * std::sin(x[1]) * Ft;
  gradState[idxGradQ(1, 2)] = -1 * std::cos(x[0]) * std::cos(x[1]) * Ft;

  // E
  gradState[idxGradQ(0, 3)] =
      std::pow(Ft, 2) *
      ((1.0 / 4.0) * (1.0 * ns.gamma - 1.0) *
           (std::sin(2 * x[0] - 2 * x[1]) + std::sin(2 * x[0] + 2 * x[1])) +
       0.5 * std::sin(2 * x[0])) /
      (ns.gamma - 1);

  gradState[idxGradQ(1, 3)] =
      std::pow(Ft, 2) *
      ((1.0 / 4.0) * (1.0 * ns.gamma - 1.0) *
           (-std::sin(2 * x[0] - 2 * x[1]) + std::sin(2 * x[0] + 2 * x[1])) +
       0.5 * std::sin(2 * x[1])) /
      (ns.gamma - 1);
}

NavierStokes::BoundaryType NavierStokes::TaylorGreen::getBoundaryType(
    int faceId) {
  return BoundaryType::analytical;
}