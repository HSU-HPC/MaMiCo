#include "Stokes.h"

void NavierStokes::Stokes::initialValues(const double* const x, const PDE& ns,
                                         Variables& vars) {
  vars.rho() = 1.0;
#if DIMENSIONS == 2
  vars.j(1.0, 0.0);
#elif DIMENSIONS == 3
  vars.j(1.0, 0.0, 0.0);
#endif
  const double pressure = 100 / ns.gamma;
  vars.E() = ns.evaluateEnergy(vars.rho(), pressure, vars.j());
}
void NavierStokes::Stokes::analyticalSolution(const double* const x,
                                              const double t, const PDE& ns,
                                              Variables& vars,
                                              double* gradState) {
  kernels::idx2 idxGradQ(DIMENSIONS, vars.SizeVariables);
  const double rho = 1;
  const double v = ns.referenceViscosity / rho;
  vars.rho() = 1.0;
  vars.j(0) = std::erf(x[1] / std::sqrt(2 * v * t));
  vars.j(1) = 0.0;
#if DIMENSIONS == 3
  vars.j(2) = 0.0;
#endif
  const double pressure = 100 / ns.gamma;
  vars.E() = ns.evaluateEnergy(vars.rho(), pressure, vars.j());

  // Assuming rho is constant.
  const double pi = std::acos(-1); // acos not constexpr for intel compiler

  gradState[idxGradQ(1, 1)] =
      std::sqrt(2) *
      exp(-1.0 / 2.0 * pow(x[1], 2) / (ns.referenceViscosity * t)) /
      (std::sqrt(pi) * sqrt(ns.referenceViscosity) * sqrt(t));
  const int EIdx = DIMENSIONS + 1;
  gradState[idxGradQ(1, EIdx)] =
      1.0 * std::sqrt(2) *
      std::exp(-1.0 / 2.0 * std::pow(x[1], 2) / (ns.referenceViscosity * t)) *
      std::erf((1.0 / 2.0) * std::sqrt(2) * x[1] /
               (std::sqrt(ns.referenceViscosity) * std::sqrt(t))) /
      (std::sqrt(pi) * std::sqrt(ns.referenceViscosity) * sqrt(t));
}

NavierStokes::BoundaryType NavierStokes::Stokes::getBoundaryType(int faceId) {
  if (faceId == 2) {
    return BoundaryType::wall;
  }
  return BoundaryType::analytical;
}
