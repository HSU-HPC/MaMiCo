#include "SmoothWave.h"
void NavierStokes::SmoothWave::initialValues(const double* const x,
                                             const PDE& ns, Variables& vars) {
#if DIMENSIONS == 2
  vars.j(0.0, 0.0);
#elif DIMENSIONS == 3
  vars.j(0.0, 0.0, 0.0);
#endif
  const auto distX = x[0] - 0.5;
  const auto distY = x[1] - 0.5;
  const auto dist = distX * distX + distY * distY;
  vars.rho() = 1 - dist;
  const double pressure = 1 - dist;

  vars.E() = ns.evaluateEnergy(vars.rho(), pressure, vars.j());
}

NavierStokes::BoundaryType NavierStokes::SmoothWave::getBoundaryType(int faceId) {
  return BoundaryType::freeSlipWall;
}

