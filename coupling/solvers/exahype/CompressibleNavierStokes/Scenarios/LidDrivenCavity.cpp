#include "LidDrivenCavity.h"
void NavierStokes::LidDrivenCavity::initialValues(const double* const x, const PDE& ns,
                                           Variables& vars) {
    vars.rho() = 1.0;
#if DIMENSIONS == 2
    vars.j(0.0, 0.0);
#else
    vars.j(0.0, 0.0, 0.0);
#endif
    // Mach = u/c
    const auto mach = 0.1;
    const auto wallSpeed = 1.0;
    const auto speedOfSound = wallSpeed/mach;
    const auto pressure = (speedOfSound * speedOfSound)/gamma;
    vars.E() = ns.evaluateEnergy(vars.rho(), pressure, vars.j());
    ns.setZ(vars.data(), 0.0); // Disable advection.
}

NavierStokes::BoundaryType NavierStokes::LidDrivenCavity::getBoundaryType(
    int faceId) {
    if (faceId == 3) {
      return BoundaryType::movingWall;
    }
    return BoundaryType::wall;
}
