#include "Couette.h"

void NavierStokes::Couette::initialValues(const double* const x, const PDE& ns,
                                           Variables& vars) {
    vars.rho() = 1.0;
    vars.j(0.0, 0.0, 0.0);

    // Mach = u/c
    const auto mach = 0.1;
    const auto wallSpeed = 0.5;
    const auto speedOfSound = wallSpeed/mach;
    const auto pressure = (speedOfSound * speedOfSound)/ns.gamma;
    vars.E() = ns.evaluateEnergy(vars.rho(), pressure, vars.j());
    ns.setZ(vars.data(), 0.0); // Disable advection.
}

void NavierStokes::Couette::analyticalSolution(const double* const x, double t, const PDE& ns,
                          Variables& vars, double* gradState) {  
    vars.rho() = 0;
    vars.j(0) = getVelocity(x[2]);
    vars.j(1) = 0.0;
    vars.j(2) = 0.0;

    // Mach = u/c
    const auto mach = 0.1;
    const auto speedOfSound = wallVelocity/mach;
    const auto pressure = (speedOfSound * speedOfSound)/ns.gamma;
    vars.E() = ns.evaluateEnergy(vars.rho(), pressure, vars.j());

    const auto gradSize = vars.variables() * 3;
    std::fill_n(gradState, gradSize, 0.0);

    gradState[idxGradQ(2, 1)] = getVelocity(x[2] + .5) - getVelocity(x[2] - .5);
}

double NavierStokes::Couette::getVelocity(double z){
    // Todo get from config
    const double kinVisc = 1.73982;
    const double channelheight = 50.0;
    const double wallVelocity = 0.5;

    double vel(0);
    vel = wallVelocity*(1.0 - z/channelheight);
    const double pi = 3.141592653589793238;
    double sum = 0.0;
    for (int k = 1; k < 30; k++){
      sum += 1.0/k * sin(k*pi*z/channelheight) * exp(-k*k * pi*pi/(channelheight*channelheight) * kinVisc * t);
    }
    vel = vel - 2.0*wallVelocity/pi * sum;
    return vel;
}

NavierStokes::BoundaryType NavierStokes::Couette::getBoundaryType(
    int faceId) {
    if (faceId == 4) {  // zmin
      return BoundaryType::movingWall;
    }
    else if (faceId == 5) {  // zmax
      return BoundaryType::wall;
    }
    return BoundaryType::analytical;
}
