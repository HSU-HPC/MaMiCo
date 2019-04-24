#ifndef COMPRESSIBLENAVIERSTOKES_ATMOSPHERE_H
#define COMPRESSIBLENAVIERSTOKES_ATMOSPHERE_H

#include "../PDE.h"

namespace NavierStokes {
    // TODO(Lukas) Rename temperature to potential temperature
double computeHydrostaticPressure(const PDE &ns, double g, double posZ,
                                  double backgroundTemperature);
double computeHydrostaticPressureGradient(const PDE &ns, double g, double posZ,
                                          double backgroundTemperature);
double computeHydrostaticTemperatureGradient(const PDE &ns, double g,
                                             double posZ,
                                             double backgroundTemperature);
double potentialTToT(const PDE &ns, double pressure, double potentialT);
}  // namespace NavierStokes

#endif  // COMPRESSIBLENAVIERSTOKES_ATMOSPHERE_H
