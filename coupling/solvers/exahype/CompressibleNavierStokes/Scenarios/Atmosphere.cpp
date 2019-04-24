#include "Atmosphere.h"

double NavierStokes::computeHydrostaticPressure(const PDE& ns, double g,
                                                double posZ,
                                                double backgroundT) {
  return std::pow(
      (1 - 1 / ns.gamma) *
          (std::pow(std::pow(ns.gamma, ns.gamma / (ns.gamma - 1)) *
                        ns.referencePressure,
                    (ns.gamma - 1) / ns.gamma) /
               (ns.gamma - 1) -
           g * std::pow(ns.referencePressure, ns.gasConstant / ns.c_p) * posZ /
               (ns.gasConstant * backgroundT)),
      ns.gamma / (ns.gamma - 1));
}

double NavierStokes::potentialTToT(const PDE& ns, double pressure,
                                   double potentialT) {
  const double potToT =
      std::pow((pressure / ns.referencePressure), ns.gasConstant / ns.c_p);
  return potentialT * potToT;
}

double NavierStokes::computeHydrostaticPressureGradient(
    const PDE& ns, double g, double posZ, double backgroundT) {
  const auto p = computeHydrostaticPressure(ns, g, posZ, backgroundT);
  return -g * std::pow(ns.referencePressure, ns.gasConstant / ns.c_p) *
         std::pow(p, 1.0 / ns.gamma) / (ns.gasConstant * backgroundT);
}

double NavierStokes::computeHydrostaticTemperatureGradient(const PDE& ns, double g,
                                             double posZ, double backgroundT) {
  const auto p = computeHydrostaticPressure(ns, g, posZ, backgroundT);
  const auto pGrad =
      computeHydrostaticPressureGradient(ns, g, posZ, backgroundT);
  return (ns.gasConstant * backgroundT *
          std::pow(p / ns.referencePressure, ns.gasConstant / ns.c_p) * pGrad) /
         (ns.c_p * p);
}
