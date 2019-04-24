#ifndef COMPRESSIBLENAVIERSTOKES_AMRSETTINGS_H
#define COMPRESSIBLENAVIERSTOKES_AMRSETTINGS_H

#include <string>

namespace NavierStokes {

enum class IndicatorVariable { rho, pressure, Z, potentialTemperature };

IndicatorVariable parseIndicatorVariable(const std::string& name);

struct AMRSettings {
  bool useAMR = false;
  bool useTotalVariation = true;
  bool correctForVolume = true;
  IndicatorVariable indicator = IndicatorVariable::pressure;
  double factorRefine = 1.5;
  double factorErase = 0.5;
};

}  // namespace NavierStokes

#endif  // COMPRESSIBLENAVIERSTOKES_AMRSETTINGS_H
