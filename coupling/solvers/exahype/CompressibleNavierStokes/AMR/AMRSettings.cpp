#include "AMRSettings.h"

NavierStokes::IndicatorVariable NavierStokes::parseIndicatorVariable(
    const std::string& name) {
  if (name == "rho") {
    return IndicatorVariable::rho;
  } else if (name == "pressure") {
    return IndicatorVariable::pressure;
  } else if (name == "Z") {
    return IndicatorVariable::Z;
  } else if (name == "potential-temperature") {
    return IndicatorVariable::potentialTemperature;
  } else {
    throw -1;
  }
}
