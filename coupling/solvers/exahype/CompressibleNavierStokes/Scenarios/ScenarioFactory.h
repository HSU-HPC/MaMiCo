#ifndef COMPRESSIBLENAVIERSTOKES_SCENARIOFACTORY_H
#define COMPRESSIBLENAVIERSTOKES_SCENARIOFACTORY_H

#include "Scenario.h"
#include <memory>
#include <string>

namespace NavierStokes {
class ScenarioFactory {
 private:
  using ScenarioPtr = std::unique_ptr<Scenario>;

 public:
  static ScenarioPtr createScenario(const std::string& scenarioName);
};
}  // namespace NavierStokes

#endif  // COMPRESSIBLENAVIERSTOKES_SCENARIOFACTORY_H
