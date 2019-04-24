#ifndef COMPRESSIBLENAVIERSTOKES_SETUPHELPER_H
#define COMPRESSIBLENAVIERSTOKES_SETUPHELPER_H

#include "AMR/AMRSettings.h"
#include "AbstractNavierStokesSolver_ADERDG.h"
#include "Scenarios/Scenario.h"
#include "Scenarios/ScenarioFactory.h"

#include <memory>
#include <tuple>

namespace NavierStokes {
struct ScenarioConfig {
  PDE ns;
  std::string scenarioName;
  std::unique_ptr<Scenario> scenario;
  AMRSettings amrSettings;
};

ScenarioConfig parseConfig(const std::vector<std::string> &cmdlineargs,
                           const exahype::parser::ParserView &constants,
                           int numberOfVariables, int numberOfParameters,
                           int numberOfGlobalObservables);

}  // namespace NavierStokes

#endif  // COMPRESSIBLENAVIERSTOKES_SETUPHELPER_H
