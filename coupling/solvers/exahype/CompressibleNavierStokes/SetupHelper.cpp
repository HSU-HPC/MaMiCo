#include "SetupHelper.h"
#include "exahype/parser/ParserView.h"

NavierStokes::ScenarioConfig NavierStokes::parseConfig(
    const std::vector<std::string>& cmdlineargs,
    const exahype::parser::ParserView& constants,
    int numberOfVariables,
    int numberOfParameters,
    int numberOfGlobalObservables) {
  // TODO(Lukas) Error handling?
  assert(constants.isValueValidString("scenario"));

  double referenceViscosity;
  assert(constants.isValueValidDouble("viscosity"));
  referenceViscosity = constants.getValueAsDouble("viscosity");

  auto scenarioName = constants.getValueAsString("scenario");
  auto scenario = ScenarioFactory::createScenario(scenarioName);

  const bool useGravity =
      constants.getValueAsBoolOrDefault("use-gravity", false);
  const bool useBackgroundState =
      constants.getValueAsBoolOrDefault("use-background-state", false);
  assert(useGravity ||
         !useBackgroundState);  // Background state only works with gravity.

  // AMR-Settings
  auto amrSettings = AMRSettings{};
  amrSettings.useAMR = constants.getValueAsBoolOrDefault("use-amr", false);
  if (amrSettings.useAMR) {
    amrSettings.useTotalVariation = constants.getValueAsBool("use-tv-amr");
    amrSettings.correctForVolume = constants.getValueAsBool("use-tv-volume-correction");
    amrSettings.indicator =
        parseIndicatorVariable(constants.getValueAsString("amr-indicator"));
    amrSettings.factorRefine = constants.getValueAsDouble("amr-factor-refine");
    amrSettings.factorErase = constants.getValueAsDouble("amr-factor-erase");
  }

  // Make sure we have the correct number of variables/parameters:
  auto numberOfNecessaryVariables = 1 + DIMENSIONS + 1;
  if (scenario->getUseAdvection()) {
    ++numberOfNecessaryVariables;
  }
  if (numberOfVariables != numberOfNecessaryVariables) {
    throw - 1;
  }

  auto numberOfNecessaryParameters = 0;
  if (useGravity) {
    ++numberOfNecessaryParameters;
  }
  if (useBackgroundState) {
    numberOfNecessaryParameters += 2;
  }
  if (numberOfParameters != numberOfNecessaryParameters) {
    throw - 1;
  }

  auto numberOfNecessaryGlobalObservables = amrSettings.useAMR ? 3 : 0;
  if (numberOfGlobalObservables < numberOfNecessaryGlobalObservables) {
    throw - 1;
  }

  auto ns = PDE(referenceViscosity, *scenario, useGravity, useBackgroundState);

  return {ns, scenarioName, std::move(scenario), amrSettings};
}