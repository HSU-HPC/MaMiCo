#include "ScenarioFactory.h"

// Convergence
#include "Scenarios/ConvergenceTest/ConvergenceTest.h"

// CFD: NS/Euler
// TODO(Lukas) Organize
#include "Scenarios/DoubleShockTube.h"
#include "Scenarios/EntropyWave.h"
#include "Scenarios/SmoothWave.h"
#include "Scenarios/SodShockTube.h"
#include "Scenarios/Stokes.h"
#include "Scenarios/TaylorGreen.h"
#include "Scenarios/LidDrivenCavity.h"
#include "Scenarios/ABCFlow.h"

// Atmospheric Flows
#include "Scenarios/TwoBubbles.h"
#include "Scenarios/DensityCurrent.h"
#include "Scenarios/CosineBubble.h"

// Coupled Scenarios
#include "Scenarios/CouplingTest.h"
#include "Scenarios/Detonation.h"

#include <stdexcept>

NavierStokes::ScenarioFactory::ScenarioPtr
NavierStokes::ScenarioFactory::createScenario(const std::string& scenarioName) {
  if (scenarioName == "sod-shock-tube") {
    return std::move(ScenarioPtr(new SodShockTube()));
  } else if (scenarioName == "double-shock-tube") {
    return std::move(ScenarioPtr(new DoubleShockTube()));
  } else if (scenarioName == "smooth-wave") {
    return std::move(ScenarioPtr(new SmoothWave()));
  } else if (scenarioName == "entropy-wave") {
    return std::move(ScenarioPtr(new EntropyWave()));
  } else if (scenarioName == "stokes") {
    return std::move(ScenarioPtr(new Stokes()));
  } else if (scenarioName == "taylor-green") {
    return std::move(ScenarioPtr(new TaylorGreen()));
  } else if (scenarioName == "lid-driven-cavity") {
    return std::move(ScenarioPtr(new LidDrivenCavity()));
  } else if (scenarioName == "abc-flow") {
    return std::move(ScenarioPtr(new ABCFlow()));
  } else if (scenarioName == "two-bubbles") {
    return std::move(ScenarioPtr(new TwoBubbles()));
  } else if (scenarioName == "density-current") {
    return std::move(ScenarioPtr(new DensityCurrent()));
  } else if (scenarioName == "cosine-bubble") {
    return std::move(ScenarioPtr(new CosineBubble()));
  } else if (scenarioName == "convergence") {
    return std::move(ScenarioPtr(new ConvergenceTest()));
  } else if (scenarioName == "coupling-test") {
    return std::move(ScenarioPtr(new CouplingTest()));
  } else if (scenarioName == "detonation") {
    return std::move(ScenarioPtr(new Detonation()));
  }

  throw std::invalid_argument(scenarioName + " is not a valid scenario!");
}
