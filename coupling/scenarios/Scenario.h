// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_SCENARIOS_SCENARIO_H_
#define _MOLECULARDYNAMICS_COUPLING_SCENARIOS_SCENARIO_H_

#include <iostream>
#include <string>
#include "tarch/logging/Logger.h"

class Scenario {
public:
  Scenario(std::string scenarioname) : _scenarioname(scenarioname), _logger(tarch::logging::Logger(scenarioname)) {
    _logger.info("Running {}", _scenarioname);
  }
  
  virtual ~Scenario() {
    _logger.info("Shutting down {}", _scenarioname);
  }

  virtual void run() = 0;

private:
  const std::string _scenarioname;
  tarch::logging::Logger _logger;
};

#endif // _MOLECULARDYNAMICS_COUPLING_SCENARIOS_SCENARIO_H_
