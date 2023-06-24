// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MAMICO_COUPLING_SCENARIO_SCENARIO_H_
#define _MAMICO_COUPLING_SCENARIO_SCENARIO_H_

#include <iostream>
#include <string>
#include "coupling/CouplingMDDefinitions.h"
#include "coupling/services/ParallelTimeIntegrationService.h"

class Scenario {
public:
  Scenario(std::string scenarioname) : _scenarioname(scenarioname) { std::cout << "Run " << scenarioname << "..." << std::endl; }
  virtual ~Scenario() { std::cout << "Shut down " << _scenarioname << std::endl; }

  virtual void run() = 0;
  virtual void init() = 0;
  virtual void runOneCouplingCycle(int cycle) = 0;

protected:
  std::unique_ptr<coupling::services::ParallelTimeIntegrationService<MY_LINKEDCELL, 3>> _timeIntegrationService;

private:
  const std::string _scenarioname;
};

#endif // _MAMICO_COUPLING_SCENARIO_SCENARIO_H_
