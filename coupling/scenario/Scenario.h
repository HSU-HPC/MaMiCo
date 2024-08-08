// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MAMICO_COUPLING_SCENARIO_SCENARIO_H_
#define _MAMICO_COUPLING_SCENARIO_SCENARIO_H_

class Scenario;

#include "coupling/CouplingMDDefinitions.h"
#include "coupling/solvers/CouetteSolver.h"
#include <iostream>
#include <string>

class Scenario {
public:
  Scenario(std::string scenarioname) : _scenarioname(scenarioname) { std::cout << "Run " << scenarioname << "..." << std::endl; }
  virtual ~Scenario() { std::cout << "Shut down " << _scenarioname << std::endl; }

  virtual void run() = 0;
  virtual void init() = 0;
  virtual void runOneCouplingCycle(int cycle) = 0;

  virtual coupling::solvers::AbstractCouetteSolver<3>* getSolver() = 0;
  const coupling::services::ParallelTimeIntegrationService<3>* getTimeIntegrationService() const { return _timeIntegrationService.get(); }

protected:
  std::unique_ptr<coupling::services::ParallelTimeIntegrationService<3>> _timeIntegrationService;

private:
  const std::string _scenarioname;
};

#endif // _MAMICO_COUPLING_SCENARIO_SCENARIO_H_
