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

#define __MAMICO_STRINGIFY_EXPAND(x) #x
#define MAMICO_STRINGIFY(x) __MAMICO_STRINGIFY_EXPAND(x)

class Scenario {
public:
  Scenario(std::string scenarioname) : _scenarioname(scenarioname) {
    getRootRank();
    if (_isRootRank) {
      std::cout << "Run " << scenarioname << "..." << std::endl;
      std::cout << "MaMiCo git commit hash = " << MAMICO_STRINGIFY(MAMICO_COMMIT_HASH) << std::endl;
    }
  }
  virtual ~Scenario() {
    if (_isRootRank) {
      std::cout << "Shut down " << _scenarioname << std::endl;
    }
  }

  virtual void run() = 0;
  virtual void init() = 0;
  virtual void runOneCouplingCycle(int cycle) = 0;
  virtual void equilibrateMicro() = 0;

  virtual coupling::solvers::AbstractCouetteSolver<3>* getSolver() = 0;
  const coupling::services::ParallelTimeIntegrationService<3>* getTimeIntegrationService() const { return _timeIntegrationService.get(); }

protected:
  std::unique_ptr<coupling::services::ParallelTimeIntegrationService<3>> _timeIntegrationService;

  /** @brief initialises all MPI variables  */
  void getRootRank() {
    int rank = 0;
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
    _isRootRank = (rank == 0);
  }

  /** @brief if this is the world global root process */
  bool _isRootRank;

private:
  const std::string _scenarioname;
};

#endif // _MAMICO_COUPLING_SCENARIO_SCENARIO_H_
