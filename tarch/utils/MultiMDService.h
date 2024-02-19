// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _TARCH_UTILS_MULTIMDSERVICE_H_
#define _TARCH_UTILS_MULTIMDSERVICE_H_

#include "tarch/TarchDefinitions.h"
#include "tarch/la/Vector.h"
#include <cstdlib>
#include <numeric>
#include <vector>
#if (TARCH_PARALLEL == TARCH_YES)
#include <mpi.h>
#endif

namespace tarch {
namespace utils {
template <unsigned int dim> class MultiMDService;
}
} // namespace tarch

/** maps a number of MD simulations onto the total number of available ranks.
 *  For each MD simulation, a regular domain decomposition into n0 x n1 x ... x
 * nD processes is assumed.
 *  We further assume that the total number of processes can be divided by the
 * number of processes required by each MD simulation.
 *  @author Philipp Neumann
 */
template <unsigned int dim> class tarch::utils::MultiMDService {
public:
  MultiMDService(const tarch::la::Vector<dim, unsigned int>& numberProcesses, const unsigned int& totalNumberMDSimulations
#if (TARCH_PARALLEL == TARCH_YES)
                 ,
                 MPI_Comm globalComm = MPI_COMM_WORLD
#endif
  );
  ~MultiMDService();

  unsigned int getGlobalNumberOfLocalMDSimulation(unsigned int localMDSimulation) const;

  int getLocalNumberOfGlobalMDSimulation(unsigned int globalMDSimulation) const;

  unsigned int getLocalNumberOfMDSimulations() const { return _thisNumberMDSimulations; }

  unsigned int getAvgNumberOfMDSimulations() const { return _avgNumberMDSimulationsPerLocalComm; }

  tarch::la::Vector<dim, unsigned int> getNumberProcessesPerMDSimulation() const { return _numberProcessesPerMDSimulation; }

#if (TARCH_PARALLEL == TARCH_YES)
  MPI_Comm getLocalCommunicator() const { return _localComm; }
#endif
  unsigned int getLocalRank() const { return _localRank; }
  unsigned int getLocalSize() const { return _localSize; }

#if (TARCH_PARALLEL == TARCH_YES)
  MPI_Comm getGlobalCommunicator() const { return _globalComm; }
#endif
  unsigned int getGlobalRank() const { return _globalRank; }
  unsigned int getGlobalSize() const { return _globalSize; }

  // TODO: unused function
  void setTotalNumberMDSimulations(unsigned int n) { _totalNumberMDSimulations = n; }
  unsigned int getTotalNumberOfMDSimulations() { return _totalNumberMDSimulations; }

  void addMDSimulationBlock() {
    _totalNumberMDSimulations += _numberLocalComms;
    _avgNumberMDSimulationsPerLocalComm = _totalNumberMDSimulations / _numberLocalComms;
    if ((unsigned int)(_globalRank / _localSize + 1) == _numberLocalComms) {
      _thisNumberMDSimulations = _totalNumberMDSimulations - _avgNumberMDSimulationsPerLocalComm * (_numberLocalComms - 1);
    } else {
      _thisNumberMDSimulations = _avgNumberMDSimulationsPerLocalComm;
    }
  }

  void removeMDSimulationBlock() {
    _totalNumberMDSimulations -= _numberLocalComms;
    _avgNumberMDSimulationsPerLocalComm = _totalNumberMDSimulations / _numberLocalComms;
    if ((unsigned int)(_globalRank / _localSize + 1) == _numberLocalComms) {
      _thisNumberMDSimulations = _totalNumberMDSimulations - _avgNumberMDSimulationsPerLocalComm * (_numberLocalComms - 1);
    } else {
      _thisNumberMDSimulations = _avgNumberMDSimulationsPerLocalComm;
    }
  }

  int getRank() const { return this->_globalRank; }
  int getSize() const { return this->_globalSize; }

  unsigned int getNumberLocalComms() const { return _numberLocalComms; }

private:
#if (TARCH_PARALLEL == TARCH_YES)
  MPI_Comm _localComm;  // communicator of "local" MD simulation
  MPI_Comm _globalComm; // normally MPI_COMM_WORLD, unless there are ranks outside all MD simulations (e.g. due to PinT)
#endif
  // number of processes used for a single MD simulation. Currently, the total
  // number of MPI processes needs to
  // be a multiple of the number of processes per MD simulation (=product of the
  // vector components)
  const tarch::la::Vector<dim, unsigned int> _numberProcessesPerMDSimulation;
  // number of local communicators
  unsigned int _numberLocalComms;
  // total number of MD simulations
  unsigned int _totalNumberMDSimulations;
  // average number of MD simulations that is processed per local communicator.
  // If we have 4 processes per MD simulation and 12 processes available and
  // want to run 26 MD simulations,
  // then this value is given by 26/(12/4) = 26/3 = 8
  unsigned int _avgNumberMDSimulationsPerLocalComm;
  // number of MD simulations for this current local communicator. Except for
  // the "last communicator", this value
  // equals _avgNumberMDSimulationsPerLocalComm. The last comunicator group
  // fills up the missing MD simulations.
  // If we have 4 processes per MD simulation and 12 processes available and
  // want to run 26 MD simulations,
  // then we have 12/4=3 communicator groups, of which group 0 and 1 handle
  // 26/3=8 MD simulations. The last group 2
  // handles 26-2*8 = 10 MD simulations.
  unsigned int _thisNumberMDSimulations;

  int _globalSize; // global number of available MPI processes
  int _globalRank; // rank in global communicator getGlobalCommunicator()

  int _localSize; // size of communicator _localComm
  int _localRank; // local rank in communicator _localComm
};

#include "tarch/utils/MultiMDService.cpph"
#endif
