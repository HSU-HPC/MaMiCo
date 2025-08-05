// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_CELLMAPPINGS_PROFILEPLOTTERMAPPING_H_
#define _MOLECULARDYNAMICS_CELLMAPPINGS_PROFILEPLOTTERMAPPING_H_

#include "simplemd/LinkedCell.h"
#include "simplemd/services/LinkedCellService.h"
#include "simplemd/services/ParallelTopologyService.h"
#include <fstream>
#include <sstream>

namespace simplemd {
namespace cellmappings {
class ProfilePlotterMapping;
}
} // namespace simplemd

/** plots a profile for a d-dimensional subspace within the molecular dynamics
 *  simulation.
 *
 *  @author Philipp Neumann
 */
class simplemd::cellmappings::ProfilePlotterMapping {
public:
  ProfilePlotterMapping(const simplemd::services::ParallelTopologyService& parallelTopologyService,
                        const simplemd::services::LinkedCellService& linkedCellService, const unsigned int& plotEveryTimestep,
                        const unsigned int& sampleEveryTimestep, const unsigned int& startAtTimestep, const double& linkedCellVolume,
                        const unsigned int& localMDSimulation)
      : _parallelTopologyService(parallelTopologyService), _linkedCellService(linkedCellService), _plotEveryTimestep(plotEveryTimestep),
        _sampleEveryTimestep(sampleEveryTimestep), _startAtTimestep(startAtTimestep), _linkedCellVolume(linkedCellVolume),
        _localMDSimulation(localMDSimulation), _currentTimestep(0), _cellCounter(0) {
    _velocityAndDensity.clear();
  }
  ~ProfilePlotterMapping() {}

  /** sets the current timestep under consideration */
  void setCurrentTimestep(const unsigned int& t) { _currentTimestep = t; }

  void beginCellIteration() {
    if (_currentTimestep < _startAtTimestep) {
      return;
    }
    // if this is the first timestep in the sampling interval, remove old entries of the vector
    if ((_currentTimestep - _startAtTimestep) % _plotEveryTimestep == 0) {
      _velocityAndDensity.clear();
      _sampleCounter = 0;
    }

    // reset cell counter
    _cellCounter = 0;
  }

  void handleCell(LinkedCell& cell, const unsigned int& cellIndex) {
    if (_currentTimestep < _startAtTimestep) {
      return;
    }

    // if this is the first timestep in the sampling interval, create vector entry for this cell and store the cell index in vector form
    if ((_currentTimestep - _startAtTimestep) % _plotEveryTimestep == 0) {
      _velocityAndDensity.push_back(tarch::la::Vector<2 * MD_DIM + 1, double>(0.0));
      const tarch::la::Vector<MD_DIM, unsigned int> localCellIndexVector = _linkedCellService.getLocalCellIndexVector(cellIndex);
      const tarch::la::Vector<MD_DIM, unsigned int> globalCellIndexVector = _parallelTopologyService.localToGlobalCellIndexVector(localCellIndexVector);
      for (unsigned int d = 0; d < MD_DIM; d++) {
        _velocityAndDensity[_cellCounter][d] = (double)globalCellIndexVector[d];
      }
    }

    // loop over all molecules and add velocity and density contributions
    if ((_currentTimestep - _startAtTimestep) % _sampleEveryTimestep == 0) {
      tarch::la::Vector<MD_DIM, double> vel(0.0);
      double mass = 0.0;
      for (auto m1 = cell.begin(); m1 != cell.end(); m1++) {
        vel += m1->getConstVelocity();
        mass += 1.0;
      }
      // add mean velocity
      if (mass != 0.0) {
        for (unsigned int d = 0; d < MD_DIM; d++) {
          _velocityAndDensity[_cellCounter][MD_DIM + d] += vel[d] / mass;
        }
        _velocityAndDensity[_cellCounter][MD_DIM * 2] += mass;
      }
    }

    // increment cell counter
    _cellCounter++;
  }

  void endCellIteration() {
    if (_currentTimestep < _startAtTimestep) {
      return;
    }

    // increment sample counter, if we just applied the sampling
    if ((_currentTimestep - _startAtTimestep) % _sampleEveryTimestep == 0) {
      _sampleCounter++;
    }

    // if this is our last frame under consideration, write to file
    if ((_currentTimestep - _startAtTimestep) % _plotEveryTimestep == _plotEveryTimestep - 1) {
      if (_sampleCounter == 0) {
        std::cout << "ERROR simplemd::cellmappings::ProfilePlotterMapping: No samples taken!" << std::endl;
      }

      const unsigned int size = (unsigned int)_velocityAndDensity.size();
      std::stringstream ss;
      ss << "Profile_" << _localMDSimulation << "_";
#if (MD_PARALLEL == MD_YES)
      ss << _parallelTopologyService.getRank() << "_";
#endif
      ss << _currentTimestep << ".profile";
      std::ofstream file(ss.str().c_str());

      if (!file.is_open()) {
        std::cout << "ERROR simplemd::cellmappings::ProfilePlotterMapping: Could not open file " << ss.str() << "!" << std::endl;
        exit(EXIT_FAILURE);
      }

      for (unsigned int i = 0; i < size; i++) {
        // write cell index in vector form
        for (unsigned int d = 0; d < MD_DIM; d++) {
          file << _velocityAndDensity[i][d] << " ";
        }

        // write averaged velocities
        for (unsigned int d = MD_DIM; d < 2 * MD_DIM; d++) {
          file << _velocityAndDensity[i][d] / (_sampleCounter) << " ";
        }
        // write number density
        file << _velocityAndDensity[i][MD_DIM * 2] / (_linkedCellVolume * _sampleCounter) << std::endl;
      }

      file.close();
    }
  }

  static const bool IsParallel = false;

private:
  const simplemd::services::ParallelTopologyService& _parallelTopologyService;
  const simplemd::services::LinkedCellService& _linkedCellService;
  const unsigned int _plotEveryTimestep;
  const unsigned int _sampleEveryTimestep;
  const unsigned int _startAtTimestep;
  const double _linkedCellVolume;
  const unsigned int _localMDSimulation;
  unsigned int _currentTimestep;
  unsigned int _sampleCounter;
  unsigned int _cellCounter;
  /** stores:
   *  0-(MD_DIM-1) - cell-index of the linked cell in vector form
   *  MD_DIM-(2*MD_DIM-1) - velocities
   *  2*MD_DIM - number density
   */
  std::vector<tarch::la::Vector<2 * MD_DIM + 1, double>> _velocityAndDensity;
};
#endif // _MOLECULARDYNAMICS_CELLMAPPINGS_PROFILEPLOTTERMAPPING_H_
