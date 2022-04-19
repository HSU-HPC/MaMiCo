// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_CELLMAPPINGS_COMPUTEMEANVELOCITYMAPPING_H_
#define _MOLECULARDYNAMICS_CELLMAPPINGS_COMPUTEMEANVELOCITYMAPPING_H_

#include "simplemd/LinkedCell.h"
#include "simplemd/MolecularDynamicsDefinitions.h"
#include "simplemd/services/ParallelTopologyService.h"
#include <fstream>
namespace simplemd {
namespace cellmappings {
class ComputeMeanVelocityMapping;
}
} // namespace simplemd

/** computes the mean velocity over all molecules within the local domain.
 *
 *  @author Philipp Neumann
 */
class simplemd::cellmappings::ComputeMeanVelocityMapping {
public:
  ComputeMeanVelocityMapping(simplemd::services::ParallelTopologyService& parallelTopologyService, const unsigned int& localMDSimulation,
                             const bool& writeToFile = true)
      : _parallelTopologyService(parallelTopologyService), _particleCounter(0), _meanVelocity(0.0), _localMDSimulation(localMDSimulation),
        _writeToFile(writeToFile) {}
  ~ComputeMeanVelocityMapping() {}

  void beginCellIteration() {
    _particleCounter = 0;
    _meanVelocity.assign(0.0);
  }
  void endCellIteration() {
#if (MD_PARALLEL == MD_YES)
    // reduce number of particles separately from mean velocity. Datatypes differ and if the number of
    // particles is extremely high, it may be inaccurately represented by floating point numbers.
    // Hence, this solution is slower, but safer.
    _parallelTopologyService.MD_Allreduce(MPI_IN_PLACE, &_particleCounter, 1, MPI_UNSIGNED, MPI_SUM);
    _parallelTopologyService.MD_Allreduce(MPI_IN_PLACE, &_meanVelocity, MD_DIM, MPI_DOUBLE, MPI_SUM);
#endif
    _meanVelocity = (1.0 / ((double)_particleCounter)) * _meanVelocity;

    if (_writeToFile == false)
      return;

#if (MD_PARALLEL == MD_YES)
    if (_parallelTopologyService.getProcessCoordinates() == tarch::la::Vector<MD_DIM, unsigned int>(0))
#endif
    {
      std::stringstream ss;
      ss << "MeanVelocity_" << _localMDSimulation << "_" << _parallelTopologyService.getRank() << ".dat";
      std::ofstream file(ss.str().c_str(), std::ios::app);
      if (!file.is_open()) {
        std::cout << "ERROR ComputeMeanVelocityMapping: Could not open file!" << std::endl;
        exit(EXIT_FAILURE);
      }
      file << _meanVelocity << std::endl;
      file.close();
    }
  }

  void handleCell(LinkedCell& cell, const unsigned int& cellIndex) {
    for (std::list<Molecule*>::const_iterator m1 = cell.begin(); m1 != cell.end(); m1++) {
      _meanVelocity += (*m1)->getConstVelocity();
      _particleCounter++;
    }
  }

  void handleCellPair(LinkedCell& cell1, LinkedCell& cell2, const unsigned int& cellIndex1, const unsigned int& cellIndex2) {}

  const tarch::la::Vector<MD_DIM, double>& getMeanVelocity() const { return _meanVelocity; }

  /** returns the global number of particles. The local number can be retrieved from MoleculeService.getNumberMolecules(). */
  const unsigned int& getGlobalNumberMolecules() const { return _particleCounter; }

private:
  simplemd::services::ParallelTopologyService& _parallelTopologyService;
  /** number of particles.
   * After endCellIteration has been called, stores the global number of particles in the MD simulation.
   */
  unsigned int _particleCounter;
  tarch::la::Vector<MD_DIM, double> _meanVelocity;

  /** number of local MD simulation */
  const unsigned int _localMDSimulation;
  /** flag whether computed value should be written in file during endCellIteration */
  const bool _writeToFile;
};

#endif // _MOLECULARDYNAMICS_CELLMAPPINGS_COMPUTEMEANVELOCITYMAPPING_H_
