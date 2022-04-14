// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_CELLMAPPINGS_COMPUTETEMPERATUREMAPPING_H_
#define _MOLECULARDYNAMICS_CELLMAPPINGS_COMPUTETEMPERATUREMAPPING_H_

#include "simplemd/LinkedCell.h"
#include "simplemd/MolecularDynamicsDefinitions.h"
#include "simplemd/services/MolecularPropertiesService.h"
#include "simplemd/services/ParallelTopologyService.h"
#include <fstream>

namespace simplemd {
namespace cellmappings {
class ComputeTemperatureMapping;
}
} // namespace simplemd

/** computes the temperature over all molecules within the local domain.
 *
 *  @author Philipp Neumann
 */
class simplemd::cellmappings::ComputeTemperatureMapping {
public:
  ComputeTemperatureMapping(simplemd::services::ParallelTopologyService& parallelTopologyService,
                            const simplemd::services::MolecularPropertiesService& molecularPropertiesService,
                            const tarch::la::Vector<MD_DIM, double>& meanVelocity, const unsigned int& localMDSimulation, const bool& writeToFile = true)
      : _parallelTopologyService(parallelTopologyService), _molecularPropertiesService(molecularPropertiesService), _meanVelocity(meanVelocity),
        _localMDSimulation(localMDSimulation), _writeToFile(writeToFile) {}
  ~ComputeTemperatureMapping() {}

  void beginCellIteration() {
    _particleCounter = 0;
    _temperature = 0.0;
  }

  void endCellIteration() {
#if (MD_PARALLEL == MD_YES)
    // reduce number of particles separately from temperature. See comment in cell-mappings/computeMeanVelocityMapping.h
    _parallelTopologyService.MD_Allreduce(MPI_IN_PLACE, &_particleCounter, 1, MPI_UNSIGNED, MPI_SUM);
    _parallelTopologyService.MD_Allreduce(MPI_IN_PLACE, &_temperature, 1, MPI_DOUBLE, MPI_SUM);
#endif
    _temperature = _temperature * _molecularPropertiesService.getMolecularProperties().getMass() /
                   (((double)MD_DIM) * _particleCounter * _molecularPropertiesService.getMolecularProperties().getKB());

    if (_writeToFile == false)
      return;

#if (MD_PARALLEL == MD_YES)
    if (_parallelTopologyService.getProcessCoordinates() == tarch::la::Vector<MD_DIM, unsigned int>(0)) {
#endif
      std::stringstream ss;
      ss << "Temperature_" << _localMDSimulation << "_" << _parallelTopologyService.getRank() << ".dat";
      std::ofstream file(ss.str().c_str(), std::ios::app);
      if (!file.is_open()) {
        std::cout << "ERROR ComputeTemperatureMapping: Could not open file!" << std::endl;
        exit(EXIT_FAILURE);
      }
      file << _temperature << std::endl;
      file.close();
#if (MD_PARALLEL == MD_YES)
    }
#endif
  }

  void handleCell(LinkedCell& cell, const unsigned int& cellIndex) {
    double buffer;
    for (std::list<Molecule*>::const_iterator m1 = cell.begin(); m1 != cell.end(); m1++) {
      buffer = tarch::la::dot(_meanVelocity - (*m1)->getConstVelocity(), _meanVelocity - (*m1)->getConstVelocity());
      _temperature += buffer;
      _particleCounter++;
    }
  }

  void handleCellPair(LinkedCell& cell1, LinkedCell& cell2, const unsigned int& cellIndex1, const unsigned int& cellIndex2) {}

  const double& getTemperature() const { return _temperature; }

private:
  simplemd::services::ParallelTopologyService& _parallelTopologyService;
  const simplemd::services::MolecularPropertiesService& _molecularPropertiesService;
  /** counts the global number of molecules */
  unsigned int _particleCounter;
  /** stores the temperature */
  double _temperature;
  /** mean velocity */
  const tarch::la::Vector<MD_DIM, double> _meanVelocity;

  /** number of local MD simulation */
  const unsigned int _localMDSimulation;
  /** flag whether computed value should be written in file during endCellIteration */
  const bool _writeToFile;
};

#endif // _MOLECULARDYNAMICS_CELLMAPPINGS_COMPUTETEMPERATUREMAPPING_H_
