// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_MOLECULARDYNAMICSSIMULATION_H_
#define _MOLECULARDYNAMICS_MOLECULARDYNAMICSSIMULATION_H_

#include "simplemd/BoundaryTreatment.h"
#include "simplemd/MolecularDynamicsDefinitions.h"
#include "simplemd/ProfilePlotter.h"
#include "simplemd/cell-mappings/ComputeMeanVelocityMapping.h"
#include "simplemd/cell-mappings/ComputeTemperatureMapping.h"
#include "simplemd/cell-mappings/EmptyLinkedListsMapping.h"
#include "simplemd/cell-mappings/LennardJonesForceMapping.h"
#include "simplemd/cell-mappings/RDFMapping.h"
#include "simplemd/configurations/MolecularDynamicsConfiguration.h"
#include "simplemd/molecule-mappings/InitialPositionAndForceUpdate.h"
#include "simplemd/molecule-mappings/UpdateLinkedCellListsMapping.h"
#include "simplemd/molecule-mappings/VTKMoleculeWriter.h"
#include "simplemd/molecule-mappings/VelocityStoermerVerletMapping.h"
#include "simplemd/services/ExternalForceService.h"
#include "simplemd/services/LinkedCellService.h"
#include "simplemd/services/MolecularPropertiesService.h"
#include "simplemd/services/MoleculeService.h"
#include "simplemd/services/ParallelTopologyService.h"
#include "tarch/utils/MultiMDService.h"
#include "tarch/utils/RandomNumberService.h"
#include <iostream>
#include <sstream>

namespace simplemd {
class MolecularDynamicsSimulation;
}

/** steers the MD simulation.
 *
 *  @author Philipp Neumann
 */
class simplemd::MolecularDynamicsSimulation {
public:
  MolecularDynamicsSimulation(
      const simplemd::configurations::MolecularDynamicsConfiguration
          &configuration);
  virtual ~MolecularDynamicsSimulation() {}

  /** initialises all services */
  void initServices();
  /** variant for multi-MD simulations */
  void initServices(const tarch::utils::MultiMDService<MD_DIM> &multiMDService,
                    unsigned int localMDSimulation);

  /** simulates one MD timestep. It's a virtual function as we (will definitely)
   * need to extend this method when going to coupled LB-MD simulations.
   */
  void simulateOneTimestep(const unsigned int &t);

  /** runs the time loop for the MD simulation */
  void runSimulation();

  /** shuts down all services and deletes pointers */
  void shutdownServices();

  void evaluateStatistics(const unsigned int &t);

private:
  /** computes the number density for molecules given per direction and a
   * certain domain size. This is only used during initialisation.
   */
  double
  getNumberDensity(unsigned int numberMolecules,
                   const tarch::la::Vector<MD_DIM, double> &domainSize) const;

protected:
  const simplemd::configurations::MolecularDynamicsConfiguration
      &_configuration;

  // molecule mappings
  simplemd::moleculemappings::VelocityStoermerVerletMapping *_timeIntegrator;
  simplemd::moleculemappings::UpdateLinkedCellListsMapping
      *_updateLinkedCellListsMapping;
  simplemd::moleculemappings::VTKMoleculeWriter *_vtkMoleculeWriter;
  std::string _vtkFilestem;

  // cell mappings
  simplemd::cellmappings::LennardJonesForceMapping *_lennardJonesForce;
  simplemd::cellmappings::EmptyLinkedListsMapping *_emptyLinkedListsMapping;
  simplemd::cellmappings::RDFMapping *_rdfMapping;

  // boundary treatment
  simplemd::BoundaryTreatment *_boundaryTreatment;
  tarch::la::Vector<MD_LINKED_CELL_NEIGHBOURS, simplemd::BoundaryType>
      _localBoundary;

  // number of this MD simulation on this local rank
  unsigned int _localMDSimulation;
  // for plotting
  simplemd::ProfilePlotter *_profilePlotter;
  // for parallel data exchange
  simplemd::services::ParallelTopologyService *_parallelTopologyService;
  // for molecule storage
  simplemd::services::MoleculeService *_moleculeService;
  std::string _checkpointFilestem;
  // for linked cell storage
  simplemd::services::LinkedCellService *_linkedCellService;
  // molecular properties (potential parameters, mass etc)
  simplemd::services::MolecularPropertiesService *_molecularPropertiesService;
  // for external forces; has default constructor, hence we use it as object
  // instead of ptr
  simplemd::services::ExternalForceService _externalForceService;
};
#endif // _MOLECULARDYNAMICS_MOLECULARDYNAMICSSIMULATION_H_
