// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_SERVICES_MOLECULESERVICE_H_
#define _MOLECULARDYNAMICS_SERVICES_MOLECULESERVICE_H_

#include "simplemd/MolecularDynamicsDefinitions.h"
#include "simplemd/Molecule.h"
#include "simplemd/MoleculeContainer.h"
#include "simplemd/molecule-mappings/ComputeMeanVelocityMapping.h"
#include "simplemd/molecule-mappings/SetMeanVelocityMapping.h"
#include "simplemd/services/MolecularPropertiesService.h"
#include "tarch/utils/RandomNumberService.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <list>
#include <sstream>
#include <vector>

#include <Kokkos_Core.hpp>

namespace simplemd {
namespace services {
class MoleculeService;

// forward declarations to remove circular dependencies
class ParallelTopologyService;
} // namespace services
} // namespace simplemd

/** data service storing and managing all molecules which are lying on one
 * process.
 *
 *  @author Philipp Neumann
 */
class simplemd::services::MoleculeService {
public:
  ~MoleculeService();

  /** initialises the molecules. Therefore, the molecules are put onto a regular Cartesian grid with moleculesPerDirection
   *  molecules in each spatial direction within the domain described by domainSize and domainOffset.
   *  meanVelocity describes a mean flow velocity, temperature a temperature (controlling fluctuations).
   *  In addition, the numberMoleculesPerAllocation can be set: In case that molecules are added to the system, we need to
   *  allocate more memory. If more memory is needed, a block of numberMoleculesPerAllocation is to be introduced.
   */
  MoleculeService(const tarch::la::Vector<MD_DIM, double>& domainSize, const tarch::la::Vector<MD_DIM, double>& domainOffset,
                  const tarch::la::Vector<MD_DIM, unsigned int>& moleculesPerDirection, const tarch::la::Vector<MD_DIM, double>& meanVelocity, const double& kB,
                  const double& temperature, const double capacityFactor, const simplemd::services::MolecularPropertiesService& molecularPropertiesService,
                  const simplemd::services::ParallelTopologyService& parallelTopologyService);

  /** initialises the MD simulation from a checkpoint-file. For a parallel simulation, this method parses
   *  checkpoint files for each rank, respectively. If multiple MD simulations are executed, make sure that the rank of the current
   *  MD simulation matches the respective rank of the checkpoint file.
   */
  MoleculeService(const std::string& checkPointFileStem, const double capacityFactor,
                  const simplemd::services::ParallelTopologyService& parallelTopologyService);
  /** initialises a potentially parallel MD simulation from a sequential checkpoint file. */
  MoleculeService(const tarch::la::Vector<MD_DIM, double>& domainSize, const tarch::la::Vector<MD_DIM, double>& domainOffset,
                  const std::string& checkPointFileStem, const double capacityFactor,
                  const simplemd::services::ParallelTopologyService& parallelTopologyService);

  /** returns the number of molecules */
  const unsigned int getNumberMolecules() const;

  /** shuts down the service */
  void shutdown();

  /** creates initial velocity for molecule from meanVelocity and given temperature and stores the result in initialVelocity */
  void getInitialVelocity(const tarch::la::Vector<MD_DIM, double>& meanVelocity, const double& kB, const double& temperature,
                          const simplemd::services::MolecularPropertiesService& molecularPropertiesService,
                          tarch::la::Vector<MD_DIM, double>& initialVelocity) const;

  /** writes a checkpoint containing:
   *  - the number of molecules and the dimension of the problem (1,2 or 3) in one line
   *  - each molecule in one line consisting of position, velocity and force_old.
   *  In parallel cases, each process writes its own checkpoint data. The file will be named
   *  filestem_t_rank.dat in any case (rank=0 in the serial case).
   *  The mapping WriteCheckPointMapping is used.
   */
  void writeCheckPoint(const simplemd::services::ParallelTopologyService& parallelTopologyService, const std::string& filestem, const unsigned int& t);

  /** resets the velocity over the whole molecule system to the mean velocity specified at the beginning */
  void resetMeanVelocity();

  simplemd::MoleculeContainer& getContainer() const { return *_moleculeContainer; }

  tarch::la::Vector<MD_DIM, double> getLocalDomainSize() { return _localDomainSize; }

  static bool tarchDebugIsOn();

private:
  void initContainer(ParallelTopologyService parallelTopologyService, size_t moleculeCount, double capacityFactor);

  /** stores the mean velocity for normalisation */
  tarch::la::Vector<MD_DIM, double> _meanVelocity;

  /** stores the spatial extent of the local domain */
  const tarch::la::Vector<MD_DIM, double> _localDomainSize;

  simplemd::MoleculeContainer* _moleculeContainer;
};

#endif // _MOLECULARDYNAMICS_SERVICES_MOLECULESERVICE_H_
