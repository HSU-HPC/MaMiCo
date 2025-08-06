// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_SERVICES_MOLECULESERVICE_H_
#define _MOLECULARDYNAMICS_SERVICES_MOLECULESERVICE_H_

#include "simplemd/MolecularDynamicsDefinitions.h"
#include "simplemd/Molecule.h"
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
class LinkedCellService;
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
                  const double& temperature, const unsigned int& blockSize, const simplemd::services::MolecularPropertiesService& molecularPropertiesService);

  /** initialises the MD simulation from a checkpoint-file. For a parallel simulation, this method parses
   *  checkpoint files for each rank, respectively. If multiple MD simulations are executed, make sure that the rank of the current
   *  MD simulation matches the respective rank of the checkpoint file.
   */
  MoleculeService(const tarch::la::Vector<MD_DIM, double>& domainSize, const tarch::la::Vector<MD_DIM, double>& domainOffset,
                  const std::string& checkPointFileStem, const unsigned int& blockSize,
                  const simplemd::services::ParallelTopologyService& parallelTopologyService);
  /** initialises a potentially parallel MD simulation from a sequential checkpoint file. */
  MoleculeService(const tarch::la::Vector<MD_DIM, double>& domainSize, const tarch::la::Vector<MD_DIM, double>& domainOffset,
                  const std::string& checkPointFileStem, const unsigned int& blockSize);

  /** adds a molecule to the system. The molecule data are copied from the const. reference to a free position within the memory field
   *  or - in case no memory is available - new memory is allocated and the molecule is put in there. Besides, the list _freeMoleculePositions
   *  is adapted, accordingly (and _numberMolecules is incremented).
   *  The function returns a pointer to the new molecule. If something goes wrong, NULL is returned.
   */
  Molecule* addMolecule(const Molecule& molecule);

  /** returns the number of molecules */
  const unsigned int& getNumberMolecules() const;

  /** deletes a molecule from the system by adding its position to the _freeMoleculePositions-list. */
  void deleteMolecule(Molecule& molecule);

  /** shuts down the service */
  void shutdown();


  /** creates initial velocity for molecule from meanVelocity and given temperature and stores the result in initialVelocity */
  void getInitialVelocity(const tarch::la::Vector<MD_DIM, double>& meanVelocity, const double& kB, const double& temperature,
                          const simplemd::services::MolecularPropertiesService& molecularPropertiesService,
                          tarch::la::Vector<MD_DIM, double>& initialVelocity) const;

  /** resets the velocity over the whole molecule system to the mean velocity specified at the beginning */
  void resetMeanVelocity();

  

private:

  /** pointer to all the molecules */
  std::vector<simplemd::Molecule*> _molecules;

  /** stores the mean velocity for normalisation */
  tarch::la::Vector<MD_DIM, double> _meanVelocity;

  /** number of molecules stored in memory */
  unsigned int _numberMolecules;

  /** positions within the _molecules array where a molecule can be inserted */
  std::list<unsigned int> _freeMoleculePositions;

  /** number of molecules that are stored in one memory block */
  unsigned int _blockSize;
};


#endif // _MOLECULARDYNAMICS_SERVICES_MOLECULESERVICE_H_
