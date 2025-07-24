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
#if (MD_OPENMP == MD_YES)
#include <omp.h>
#endif

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

  /** can be used to apply a molecule-mapping which is iterated over all molecules of this process
   *  (e.g. time integration)
   */
  template <class A> void iterateMolecules(A& a, const bool& useOpenMP);

  /** creates initial velocity for molecule from meanVelocity and given temperature and stores the result in initialVelocity */
  void getInitialVelocity(const tarch::la::Vector<MD_DIM, double>& meanVelocity, const double& kB, const double& temperature,
                          const simplemd::services::MolecularPropertiesService& molecularPropertiesService,
                          tarch::la::Vector<MD_DIM, double>& initialVelocity) const;

  /** reorganises the storage of the molecules. If a simulation requires a big number of
   *  molecule deletions/ insertions, e.g., this might be useful to speed up the simulation.
   *  Besides, the molecules are stored in memory such that those molecules belonging to the same linked cell are
   *  located very closely in memory (actually, we sort the molecules according to their linked cell position, that is
   *  lexicographically w.r.t. to the linked cell index, and store them in this sequence).
   */
  void reorganiseMemory(const simplemd::services::ParallelTopologyService& parallelTopologyService, simplemd::services::LinkedCellService& linkedCellService);

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

#pragma region TODO replace with Kokkos datastructure
  Molecule* getCellMolecule(const unsigned int cellIndex, const unsigned int moleculeIndex);
  void addCellMolecule(Molecule& molecule, const unsigned int cellIndex);
  void clearCellMolecules(const unsigned int cellIndex);

private:
  // TODO (temporary): Will be replaced by datastructure wrapping Kokkos::View<simplemd::Molecule**>
  std::vector<std::vector<simplemd::Molecule*>> _linkedCellsMolecules;
#pragma endregion TODO replace with Kokkos datastructure

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

template <class A> void simplemd::services::MoleculeService::iterateMolecules(A& a, const bool& useOpenMP) {

  const unsigned int blockSize = _blockSize;
  const unsigned int freeMoleculePositions = (const unsigned int)_freeMoleculePositions.size();
  const unsigned int numberMolecules = _numberMolecules;
  const unsigned int freeMoleculePositionsAndNumberMolecules = (const unsigned int)(numberMolecules + freeMoleculePositions);
  // start iteration();
  a.beginMoleculeIteration();

// open MP
#if (MD_OPENMP == MD_YES)
  if (useOpenMP) {

    // sort empty positions list
    if (!_freeMoleculePositions.empty()) {
      _freeMoleculePositions.sort();
      std::list<unsigned int>::iterator myIt = _freeMoleculePositions.begin();
      unsigned int start = 0;

      // loop over all intervals, starting at a certain point and ranging up to a deleted position
      for (unsigned int i = 0; i < freeMoleculePositions; i++) {
        const unsigned int end = (*myIt);
#pragma omp parallel for
        for (unsigned int j = start; j < end; j++) {
#if (MD_DEBUG == MD_YES)
          std::cout << "Handle molecule " << j << std::endl;
#endif
          a.handleMolecule(_molecules[j / blockSize][j % blockSize]);
        }

        // go to next possible start position (one position after *myIt) and increment myIt
        // -> remark: It may happen that *(myIt++) == (*myIt)+1. Then, the upper inner loop degenerates to
        // an empty loop...
        start = end + 1;
        myIt++;
      }

// do final loop (from last deleted molecule to last existing molecule)
#pragma omp parallel for
      for (unsigned int i = start; i < freeMoleculePositionsAndNumberMolecules; i++) {
#if (MD_DEBUG == MD_YES)
        std::cout << "Handle molecule " << i << std::endl;
#endif
        a.handleMolecule(_molecules[i / blockSize][i % blockSize]);
      }
    } else {
#pragma omp parallel for
      for (unsigned int i = 0; i < numberMolecules; i++) {
#if (MD_DEBUG == MD_YES)
        std::cout << "Handle molecule " << i << std::endl;
#endif
        a.handleMolecule(_molecules[i / blockSize][i % blockSize]);
      }
    }

    // no Open MP
  } else {
#endif

    // sort empty positions list
    if (!_freeMoleculePositions.empty()) {
      _freeMoleculePositions.sort();
      std::list<unsigned int>::iterator myIt = _freeMoleculePositions.begin();
      unsigned int start = 0;

      // loop over all intervals, starting at a certain point and ranging up to a deleted position
      for (unsigned int i = 0; i < freeMoleculePositions; i++) {
        const unsigned int end = *myIt;
        for (unsigned int j = start; j < end; j++) {
#if (MD_DEBUG == MD_YES)
          std::cout << "Handle molecule " << j << std::endl;
#endif
          a.handleMolecule(_molecules[j / blockSize][j % blockSize]);
        }

        // go to next possible start position (one position after *myIt) and increment myIt
        // -> remark: It may happen that *(myIt++) == (*myIt)+1. Then, the upper inner loop degenerates to
        // an empty loop...
        start = (*myIt) + 1;
        myIt++;
      }

      // do final loop (from last deleted molecule to last existing molecule)
      for (unsigned int i = start; i < freeMoleculePositionsAndNumberMolecules; i++) {
#if (MD_DEBUG == MD_YES)
        std::cout << "Handle molecule " << i << std::endl;
#endif
        a.handleMolecule(_molecules[i / blockSize][i % blockSize]);
      }
    } else {
      for (unsigned int i = 0; i < numberMolecules; i++) {
#if (MD_DEBUG == MD_YES)
        std::cout << "Handle molecule " << i << std::endl;
#endif
        a.handleMolecule(_molecules[i / blockSize][i % blockSize]);
      }
    }

#if (MD_OPENMP == MD_YES)
  }
#endif

  // end iteration();
  a.endMoleculeIteration();
}

#endif // _MOLECULARDYNAMICS_SERVICES_MOLECULESERVICE_H_
