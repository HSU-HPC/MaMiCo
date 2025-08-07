// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_SERVICES_LINKEDCELLSERVICE_H_
#define _MOLECULARDYNAMICS_SERVICES_LINKEDCELLSERVICE_H_

#include "simplemd/LinkedCell.h"
#include "simplemd/Molecule.h"
#include "simplemd/services/ParallelTopologyService.h"
#include "simplemd/services/MoleculeService.h"
#include "tarch/la/Vector.h"
#include <iostream>

#include <Kokkos_Core.hpp>

namespace simplemd {
namespace services {
class LinkedCellService;
}
} // namespace simplemd

/** manages the linked cell data structure of the simulation.
 *  @author Philipp Neumann
 */
class simplemd::services::LinkedCellService {
public:
  /** initialises the linked cell service:
   *  domainSize - size of local domain
   *  domainOffset - starting coordinate of local domain (lowerLeftFront point)
   *  numberOfCells - local number of cells
   */
  LinkedCellService(const tarch::la::Vector<MD_DIM, double>& domainSize, const tarch::la::Vector<MD_DIM, double>& domainOffset,
                    const simplemd::services::ParallelTopologyService& parallelTopologyService, simplemd::services::MoleculeService& moleculeService);

  /** shuts down the service, frees memory and resets all variables */
  void shutdown();

  /** puts the molecule into the cell defined by the local index (vector)forward declarations to circumvent circular includes
   * coordinates localCellIndex */
  void addMoleculeToLinkedCell(Molecule& molecule, const tarch::la::Vector<MD_DIM, unsigned int>& localCellIndex);
  /** puts the molecule into the cell defined by the local index (scalar)
   * coordinates localCellIndex */
  void addMoleculeToLinkedCell(Molecule& molecule, const unsigned int& localCellIndex);

  void deleteMoleculeFromLinkedCell(Molecule& molecule, const tarch::la::Vector<MD_DIM, unsigned int>& localCellIndex);

  /** returns the linked cell at the respective coordinates */
  LinkedCell& getLinkedCell(const tarch::la::Vector<MD_DIM, unsigned int>& localCellIndex);

  // TODO REMOVE BELOW
  /** returns the index of the first (non-ghost) cell */
  const tarch::la::Vector<MD_DIM, unsigned int>& getLocalIndexOfFirstCell() const { return *((tarch::la::Vector<MD_DIM, unsigned int>*)NULL); };

  /** returns the number of (non-ghost) cells */
  const tarch::la::Vector<MD_DIM, unsigned int>& getLocalNumberOfCells() const { return *((tarch::la::Vector<MD_DIM, unsigned int>*)NULL); };
  // TODO REMOVE ABOVE

  /** returns the mesh width */
  const tarch::la::Vector<MD_DIM, double>& getMeshWidth() const;

  /** returns the local domain offset (for the domain of this process) */
  const tarch::la::Vector<MD_DIM, double>& getLocalDomainOffset() const;

  /** returns the local domain size (for the domain of this process) */
  const tarch::la::Vector<MD_DIM, double>& getLocalDomainSize() const;

  /** returns the local cell index from the local cell index vector */
  unsigned int getLocalCellIndex(const tarch::la::Vector<MD_DIM, unsigned int>& cellIndexVector) const;

  ~LinkedCellService() {
    if (_cells != NULL) {
      delete[] _cells;
      _cells = NULL;
    }
  }

private:
  /** initialise linked-cell structure for local process.
   *  indexOffset denotes the integer coordinates of the first cell
   *  within the local cell structure; the grid has a total of numberOfCells
   * cells. globalIndexFirstCell denotes the global index of the lower left cell
   * of the simulation. This is important for parallel computations, only.
   */
  void initCellStructure();

  /** returns local index from (local) coordinate vector */
  unsigned int getLocalIndexFromLocalVector(const tarch::la::Vector<MD_DIM, unsigned int>& coords) const;

  /** computes the mesh width from domain size and local number of grid cells */
  tarch::la::Vector<MD_DIM, double> getMeshwidth(const tarch::la::Vector<MD_DIM, double>& domainSize,
                                                 const tarch::la::Vector<MD_DIM, unsigned int>& localNumberCells) const {
    tarch::la::Vector<MD_DIM, double> meshWidth(0.0);
    for (unsigned int d = 0; d < MD_DIM; d++) {
      meshWidth[d] = domainSize[d] / localNumberCells[d];
    }
    return meshWidth;
  }

  /** contains all (local) linked cells */
  LinkedCell* _cells;
  /** size of global domain */
  const tarch::la::Vector<MD_DIM, double> _domainSize;
  /** offset of local domain */
  const tarch::la::Vector<MD_DIM, double> _domainOffset;
  /** mesh width of linked cells */
  const tarch::la::Vector<MD_DIM, double> _meshWidth;
  /** number of cells of local domain, without ghost layer
   */
  const tarch::la::Vector<MD_DIM, unsigned int> _numberOfCells;
  /** index of first cell under consideration. It is 1,1,1, due to a ghost cell
   * layer around the domain.
   */
  const tarch::la::Vector<MD_DIM, unsigned int> _indexOffset;
  /** number of cells of local domain, including ghost layer */
  const tarch::la::Vector<MD_DIM, unsigned int> _totalNumberOfCells;
/** _totalNumberOfCells(0)*_totalNumberOfCells(1); only stored for performance
 * reasons */
#if (MD_DIM > 2)
  const unsigned int _totalNumberOfCells_X_By_totalNumberOfCells_Y;
#endif
};

#endif // _MOLECULARDYNAMICS_SERVICES_LINKEDCELLSERVICE_H_
