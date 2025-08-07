// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_BOUNDARYTREATMENT_H_
#define _MOLECULARDYNAMICS_BOUNDARYTREATMENT_H_

#include "simplemd/cell-mappings/DeleteMoleculesMapping.h"
#include "simplemd/cell-mappings/LennardJonesForceMapping.h"
#include "simplemd/cell-mappings/PeriodicAndParallelBoundaryFillCellsMapping.h"
#include "simplemd/cell-mappings/PeriodicBoundaryEmptyCellsMapping.h"
#if (MD_PARALLEL == MD_YES)
#include "simplemd/cell-mappings/ParallelBoundaryEmptyCellsMapping.h"
#endif
#include "simplemd/MolecularDynamicsDefinitions.h"
#include "simplemd/services/LinkedCellService.h"
#include "simplemd/services/ParallelTopologyService.h"
#include "simplemd/MoleculeContainer.h"

/** This file comprises some functions that are triggered in order to have consistent
 *  boundary cells in the ghost layer.
 *  @author Philipp Neumann
 */

namespace simplemd {
class BoundaryTreatment;
}

class simplemd::BoundaryTreatment {
public:
  BoundaryTreatment(simplemd::services::ParallelTopologyService& parallelTopologyService, simplemd::MoleculeContainer& moleculeContainer,
                    simplemd::services::LinkedCellService& linkedCellService)
      : _moleculeContainer(moleculeContainer), _linkedCellService(linkedCellService), _periodicBoundaryMapping(parallelTopologyService, moleculeContainer),
        _deleteMoleculesMapping(),
#if (MD_PARALLEL == MD_YES)
        _parallelBoundaryMapping(parallelTopologyService, moleculeContainer),
#endif
        _fillCellsMapping(parallelTopologyService, moleculeContainer, linkedCellService) {
  }
  ~BoundaryTreatment() {}

  /** remove all molecules from ghost cells. This function is triggered for all ghost cells
   *  right after the force computation, since afterwards, the boundary particles
   *  are not needed anymore.
   */
  void emptyGhostBoundaryCells();

  /** remove all molecules from ghost cells which belong to periodic/ parallel conditions and
   *  put these molecules into the respective inner cell. This function is typically
   *  triggered after the time integration. When advancing in time, some molecules
   *  might enter the ghost layer. So, we need to send those molecules back into the
   *  original domain.
   */
  void putBoundaryParticlesToInnerCells(const tarch::la::Vector<MD_LINKED_CELL_NEIGHBOURS, simplemd::BoundaryType>& boundary,
                                        simplemd::services::ParallelTopologyService& parallelTopologyService);

  /** take all molecules from an inner cell and copy them to each periodic/ parallel ghost cell for the next
   *  timestep. This function is triggered right before the force evaluation between all particle
   *  pairs.
   */
  void fillBoundaryCells(const tarch::la::Vector<MD_LINKED_CELL_NEIGHBOURS, simplemd::BoundaryType>& boundary,
                         simplemd::services::ParallelTopologyService& parallelTopologyService);

  /** combined version of the two forementioned functions. Needed to reduce communication
   *  calls in half in the parallel case.
   *  @see putBoundaryParticlesToInnerCells
   *  @see fillBoundaryCells
   */
  void putBoundaryParticlesToInnerCellsAndFillBoundaryCells(const tarch::la::Vector<MD_LINKED_CELL_NEIGHBOURS, simplemd::BoundaryType>& boundary,
                                                            simplemd::services::ParallelTopologyService& parallelTopologyService);

  /** overlaps waiting for communication requests to be fulfilled with force computations
   *  on inner part of domain.
   *  @see putBoundaryParticlesToInnerCellsAndFillBoundaryCells
   */
  void putBoundaryParticlesToInnerCellsFillBoundaryCellsAndOverlapWithForceComputations(
      const tarch::la::Vector<MD_LINKED_CELL_NEIGHBOURS, simplemd::BoundaryType>& boundary,
      simplemd::services::ParallelTopologyService& parallelTopologyService, simplemd::cellmappings::LennardJonesForceMapping& lennardJonesForce);

private:
  /** applies the mapping myMapping to all boundaries of the domain which are of type
   *  boundaryType. Here, all cells in the respective ghost layer are traversed.
   */
  template <class Mapping>
  void applyMappingToBoundaryCells(const tarch::la::Vector<MD_LINKED_CELL_NEIGHBOURS, simplemd::BoundaryType>& boundary,
                                   const simplemd::BoundaryType boundaryType, Mapping& myMapping) const;

  /** applies the mapping myMapping to all cells which lie on the outermost layer of inner cells
   */
  template <class Mapping> void applyMappingToOutermostNonBoundaryCells(Mapping& myMapping) const;

  /** applies the mapping myMapping to all cells which are not directly influenced by communication
   *  at the current iteration
   */
  template <class Mapping> void applyMappingToCommunicationIndependentCells(Mapping& myMapping) const;

  /** applies the mapping myMapping to all cells which are directly influenced by communication
   */
  template <class Mapping> void applyMappingToCommunicationDependentCells(Mapping& myMapping) const;

  simplemd::MoleculeContainer& _moleculeContainer;
  simplemd::services::LinkedCellService& _linkedCellService;
  simplemd::cellmappings::PeriodicBoundaryEmptyCellsMapping _periodicBoundaryMapping;
  simplemd::cellmappings::DeleteMoleculesMapping _deleteMoleculesMapping;
#if (MD_PARALLEL == MD_YES)
  simplemd::cellmappings::ParallelBoundaryEmptyCellsMapping _parallelBoundaryMapping;
#endif
  simplemd::cellmappings::PeriodicAndParallelBoundaryFillCellsMapping _fillCellsMapping;
};

template <class Mapping>
void simplemd::BoundaryTreatment::applyMappingToBoundaryCells(const tarch::la::Vector<MD_LINKED_CELL_NEIGHBOURS, simplemd::BoundaryType>& boundary,
                                                              const simplemd::BoundaryType boundaryType, Mapping& myMapping) const {
  tarch::la::Vector<MD_DIM, unsigned int> startOuter;
  tarch::la::Vector<MD_DIM, unsigned int> numberCellsOuter;

#if (MD_DIM == 1)
  if (boundary[0] == boundaryType) {
    // left outer corner
    numberCellsOuter[0] = 1;
    startOuter[0] = 0;
    _moleculeContainer.iterateCells(myMapping, startOuter, numberCellsOuter);
  }
  if (boundary[1] == boundaryType) {
    // right outer corner
    numberCellsOuter[0] = 1;
    startOuter[0] = _moleculeContainer.getLocalNumberOfCells()[0] + 2 * _moleculeContainer.getLocalIndexOfFirstCell()[0] - 1;
    _moleculeContainer.iterateCells(myMapping, startOuter, numberCellsOuter);
  }
#endif

#if (MD_DIM == 2)
  // corners -----------------------------------------------
  if (boundary[0] == boundaryType) {
    // left front corner
    numberCellsOuter[0] = 1;
    numberCellsOuter[1] = 1;
    startOuter[0] = 0;
    startOuter[1] = 0;
    _moleculeContainer.iterateCells(myMapping, startOuter, numberCellsOuter);
  }
  if (boundary[7] == boundaryType) {
    // right back corner
    numberCellsOuter[0] = 1;
    numberCellsOuter[1] = 1;
    startOuter[0] = _moleculeContainer.getLocalNumberOfCells()[0] + 2 * _moleculeContainer.getLocalIndexOfFirstCell()[0] - 1;
    startOuter[1] = _moleculeContainer.getLocalNumberOfCells()[1] + 2 * _moleculeContainer.getLocalIndexOfFirstCell()[1] - 1;
    _moleculeContainer.iterateCells(myMapping, startOuter, numberCellsOuter);
  }
  if (boundary[2] == boundaryType) {
    // right front corner
    numberCellsOuter[0] = 1;
    numberCellsOuter[1] = 1;
    startOuter[0] = _moleculeContainer.getLocalNumberOfCells()[0] + 2 * _moleculeContainer.getLocalIndexOfFirstCell()[0] - 1;
    startOuter[1] = 0;
    _moleculeContainer.iterateCells(myMapping, startOuter, numberCellsOuter);
  }
  if (boundary[5] == boundaryType) {
    // left back corner
    numberCellsOuter[0] = 1;
    numberCellsOuter[1] = 1;
    startOuter[0] = 0;
    startOuter[1] = _moleculeContainer.getLocalNumberOfCells()[1] + 2 * _moleculeContainer.getLocalIndexOfFirstCell()[1] - 1;
    _moleculeContainer.iterateCells(myMapping, startOuter, numberCellsOuter);
  }
  // edges -------------------------------------------------
  if (boundary[1] == boundaryType) {
    // front edge
    numberCellsOuter[0] = _moleculeContainer.getLocalNumberOfCells()[0];
    numberCellsOuter[1] = 1;
    startOuter[0] = 1;
    startOuter[1] = 0;
    _moleculeContainer.iterateCells(myMapping, startOuter, numberCellsOuter);
  }
  if (boundary[6] == boundaryType) {
    // back edge
    numberCellsOuter[0] = _moleculeContainer.getLocalNumberOfCells()[0];
    numberCellsOuter[1] = 1;
    startOuter[0] = 1;
    startOuter[1] = _moleculeContainer.getLocalNumberOfCells()[1] + 2 * _moleculeContainer.getLocalIndexOfFirstCell()[1] - 1;
    _moleculeContainer.iterateCells(myMapping, startOuter, numberCellsOuter);
  }
  if (boundary[3] == boundaryType) {
    // left edge
    numberCellsOuter[0] = 1;
    numberCellsOuter[1] = _moleculeContainer.getLocalNumberOfCells()[1];
    startOuter[0] = 0;
    startOuter[1] = 1;
    _moleculeContainer.iterateCells(myMapping, startOuter, numberCellsOuter);
  }
  if (boundary[4] == boundaryType) {
    // right edge
    numberCellsOuter[0] = 1;
    numberCellsOuter[1] = _moleculeContainer.getLocalNumberOfCells()[1];
    startOuter[0] = _moleculeContainer.getLocalNumberOfCells()[0] + 2 * _moleculeContainer.getLocalIndexOfFirstCell()[0] - 1;
    startOuter[1] = 1;
    _moleculeContainer.iterateCells(myMapping, startOuter, numberCellsOuter);
  }
#endif

#if (MD_DIM == 3)
  // corners ---------------------------------------------
  if (boundary[0] == boundaryType) {
    // lower,left,front corner
    numberCellsOuter = tarch::la::Vector<MD_DIM, unsigned int>(1);
    startOuter = tarch::la::Vector<MD_DIM, unsigned int>(0);
    _moleculeContainer.iterateCells(myMapping, startOuter, numberCellsOuter);
  }
  if (boundary[25] == boundaryType) {
    // upper,right,back corner
    numberCellsOuter = tarch::la::Vector<MD_DIM, unsigned int>(1);
    startOuter = _moleculeContainer.getLocalNumberOfCells() + 2u * _moleculeContainer.getLocalIndexOfFirstCell() - tarch::la::Vector<MD_DIM, unsigned int>(1);
    _moleculeContainer.iterateCells(myMapping, startOuter, numberCellsOuter);
  }
  if (boundary[2] == boundaryType) {
    // lower,right,front corner
    numberCellsOuter = tarch::la::Vector<MD_DIM, unsigned int>(1);
    startOuter[0] = _moleculeContainer.getLocalNumberOfCells()[0] + 2 * _moleculeContainer.getLocalIndexOfFirstCell()[0] - 1;
    startOuter[1] = 0;
    startOuter[2] = 0;
    _moleculeContainer.iterateCells(myMapping, startOuter, numberCellsOuter);
  }
  if (boundary[23] == boundaryType) {
    // upper,left,back corner
    numberCellsOuter = tarch::la::Vector<MD_DIM, unsigned int>(1);
    startOuter[0] = 0;
    startOuter[1] = _moleculeContainer.getLocalNumberOfCells()[1] + 2 * _moleculeContainer.getLocalIndexOfFirstCell()[1] - 1;
    startOuter[2] = _moleculeContainer.getLocalNumberOfCells()[2] + 2 * _moleculeContainer.getLocalIndexOfFirstCell()[2] - 1;
    _moleculeContainer.iterateCells(myMapping, startOuter, numberCellsOuter);
  }
  if (boundary[6] == boundaryType) {
    // lower,left,back corner
    numberCellsOuter = tarch::la::Vector<MD_DIM, unsigned int>(1);
    startOuter[0] = 0;
    startOuter[1] = _moleculeContainer.getLocalNumberOfCells()[1] + 2 * _moleculeContainer.getLocalIndexOfFirstCell()[1] - 1;
    startOuter[2] = 0;
    _moleculeContainer.iterateCells(myMapping, startOuter, numberCellsOuter);
  }
  if (boundary[19] == boundaryType) {
    // upper,right,front corner
    numberCellsOuter = tarch::la::Vector<MD_DIM, unsigned int>(1);
    startOuter[0] = _moleculeContainer.getLocalNumberOfCells()[0] + 2 * _moleculeContainer.getLocalIndexOfFirstCell()[0] - 1;
    startOuter[1] = 0;
    startOuter[2] = _moleculeContainer.getLocalNumberOfCells()[2] + 2 * _moleculeContainer.getLocalIndexOfFirstCell()[2] - 1;
    _moleculeContainer.iterateCells(myMapping, startOuter, numberCellsOuter);
  }
  if (boundary[8] == boundaryType) {
    // lower,right,back corner
    numberCellsOuter = tarch::la::Vector<MD_DIM, unsigned int>(1);
    startOuter[0] = _moleculeContainer.getLocalNumberOfCells()[0] + 2 * _moleculeContainer.getLocalIndexOfFirstCell()[0] - 1;
    startOuter[1] = _moleculeContainer.getLocalNumberOfCells()[1] + 2 * _moleculeContainer.getLocalIndexOfFirstCell()[1] - 1;
    startOuter[2] = 0;
    _moleculeContainer.iterateCells(myMapping, startOuter, numberCellsOuter);
  }
  if (boundary[17] == boundaryType) {
    // upper,left,front corner
    numberCellsOuter = tarch::la::Vector<MD_DIM, unsigned int>(1);
    startOuter[0] = 0;
    startOuter[1] = 0;
    startOuter[2] = _moleculeContainer.getLocalNumberOfCells()[2] + 2 * _moleculeContainer.getLocalIndexOfFirstCell()[2] - 1;
    _moleculeContainer.iterateCells(myMapping, startOuter, numberCellsOuter);
  }

  // edges -------------------------------------------
  if (boundary[1] == boundaryType) {
    // x-axis aligned: lower,front edge
    numberCellsOuter[0] = _moleculeContainer.getLocalNumberOfCells()[0];
    numberCellsOuter[1] = 1;
    numberCellsOuter[2] = 1;
    startOuter[0] = _moleculeContainer.getLocalIndexOfFirstCell()[0];
    startOuter[1] = 0;
    startOuter[2] = 0;
    _moleculeContainer.iterateCells(myMapping, startOuter, numberCellsOuter);
  }
  if (boundary[24] == boundaryType) {
    // upper,back edge
    numberCellsOuter[0] = _moleculeContainer.getLocalNumberOfCells()[0];
    numberCellsOuter[1] = 1;
    numberCellsOuter[2] = 1;
    startOuter[0] = _moleculeContainer.getLocalIndexOfFirstCell()[0];
    startOuter[1] = _moleculeContainer.getLocalNumberOfCells()[1] + 2 * _moleculeContainer.getLocalIndexOfFirstCell()[1] - 1;
    startOuter[2] = _moleculeContainer.getLocalNumberOfCells()[2] + 2 * _moleculeContainer.getLocalIndexOfFirstCell()[2] - 1;
    _moleculeContainer.iterateCells(myMapping, startOuter, numberCellsOuter);
  }
  if (boundary[7] == boundaryType) {
    // x-axis aligned: lower,back edge
    numberCellsOuter[0] = _moleculeContainer.getLocalNumberOfCells()[0];
    numberCellsOuter[1] = 1;
    numberCellsOuter[2] = 1;
    startOuter[0] = _moleculeContainer.getLocalIndexOfFirstCell()[0];
    startOuter[1] = _moleculeContainer.getLocalNumberOfCells()[1] + 2 * _moleculeContainer.getLocalIndexOfFirstCell()[1] - 1;
    startOuter[2] = 0;
    _moleculeContainer.iterateCells(myMapping, startOuter, numberCellsOuter);
  }
  if (boundary[18] == boundaryType) {
    // upper,front edge
    numberCellsOuter[0] = _moleculeContainer.getLocalNumberOfCells()[0];
    numberCellsOuter[1] = 1;
    numberCellsOuter[2] = 1;
    startOuter[0] = _moleculeContainer.getLocalIndexOfFirstCell()[0];
    startOuter[1] = 0;
    startOuter[2] = _moleculeContainer.getLocalNumberOfCells()[2] + 2 * _moleculeContainer.getLocalIndexOfFirstCell()[2] - 1;
    _moleculeContainer.iterateCells(myMapping, startOuter, numberCellsOuter);
  }
  if (boundary[3] == boundaryType) {
    // y-axis aligned: lower,left edge
    numberCellsOuter[0] = 1;
    numberCellsOuter[1] = _moleculeContainer.getLocalNumberOfCells()[1];
    numberCellsOuter[2] = 1;
    startOuter[0] = 0;
    startOuter[1] = _moleculeContainer.getLocalIndexOfFirstCell()[1];
    startOuter[2] = 0;
    _moleculeContainer.iterateCells(myMapping, startOuter, numberCellsOuter);
  }
  if (boundary[22] == boundaryType) {
    // upper,right edge
    numberCellsOuter[0] = 1;
    numberCellsOuter[1] = _moleculeContainer.getLocalNumberOfCells()[1];
    numberCellsOuter[2] = 1;
    startOuter[0] = _moleculeContainer.getLocalNumberOfCells()[0] + 2 * _moleculeContainer.getLocalIndexOfFirstCell()[0] - 1;
    startOuter[1] = _moleculeContainer.getLocalIndexOfFirstCell()[1];
    startOuter[2] = _moleculeContainer.getLocalNumberOfCells()[2] + 2 * _moleculeContainer.getLocalIndexOfFirstCell()[2] - 1;
    _moleculeContainer.iterateCells(myMapping, startOuter, numberCellsOuter);
  }
  if (boundary[5] == boundaryType) {
    // y-axis aligned: lower,right edge
    numberCellsOuter[0] = 1;
    numberCellsOuter[1] = _moleculeContainer.getLocalNumberOfCells()[1];
    numberCellsOuter[2] = 1;
    startOuter[0] = _moleculeContainer.getLocalNumberOfCells()[0] + 2 * _moleculeContainer.getLocalIndexOfFirstCell()[0] - 1;
    startOuter[1] = _moleculeContainer.getLocalIndexOfFirstCell()[1];
    startOuter[2] = 0;
    _moleculeContainer.iterateCells(myMapping, startOuter, numberCellsOuter);
  }
  if (boundary[20] == boundaryType) {
    // upper,left edge
    numberCellsOuter[0] = 1;
    numberCellsOuter[1] = _moleculeContainer.getLocalNumberOfCells()[1];
    numberCellsOuter[2] = 1;
    startOuter[0] = 0;
    startOuter[1] = _moleculeContainer.getLocalIndexOfFirstCell()[1];
    startOuter[2] = _moleculeContainer.getLocalNumberOfCells()[2] + 2 * _moleculeContainer.getLocalIndexOfFirstCell()[2] - 1;
    _moleculeContainer.iterateCells(myMapping, startOuter, numberCellsOuter);
  }
  if (boundary[9] == boundaryType) {
    // z-axis aligned: left,front edge
    numberCellsOuter[0] = 1;
    numberCellsOuter[1] = 1;
    numberCellsOuter[2] = _moleculeContainer.getLocalNumberOfCells()[2];
    startOuter[0] = 0;
    startOuter[1] = 0;
    startOuter[2] = _moleculeContainer.getLocalIndexOfFirstCell()[2];
    _moleculeContainer.iterateCells(myMapping, startOuter, numberCellsOuter);
  }
  if (boundary[16] == boundaryType) {
    // right,back edge
    numberCellsOuter[0] = 1;
    numberCellsOuter[1] = 1;
    numberCellsOuter[2] = _moleculeContainer.getLocalNumberOfCells()[2];
    startOuter[0] = _moleculeContainer.getLocalNumberOfCells()[0] + 2 * _moleculeContainer.getLocalIndexOfFirstCell()[0] - 1;
    startOuter[1] = _moleculeContainer.getLocalNumberOfCells()[1] + 2 * _moleculeContainer.getLocalIndexOfFirstCell()[1] - 1;
    startOuter[2] = _moleculeContainer.getLocalIndexOfFirstCell()[2];
    _moleculeContainer.iterateCells(myMapping, startOuter, numberCellsOuter);
  }
  if (boundary[11] == boundaryType) {
    // z-axis aligned: right,front edge
    numberCellsOuter[0] = 1;
    numberCellsOuter[1] = 1;
    numberCellsOuter[2] = _moleculeContainer.getLocalNumberOfCells()[2];
    startOuter[0] = _moleculeContainer.getLocalNumberOfCells()[0] + 2 * _moleculeContainer.getLocalIndexOfFirstCell()[0] - 1;
    startOuter[1] = 0;
    startOuter[2] = _moleculeContainer.getLocalIndexOfFirstCell()[2];
    _moleculeContainer.iterateCells(myMapping, startOuter, numberCellsOuter);
  }
  if (boundary[14] == boundaryType) {
    // left,back edge
    numberCellsOuter[0] = 1;
    numberCellsOuter[1] = 1;
    numberCellsOuter[2] = _moleculeContainer.getLocalNumberOfCells()[2];
    startOuter[0] = 0;
    startOuter[1] = _moleculeContainer.getLocalNumberOfCells()[1] + 2 * _moleculeContainer.getLocalIndexOfFirstCell()[1] - 1;
    startOuter[2] = _moleculeContainer.getLocalIndexOfFirstCell()[2];
    _moleculeContainer.iterateCells(myMapping, startOuter, numberCellsOuter);
  }
  // faces ---------------------------------------
  if (boundary[4] == boundaryType) {
    // bottom face
    numberCellsOuter[0] = _moleculeContainer.getLocalNumberOfCells()[0];
    numberCellsOuter[1] = _moleculeContainer.getLocalNumberOfCells()[1];
    numberCellsOuter[2] = 1;
    startOuter[0] = _moleculeContainer.getLocalIndexOfFirstCell()[0];
    startOuter[1] = _moleculeContainer.getLocalIndexOfFirstCell()[1];
    startOuter[2] = 0;
    _moleculeContainer.iterateCells(myMapping, startOuter, numberCellsOuter);
  }
  if (boundary[21] == boundaryType) {
    // top face
    numberCellsOuter[0] = _moleculeContainer.getLocalNumberOfCells()[0];
    numberCellsOuter[1] = _moleculeContainer.getLocalNumberOfCells()[1];
    numberCellsOuter[2] = 1;
    startOuter[0] = _moleculeContainer.getLocalIndexOfFirstCell()[0];
    startOuter[1] = _moleculeContainer.getLocalIndexOfFirstCell()[1];
    startOuter[2] = _moleculeContainer.getLocalNumberOfCells()[2] + 2 * _moleculeContainer.getLocalIndexOfFirstCell()[2] - 1;
    _moleculeContainer.iterateCells(myMapping, startOuter, numberCellsOuter);
  }
  if (boundary[12] == boundaryType) {
    // left face
    numberCellsOuter[0] = 1;
    numberCellsOuter[1] = _moleculeContainer.getLocalNumberOfCells()[1];
    numberCellsOuter[2] = _moleculeContainer.getLocalNumberOfCells()[2];
    startOuter[0] = 0;
    startOuter[1] = _moleculeContainer.getLocalIndexOfFirstCell()[1];
    startOuter[2] = _moleculeContainer.getLocalIndexOfFirstCell()[2];
    _moleculeContainer.iterateCells(myMapping, startOuter, numberCellsOuter);
  }
  if (boundary[13] == boundaryType) {
    // right face
    numberCellsOuter[0] = 1;
    numberCellsOuter[1] = _moleculeContainer.getLocalNumberOfCells()[1];
    numberCellsOuter[2] = _moleculeContainer.getLocalNumberOfCells()[2];
    startOuter[0] = _moleculeContainer.getLocalNumberOfCells()[0] + 2 * _moleculeContainer.getLocalIndexOfFirstCell()[0] - 1;
    startOuter[1] = _moleculeContainer.getLocalIndexOfFirstCell()[1];
    startOuter[2] = _moleculeContainer.getLocalIndexOfFirstCell()[2];
    _moleculeContainer.iterateCells(myMapping, startOuter, numberCellsOuter);
  }
  if (boundary[10] == boundaryType) {
    // front face
    numberCellsOuter[0] = _moleculeContainer.getLocalNumberOfCells()[0];
    numberCellsOuter[1] = 1;
    numberCellsOuter[2] = _moleculeContainer.getLocalNumberOfCells()[2];
    startOuter[0] = _moleculeContainer.getLocalIndexOfFirstCell()[0];
    startOuter[1] = 0;
    startOuter[2] = _moleculeContainer.getLocalIndexOfFirstCell()[2];
    _moleculeContainer.iterateCells(myMapping, startOuter, numberCellsOuter);
  }
  if (boundary[15] == boundaryType) {
    // back face
    numberCellsOuter[0] = _moleculeContainer.getLocalNumberOfCells()[0];
    numberCellsOuter[1] = 1;
    numberCellsOuter[2] = _moleculeContainer.getLocalNumberOfCells()[2];
    startOuter[0] = _moleculeContainer.getLocalIndexOfFirstCell()[0];
    startOuter[1] = _moleculeContainer.getLocalNumberOfCells()[1] + 2 * _moleculeContainer.getLocalIndexOfFirstCell()[1] - 1;
    startOuter[2] = _moleculeContainer.getLocalIndexOfFirstCell()[2];
    _moleculeContainer.iterateCells(myMapping, startOuter, numberCellsOuter);
  }
#endif
}

template <class Mapping> void simplemd::BoundaryTreatment::applyMappingToOutermostNonBoundaryCells(Mapping& myMapping) const {
  tarch::la::Vector<MD_DIM, unsigned int> startOuter;
  tarch::la::Vector<MD_DIM, unsigned int> numberCellsOuter;
  bool emptySweep = true;

#if (MD_DIM == 1)
  startOuter[0] = 1;
  numberCellsOuter[0] = 1;
  _moleculeContainer.iterateCells(myMapping, startOuter, numberCellsOuter);
  startOuter[0] = _moleculeContainer.getLocalNumberOfCells()[0] + _moleculeContainer.getLocalIndexOfFirstCell()[0] - 1;
  _moleculeContainer.iterateCells(myMapping, startOuter, numberCellsOuter);
#endif

#if (MD_DIM == 2)
  // lower edge
  startOuter[0] = 1;
  startOuter[1] = 1;
  numberCellsOuter[0] = _moleculeContainer.getLocalNumberOfCells()[0];
  numberCellsOuter[1] = 1;
  _moleculeContainer.iterateCells(myMapping, startOuter, numberCellsOuter);
  // upper edge
  startOuter[0] = 1;
  startOuter[1] = _moleculeContainer.getLocalNumberOfCells()[1] + _moleculeContainer.getLocalIndexOfFirstCell()[1] - 1;
  _moleculeContainer.iterateCells(myMapping, startOuter, numberCellsOuter);
  // left edge
  startOuter[0] = 1;
  startOuter[1] = 2;
  numberCellsOuter[0] = 1;
  numberCellsOuter[1] = _moleculeContainer.getLocalNumberOfCells()[1] - 2;
  emptySweep = (numberCellsOuter[1] < 1);
  if (!emptySweep) {
    _moleculeContainer.iterateCells(myMapping, startOuter, numberCellsOuter);
  }
  // right edge
  startOuter[0] = _moleculeContainer.getLocalNumberOfCells()[0] + _moleculeContainer.getLocalIndexOfFirstCell()[0] - 1;
  startOuter[1] = 2;
  if (!emptySweep) {
    _moleculeContainer.iterateCells(myMapping, startOuter, numberCellsOuter);
  }
#endif

#if (MD_DIM == 3)
  // lower face
  numberCellsOuter[0] = _moleculeContainer.getLocalNumberOfCells()[0];
  numberCellsOuter[1] = _moleculeContainer.getLocalNumberOfCells()[1];
  numberCellsOuter[2] = 1;
  startOuter = tarch::la::Vector<MD_DIM, unsigned int>(1);
  _moleculeContainer.iterateCells(myMapping, startOuter, numberCellsOuter);
  // top face
  startOuter[0] = 1;
  startOuter[1] = 1;
  startOuter[2] = _moleculeContainer.getLocalNumberOfCells()[2] + _moleculeContainer.getLocalIndexOfFirstCell()[2] - 1;
  _moleculeContainer.iterateCells(myMapping, startOuter, numberCellsOuter);
  // left face
  numberCellsOuter[0] = 1;
  numberCellsOuter[1] = _moleculeContainer.getLocalNumberOfCells()[1];
  numberCellsOuter[2] = _moleculeContainer.getLocalNumberOfCells()[2] - 2;
  emptySweep = (numberCellsOuter[2] < 1);
  startOuter[0] = 1;
  startOuter[1] = 1;
  startOuter[2] = 2;
  if (!emptySweep) {
    _moleculeContainer.iterateCells(myMapping, startOuter, numberCellsOuter);
  }
  // right face
  startOuter[0] = _moleculeContainer.getLocalNumberOfCells()[0] + _moleculeContainer.getLocalIndexOfFirstCell()[0] - 1;
  startOuter[1] = 1;
  startOuter[2] = 2;
  if (!emptySweep) {
    _moleculeContainer.iterateCells(myMapping, startOuter, numberCellsOuter);
  }
  // front face
  numberCellsOuter[0] = _moleculeContainer.getLocalNumberOfCells()[0] - 2;
  numberCellsOuter[1] = 1;
  numberCellsOuter[2] = _moleculeContainer.getLocalNumberOfCells()[2] - 2;
  emptySweep = (numberCellsOuter[0] < 1) || (numberCellsOuter[2] < 1);
  startOuter[0] = 2;
  startOuter[1] = 1;
  startOuter[2] = 2;
  if (!emptySweep) {
    _moleculeContainer.iterateCells(myMapping, startOuter, numberCellsOuter);
  }
  // back face
  startOuter[0] = 2;
  startOuter[1] = _moleculeContainer.getLocalNumberOfCells()[1] + _moleculeContainer.getLocalIndexOfFirstCell()[1] - 1;
  startOuter[2] = 2;
  if (!emptySweep) {
    _moleculeContainer.iterateCells(myMapping, startOuter, numberCellsOuter);
  }
#endif
}

template <class Mapping> void simplemd::BoundaryTreatment::applyMappingToCommunicationIndependentCells(Mapping& myMapping) const {
  // starting point and range of cells, on which we can (for example) compute forces before the messages carrying boundary and process-leaving particles
  //  have arrived. Due to handling of iterateCellParis, we need to leave 1 inner cell on the "left" and 2 on the "right".
  //  Or, counting from ghost cells, we leave 2 cells on the "left" and 3 on the "right".
  tarch::la::Vector<MD_DIM, unsigned int> pairIterationStart(_moleculeContainer.getLocalIndexOfFirstCell());
  tarch::la::Vector<MD_DIM, unsigned int> pairIterationLength(_moleculeContainer.getLocalNumberOfCells());
  for (unsigned int d = 0; d < MD_DIM; d++) {
    pairIterationStart[d] += 1;
    pairIterationLength[d] -= 3;
  }

// apply force mapping
#if (MD_DEBUG == MD_YES)
  std::cout << "applying mapping on communication independent cells on pairIterationStart: " << pairIterationStart
            << "\n with pairIterationLength: " << pairIterationLength << std::endl;
#endif
  _moleculeContainer.iterateCellPairs(myMapping, pairIterationStart, pairIterationLength);
}

template <class Mapping> void simplemd::BoundaryTreatment::applyMappingToCommunicationDependentCells(Mapping& myMapping) const {
  // apply mapping to all cell-pairs, not processed by the above function

  tarch::la::Vector<MD_DIM, unsigned int> pairIterationStart(_moleculeContainer.getLocalIndexOfFirstCell());
  tarch::la::Vector<MD_DIM, unsigned int> pairIterationLength(_moleculeContainer.getLocalNumberOfCells());

#if (MD_DIM == 1)
  // left side
  pairIterationStart[0] = 0;
  pairIterationLength[0] = _moleculeContainer.getLocalIndexOfFirstCell()[0] + 1;
#if (MD_DEBUG == MD_YES)
  std::cout << "applying mapping on communication dependent cells with pairIterationStart: " << pairIterationStart
            << "\n with pairIterationLength: " << pairIterationLength << std::endl;
#endif
  _moleculeContainer.iterateCellPairs(myMapping, pairIterationStart, pairIterationLength);

  // right side
  // changing only necessary values:
  pairIterationStart[0] = _moleculeContainer.getLocalIndexOfFirstCell()[0] + _moleculeContainer.getLocalNumberOfCells()[0] - 2;
#if (MD_DEBUG == MD_YES)
  std::cout << "applying mapping on communication dependent cells with pairIterationStart: " << pairIterationStart
            << "\n with pairIterationLength: " << pairIterationLength << std::endl;
#endif
  _moleculeContainer.iterateCellPairs(myMapping, pairIterationStart, pairIterationLength);
#endif
#if (MD_DIM == 2)
  // whole lower edge
  pairIterationStart[0] = 0;
  pairIterationStart[1] = 0;
  pairIterationLength[0] = _moleculeContainer.getLocalIndexOfFirstCell()[0] + _moleculeContainer.getLocalNumberOfCells()[0];
  pairIterationLength[1] = _moleculeContainer.getLocalIndexOfFirstCell()[1] + 1;

#if (MD_DEBUG == MD_YES)
  std::cout << "applying mapping on communication dependent cells with pairIterationStart: " << pairIterationStart
            << "\n with pairIterationLength: " << pairIterationLength << std::endl;
#endif
  _moleculeContainer.iterateCellPairs(myMapping, pairIterationStart, pairIterationLength);

  // whole upper edge
  // changing only necessary values:
  pairIterationStart[1] = _moleculeContainer.getLocalIndexOfFirstCell()[1] + _moleculeContainer.getLocalNumberOfCells()[1] - 2;
#if (MD_DEBUG == MD_YES)
  std::cout << "applying mapping on communication dependent cells with pairIterationStart: " << pairIterationStart
            << "\n with pairIterationLength: " << pairIterationLength << std::endl;
#endif
  _moleculeContainer.iterateCellPairs(myMapping, pairIterationStart, pairIterationLength);

  // remaining of left edge
  pairIterationStart[0] = 0;
  pairIterationStart[1] = _moleculeContainer.getLocalIndexOfFirstCell()[1] + 1;
  pairIterationLength[0] = _moleculeContainer.getLocalIndexOfFirstCell()[0] + 1;
  pairIterationLength[1] = _moleculeContainer.getLocalNumberOfCells()[1] - 3;
#if (MD_DEBUG == MD_YES)
  std::cout << "applying mapping on communication dependent cells with pairIterationStart: " << pairIterationStart
            << "\n with pairIterationLength: " << pairIterationLength << std::endl;
#endif
  _moleculeContainer.iterateCellPairs(myMapping, pairIterationStart, pairIterationLength);

  // remaining of right edge
  // changing only necessary values:
  pairIterationStart[0] = _moleculeContainer.getLocalIndexOfFirstCell()[0] + _moleculeContainer.getLocalNumberOfCells()[0] - 2;
#if (MD_DEBUG == MD_YES)
  std::cout << "applying mapping on communication dependent cells with pairIterationStart: " << pairIterationStart
            << "\n with pairIterationLength: " << pairIterationLength << std::endl;
#endif
  _moleculeContainer.iterateCellPairs(myMapping, pairIterationStart, pairIterationLength);

#endif
#if (MD_DIM == 3)
  // whole lower face
  pairIterationStart[0] = 0;
  pairIterationStart[1] = 0;
  pairIterationStart[2] = 0;

  pairIterationLength[0] = _moleculeContainer.getLocalIndexOfFirstCell()[0] + _moleculeContainer.getLocalNumberOfCells()[0];
  pairIterationLength[1] = _moleculeContainer.getLocalIndexOfFirstCell()[1] + _moleculeContainer.getLocalNumberOfCells()[1];
  pairIterationLength[2] = _moleculeContainer.getLocalIndexOfFirstCell()[2] + 1;

#if (MD_DEBUG == MD_YES)
  std::cout << "applying mapping on communication dependent cells with pairIterationStart: " << pairIterationStart
            << "\n with pairIterationLength: " << pairIterationLength << std::endl;
#endif
  _moleculeContainer.iterateCellPairs(myMapping, pairIterationStart, pairIterationLength);

  // whole top face
  // changing only necessary values
  pairIterationStart[2] = _moleculeContainer.getLocalIndexOfFirstCell()[2] + _moleculeContainer.getLocalNumberOfCells()[2] - 2;

#if (MD_DEBUG == MD_YES)
  std::cout << "applying mapping on communication dependent cells with pairIterationStart: " << pairIterationStart
            << "\n with pairIterationLength: " << pairIterationLength << std::endl;
#endif
  _moleculeContainer.iterateCellPairs(myMapping, pairIterationStart, pairIterationLength);

  // remaining of left face
  pairIterationStart[0] = 0;
  pairIterationStart[1] = 0;
  pairIterationStart[2] = _moleculeContainer.getLocalIndexOfFirstCell()[2] + 1;

  pairIterationLength[0] = _moleculeContainer.getLocalIndexOfFirstCell()[0] + 1;
  pairIterationLength[1] = _moleculeContainer.getLocalIndexOfFirstCell()[1] + _moleculeContainer.getLocalNumberOfCells()[1];
  pairIterationLength[2] = _moleculeContainer.getLocalNumberOfCells()[2] - 3;

#if (MD_DEBUG == MD_YES)
  std::cout << "applying mapping on communication dependent cells with pairIterationStart: " << pairIterationStart
            << "\n with pairIterationLength: " << pairIterationLength << std::endl;
#endif
  _moleculeContainer.iterateCellPairs(myMapping, pairIterationStart, pairIterationLength);

  // remaining of right face
  pairIterationStart[0] = _moleculeContainer.getLocalIndexOfFirstCell()[0] + _moleculeContainer.getLocalNumberOfCells()[0] - 2;

#if (MD_DEBUG == MD_YES)
  std::cout << "applying mapping on communication dependent cells with pairIterationStart: " << pairIterationStart
            << "\n with pairIterationLength: " << pairIterationLength << std::endl;
#endif
  _moleculeContainer.iterateCellPairs(myMapping, pairIterationStart, pairIterationLength);

  // remaining of front face
  pairIterationStart[0] = _moleculeContainer.getLocalIndexOfFirstCell()[0] + 1;
  pairIterationStart[1] = 0;
  pairIterationStart[2] = _moleculeContainer.getLocalIndexOfFirstCell()[2] + 1;

  pairIterationLength[0] = _moleculeContainer.getLocalNumberOfCells()[0] - 3;
  pairIterationLength[1] = _moleculeContainer.getLocalIndexOfFirstCell()[1] + 1;
  pairIterationLength[2] = _moleculeContainer.getLocalNumberOfCells()[2] - 3;

#if (MD_DEBUG == MD_YES)
  std::cout << "applying mapping on communication dependent cells with pairIterationStart: " << pairIterationStart
            << "\n with pairIterationLength: " << pairIterationLength << std::endl;
#endif
  _moleculeContainer.iterateCellPairs(myMapping, pairIterationStart, pairIterationLength);

  // remaining of back face
  pairIterationStart[1] = _moleculeContainer.getLocalIndexOfFirstCell()[1] + _moleculeContainer.getLocalNumberOfCells()[1] - 2;

#if (MD_DEBUG == MD_YES)
  std::cout << "applying mapping on communication dependent cells with pairIterationStart: " << pairIterationStart
            << "\n with pairIterationLength: " << pairIterationLength << std::endl;
#endif
  _moleculeContainer.iterateCellPairs(myMapping, pairIterationStart, pairIterationLength);
#endif
}

#endif // _MOLECULARDYNAMICS_BOUNDARYTREATMENT_H_
