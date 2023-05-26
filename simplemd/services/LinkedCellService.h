// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_SERVICES_LINKEDCELLSERVICE_H_
#define _MOLECULARDYNAMICS_SERVICES_LINKEDCELLSERVICE_H_

#include "simplemd/LinkedCell.h"
#include "simplemd/Molecule.h"
#include "simplemd/molecule-mappings/UpdateLinkedCellListsMapping.h"
#include "simplemd/services/ParallelTopologyService.h"
#include "tarch/la/Vector.h"
#include <iostream>

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

  /** puts the molecule into the cell defined by the local index (vector)
   * coordinates localCellIndex */
  void addMoleculeToLinkedCell(Molecule& molecule, const tarch::la::Vector<MD_DIM, unsigned int>& localCellIndex);
  /** puts the molecule into the cell defined by the local index (scalar)
   * coordinates localCellIndex */
  void addMoleculeToLinkedCell(Molecule& molecule, const unsigned int& localCellIndex);

  void deleteMoleculeFromLinkedCell(Molecule& molecule, const tarch::la::Vector<MD_DIM, unsigned int>& localCellIndex);

  /** returns the linked cell at the respective coordinates */
  LinkedCell& getLinkedCell(const tarch::la::Vector<MD_DIM, unsigned int>& localCellIndex);

  /** iterates over all cells in the range defined by the lower left front
   * corner cell lowerLeftFrontCell and the size of the domain cellRange.
   * cellRange defines a number of cells in each spatial direction that the
   * class A shall be applied to. lowerLeftFrontCell needs to be given in local
   * coordinates.
   */
  template <class A>
  void iterateCells(A& a, const tarch::la::Vector<MD_DIM, unsigned int>& lowerLeftFrontCell, const tarch::la::Vector<MD_DIM, unsigned int>& cellRange,
                    const bool& useOpenMP);

  /** iterates over all cells in the inner part (i.e. does not consider the
   * ghost layer) */
  template <class A> void iterateCells(A& a, const bool& useOpenMP);

  /** iterates over all cell pairs for the cells in the inner region of each
   * local process */
  template <class A> void iterateCellPairs(A& a, const bool& useOpenMP) const;

  /** iterates over all cell pairs cell1 and cell2 with cell1 in the range
   * described by lowerLeftFrontCell and cellRange; cell2 does not need to lie
   * within the range (example: iterate only over lowerLeftFrontCell=(1,1,1) and
   * cellRange=(1,1,1). Then, we will consider amongst others the pair
   * (0,0,0),(1,1,1)).
   */
  template <class A>
  void iterateCellPairs(A& a, const tarch::la::Vector<MD_DIM, unsigned int>& lowerLeftFrontCell, const tarch::la::Vector<MD_DIM, unsigned int>& cellRange,
                        const bool& useOpenMP) const;

  /** iterates over all cell triplets for the cells in the inner region of each
   * local process */
  template <class A> void iterateCellTriplets(A& a, const bool& useOpenMP) const;

  /** iterates over all cell triplets cell1, cell2 and cell3 with cell1 in the range
   * described by lowerLeftFrontCell and cellRange; cell2 and cell3 do not need to lie
   * within the range (example: iterate only over lowerLeftFrontCell=(1,1,1) and
   * cellRange=(1,1,1). Then, we will consider amongst others the triplet
   * (0,0,0),(0,0,1),(1,1,1)).
   */
  template <class A>
  void iterateCellTriplets(A& a, const tarch::la::Vector<MD_DIM, unsigned int>& lowerLeftFrontCell, const tarch::la::Vector<MD_DIM, unsigned int>& cellRange,
                        const bool& useOpenMP) const;

  /** returns the index of the first (non-ghost) cell */
  const tarch::la::Vector<MD_DIM, unsigned int>& getLocalIndexOfFirstCell() const;

  /** returns the number of (non-ghost) cells */
  const tarch::la::Vector<MD_DIM, unsigned int>& getLocalNumberOfCells() const;

  /** returns the mesh width */
  const tarch::la::Vector<MD_DIM, double>& getMeshWidth() const;

  /** returns the local domain offset (for the domain of this process) */
  const tarch::la::Vector<MD_DIM, double>& getLocalDomainOffset() const;

  /** returns the local domain size (for the domain of this process) */
  const tarch::la::Vector<MD_DIM, double>& getLocalDomainSize() const;

  /** returns true if the local cell index cellIndex describes a linked cell
   * within the ghost layer */
  bool isGhostCell(const unsigned int& cellIndex) const;

  /** returns the local cell index vector for the local cell index cellIndex */
  tarch::la::Vector<MD_DIM, unsigned int> getLocalCellIndexVector(const unsigned int cellIndex) const;

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

template <class A>
void simplemd::services::LinkedCellService::iterateCells(A& a, const tarch::la::Vector<MD_DIM, unsigned int>& lowerLeftFrontCell,
                                                         const tarch::la::Vector<MD_DIM, unsigned int>& cellRange, const bool& useOpenMP) {
  unsigned int index = 0;
#if (MD_ERROR == MD_YES)
  for (unsigned int d = 0; d < MD_DIM; d++) {
    if (cellRange[d] == 0) {
      std::cout << "ERROR simplemd::services::LinkedCellService::iterateCells: "
                   "cellRange("
                << d << ")==0!" << std::endl;
      exit(EXIT_FAILURE);
    }
    if (lowerLeftFrontCell[d] + cellRange[d] > 2 * _indexOffset[d] + _numberOfCells[d]) {
      std::cout << "ERROR simplemd::services::LinkedCellService::iterateCells(): "
                   "defined Range does not fit into local sub-domain!"
                << std::endl;
      exit(EXIT_FAILURE);
    }
  }
#endif

  // start iteration();
  a.beginCellIteration();

#if (MD_OPENMP == MD_YES)
  if (useOpenMP) {
    const tarch::la::Vector<MD_DIM, unsigned int> size(simplemd::services::LinkedCellService::getInstance().getLocalNumberOfCells() +
                                                       2 * simplemd::services::LinkedCellService::getInstance().getLocalIndexOfFirstCell());
    const int length = cellRange(0)
#if (MD_DIM > 1)
                       * cellRange(1)
#endif
#if (MD_DIM > 2)
                       * cellRange(2)
#endif
        ;
// loop over domain, but with a single loop
#pragma omp parallel for
    for (int i = 0; i < length; i++) {
// compute index of the current cell
#if (MD_DIM > 1)
      int helpIndex1 = i;
      int helpIndex2 = 0;
#endif
      unsigned int index = 0;

#if (MD_DIM > 2)
      // determine plane within traversed block
      helpIndex2 = helpIndex1 / (cellRange(0) * cellRange(1));
      // save rest of index in helpIndex1
      helpIndex1 = helpIndex1 - helpIndex2 * (cellRange(0) * cellRange(1));
      // compute contribution to index
      index += (lowerLeftFrontCell(2) + helpIndex2) * size(0) * size(1);
#endif
#if (MD_DIM > 1)
      // determine plane within traversed block
      helpIndex2 = helpIndex1 / cellRange(0);
      // save rest of index in helpIndex1
      helpIndex1 = helpIndex1 - helpIndex2 * cellRange(0);
      // compute contribution to index
      index += (lowerLeftFrontCell(1) + helpIndex2) * size(0);
      // compute contribution for last dimension
      index += (lowerLeftFrontCell(0) + helpIndex1);
#else
      index = lowerLeftFrontCell(0) + i;
#endif
#if (MD_DEBUG == MD_YES)
      std::cout << "Handle cell " << index << std::endl;
#endif

      // handle cell
      a.handleCell(_cells[index], index);
    }
  } else {
#endif

    tarch::la::Vector<MD_DIM, unsigned int> coords(0);
// loop over domain
#if (MD_DIM > 2)
    for (coords[2] = lowerLeftFrontCell[2]; coords[2] < lowerLeftFrontCell[2] + cellRange[2]; coords[2]++) {
#endif
#if (MD_DIM > 1)
      for (coords[1] = lowerLeftFrontCell[1]; coords[1] < lowerLeftFrontCell[1] + cellRange[1]; coords[1]++) {
#endif
        for (coords[0] = lowerLeftFrontCell[0]; coords[0] < lowerLeftFrontCell[0] + cellRange[0]; coords[0]++) {
#if (MD_DEBUG == MD_YES)
          std::cout << "Handle cell " << coords << std::endl;
#endif
          index = getLocalIndexFromLocalVector(coords);
          a.handleCell(_cells[index], index);
        }
#if (MD_DIM > 1)
      }
#endif
#if (MD_DIM > 2)
    }
#endif

#if (MD_OPENMP == MD_YES)
  }
#endif

  // end iteration();
  a.endCellIteration();
}

template <class A> void simplemd::services::LinkedCellService::iterateCells(A& a, const bool& useOpenMP) {
  iterateCells(a, _indexOffset, _numberOfCells, useOpenMP);
}

template <class A>
void simplemd::services::LinkedCellService::iterateCellPairs(A& a, const tarch::la::Vector<MD_DIM, unsigned int>& lowerLeftFrontCell,
                                                             const tarch::la::Vector<MD_DIM, unsigned int>& cellRange, const bool& useOpenMP) const {
  tarch::la::Vector<MD_LINKED_CELL_NEIGHBOURS / 2, unsigned int> neighbourOffset;
  tarch::la::Vector<MD_LINKED_CELL_NEIGHBOURS / 2, unsigned int> indexOffset;
#if (MD_DIM == 1)
  neighbourOffset[0] = 1;
  indexOffset[0] = 0;
#elif (MD_DIM == 2)
  indexOffset[0] = 0;
  neighbourOffset[0] = 1;
  indexOffset[1] = 0;
  neighbourOffset[1] = _numberOfCells[0] + 2;
  indexOffset[2] = 0;
  neighbourOffset[2] = _numberOfCells[0] + 3;
  indexOffset[3] = 1;
  neighbourOffset[3] = _numberOfCells[0] + 2;
#elif (MD_DIM == 3)
  indexOffset[0] = 0;
  neighbourOffset[0] = 1;
  indexOffset[1] = 0;
  neighbourOffset[1] = _numberOfCells[0] + 2;
  indexOffset[2] = 0;
  neighbourOffset[2] = _numberOfCells[0] + 3;
  indexOffset[3] = 0;
  neighbourOffset[3] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2);
  indexOffset[4] = 0;
  neighbourOffset[4] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2) + 1;
  indexOffset[5] = 0;
  neighbourOffset[5] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2) + (_numberOfCells[0] + 2);
  indexOffset[6] = 0;
  neighbourOffset[6] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2) + (_numberOfCells[0] + 2) + 1;

  indexOffset[7] = 1;
  neighbourOffset[7] = _numberOfCells[0] + 2;
  indexOffset[8] = 1;
  neighbourOffset[8] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2);
  indexOffset[9] = 1;
  neighbourOffset[9] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2) + (_numberOfCells[0] + 2);

  indexOffset[10] = _numberOfCells[0] + 2;
  neighbourOffset[10] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2);
  indexOffset[11] = _numberOfCells[0] + 2;
  neighbourOffset[11] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2) + 1;

  indexOffset[12] = (_numberOfCells[0] + 2) + 1;
  neighbourOffset[12] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2);
#endif

#if (MD_ERROR == MD_YES)
  for (unsigned int d = 0; d < MD_DIM; d++) {
    if (cellRange[d] > 2 * _indexOffset[d] + _numberOfCells[d] - 1) {
      std::cout << "ERROR simplemd::services::LinkedCellService::iterateCellPairs(): "
                   "defined Range does not fit into local sub-domain!"
                << std::endl;
      exit(EXIT_FAILURE);
    }
  }
#endif

  // start iteration();
  a.beginCellIteration();

#if (MD_OPENMP == MD_YES)
  if (useOpenMP) {
    const tarch::la::Vector<MD_DIM, unsigned int> size(simplemd::services::LinkedCellService::getInstance().getLocalNumberOfCells() +
                                                       2 * simplemd::services::LinkedCellService::getInstance().getLocalIndexOfFirstCell());

// iterate over the domain in a red-black manner
#if (MD_DIM > 2)
    for (unsigned int z = 0; z < 2; z++) {
#endif
#if (MD_DIM > 1)
      for (unsigned int y = 0; y < 2; y++) {
#endif
        for (unsigned int x = 0; x < 2; x++) {
          // determine range/ length of blocks for red-black traversal.
          // For odd block sizes, we need to do some more work in the
          // x/y/z==0-traversals. The second x/y/z==1-traversals are reduced by
          // the normal integer-rounding in this case.
          const tarch::la::Vector<MD_DIM, unsigned int> lengthVector((cellRange[0] + (cellRange[0] % 2) * (x == 0)) / 2
#if (MD_DIM > 1)
                                                                     ,
                                                                     (cellRange[1] + (cellRange[1] % 2) * (y == 0)) / 2
#endif
#if (MD_DIM > 2)
                                                                     ,
                                                                     (cellRange[2] + (cellRange[2] % 2) * (z == 0)) / 2
#endif
          );
          const int length = lengthVector[0]
#if (MD_DIM > 1)
                             * lengthVector[1]
#endif
#if (MD_DIM > 2)
                             * lengthVector[2]
#endif
              ;

// parallelise loop for all cells that are to be traversed in this way
#pragma omp parallel for
          for (int j = 0; j < length; j++) {
            // compute index of the current cell
            unsigned int index = 0;
#if (MD_DIM > 1)
            int helpIndex1 = j;
            int helpIndex2 = 0;
#endif
            unsigned int coordsCell1Buffer;
            unsigned int coordsCell2Buffer;

#if (MD_DIM > 2)
            // determine plane within traversed block
            helpIndex2 = helpIndex1 / (lengthVector[0] * lengthVector[1]);
            // save rest of index in helpIndex1
            helpIndex1 = helpIndex1 - helpIndex2 * (lengthVector[0] * lengthVector[1]);
            // compute contribution to index
            index += (lowerLeftFrontCell[2] + 2 * helpIndex2 + z) * size[0] * size[1];
#endif
#if (MD_DIM > 1)
            // determine plane within traversed block
            helpIndex2 = helpIndex1 / lengthVector[0];
            // save rest of index in helpIndex1
            helpIndex1 = helpIndex1 - helpIndex2 * lengthVector[0];
            // compute contribution to index
            index += (lowerLeftFrontCell[1] + 2 * helpIndex2 + y) * size[0];
            // compute contribution for last dimension
            index += (lowerLeftFrontCell[0] + 2 * helpIndex1 + x);
#else
        index = lowerLeftFrontCell[0] + 2 * j + x;
#endif
#if (MD_DEBUG == MD_YES)
            std::cout << "Handle cell " << index << std::endl;
#endif

            a.handleCell(_cells[index], index);
            // handle pairs (lower,left,back-oriented cells)
            for (unsigned int i = 0; i < MD_LINKED_CELL_NEIGHBOURS / 2; i++) {
#if (MD_DEBUG == MD_YES)
              std::cout << "iterateCellPairs: Pair index " << index + indexOffset[i] << "," << index + neighbourOffset[i] << std::endl;
#endif
              coordsCell1Buffer = index + indexOffset[i];
              coordsCell2Buffer = index + neighbourOffset[i];
              a.handleCellPair(_cells[coordsCell1Buffer], _cells[coordsCell2Buffer], coordsCell1Buffer, coordsCell2Buffer);
            }
          } // j
        }   // x
#if (MD_DIM > 1)
      } // y
#endif
#if (MD_DIM > 2)
    } // z
#endif
    // now: no open mp
  } else {
#endif

    tarch::la::Vector<MD_DIM, unsigned int> coords(0);
    unsigned int coordsCell1Buffer(0);
    unsigned int coordsCell2Buffer(0);
    unsigned int index;

// loop over domain
#if (MD_DIM > 2)
    for (coords[2] = lowerLeftFrontCell[2]; coords[2] < lowerLeftFrontCell[2] + cellRange[2]; coords[2]++) {
#endif
#if (MD_DIM > 1)
      for (coords[1] = lowerLeftFrontCell[1]; coords[1] < lowerLeftFrontCell[1] + cellRange[1]; coords[1]++) {
#endif
        for (coords[0] = lowerLeftFrontCell[0]; coords[0] < lowerLeftFrontCell[0] + cellRange[0]; coords[0]++) {
          // handle cell itself
          index = getLocalIndexFromLocalVector(coords);
#if (MD_DEBUG == MD_YES)
          std::cout << "iterateCellPairs: Single index " << index << std::endl;
#endif

          a.handleCell(_cells[index], index);
          // handle pairs (lower,left,back-oriented cells)
          for (unsigned int i = 0; i < MD_LINKED_CELL_NEIGHBOURS / 2; i++) {
#if (MD_DEBUG == MD_YES)
            std::cout << "iterateCellPairs: Pair index " << index + indexOffset[i] << "," << index + neighbourOffset[i] << std::endl;
#endif

            coordsCell1Buffer = index + indexOffset[i];
            coordsCell2Buffer = index + neighbourOffset[i];
            a.handleCellPair(_cells[coordsCell1Buffer], _cells[coordsCell2Buffer], coordsCell1Buffer, coordsCell2Buffer);
          }
        } // coords(0)
#if (MD_DIM > 1)
      }
#endif
#if (MD_DIM > 2)
    }
#endif

#if (MD_OPENMP == MD_YES)
  }
#endif

  // end iteration();
  a.endCellIteration();
}

template <class A> void simplemd::services::LinkedCellService::iterateCellPairs(A& a, const bool& useOpenMP) const {
  const tarch::la::Vector<MD_DIM, unsigned int> pairIterationStart(0);
  const tarch::la::Vector<MD_DIM, unsigned int> pairIterationLength(getLocalNumberOfCells() + getLocalIndexOfFirstCell());
  iterateCellPairs(a, pairIterationStart, pairIterationLength, useOpenMP);
}

template <class A>
void simplemd::services::LinkedCellService::iterateCellTriplets(A& a, const tarch::la::Vector<MD_DIM, unsigned int>& lowerLeftFrontCell,
                                                             const tarch::la::Vector<MD_DIM, unsigned int>& cellRange, const bool& useOpenMP) const {
  tarch::la::Vector<MD_LINKED_CELL_NEIGHBOURS / 2, unsigned int> neighbourOffset;
  tarch::la::Vector<MD_LINKED_CELL_NEIGHBOURS / 2, unsigned int> indexOffset;
#if (MD_DIM == 2)
  tarch::la::Vector<4, unsigned int> tripletOffset;
#elif (MD_DIM == 3)
  tarch::la::Vector<44, unsigned int> tripletIndexOffset;
  tarch::la::Vector<44, unsigned int> tripletNeighbourOffset1;
  tarch::la::Vector<44, unsigned int> tripletNeighbourOffset2;
#endif

#if (MD_DIM == 1)
  neighbourOffset[0] = 1;
  indexOffset[0] = 0;
#elif (MD_DIM == 2)
  indexOffset[0] = 0;
  neighbourOffset[0] = 1;
  tripletOffset[0] = _numberOfCells[0] + 2;
  indexOffset[1] = 0;
  neighbourOffset[1] = _numberOfCells[0] + 2;
  tripletOffset[1] = _numberOfCells[0] + 3;
  indexOffset[2] = 0;
  neighbourOffset[2] = _numberOfCells[0] + 3;
  tripletOffset[2] = 1;
  indexOffset[3] = 1;
  neighbourOffset[3] = _numberOfCells[0] + 2;
  tripletOffset[3] = _numberOfCells[0] + 3;
#elif (MD_DIM == 3)
  indexOffset[0] = 0;
  neighbourOffset[0] = 1;
  indexOffset[1] = 0;
  neighbourOffset[1] = _numberOfCells[0] + 2;
  indexOffset[2] = 0;
  neighbourOffset[2] = _numberOfCells[0] + 3;
  indexOffset[3] = 0;
  neighbourOffset[3] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2);
  indexOffset[4] = 0;
  neighbourOffset[4] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2) + 1;
  indexOffset[5] = 0;
  neighbourOffset[5] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2) + (_numberOfCells[0] + 2);
  indexOffset[6] = 0;
  neighbourOffset[6] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2) + (_numberOfCells[0] + 2) + 1;

  indexOffset[7] = 1;
  neighbourOffset[7] = _numberOfCells[0] + 2;
  indexOffset[8] = 1;
  neighbourOffset[8] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2);
  indexOffset[9] = 1;
  neighbourOffset[9] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2) + (_numberOfCells[0] + 2);

  indexOffset[10] = _numberOfCells[0] + 2;
  neighbourOffset[10] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2);
  indexOffset[11] = _numberOfCells[0] + 2;
  neighbourOffset[11] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2) + 1;

  indexOffset[12] = (_numberOfCells[0] + 2) + 1;
  neighbourOffset[12] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2);

  tripletIndexOffset[0] = 0;
  tripletNeighbourOffset1[0] = 1;
  tripletNeighbourOffset2[0] = _numberOfCells[0] + 2;
  tripletIndexOffset[1] = 0;
  tripletNeighbourOffset1[1] = 1;
  tripletNeighbourOffset2[1] = _numberOfCells[0] + 3;
  tripletIndexOffset[2] = 0;
  tripletNeighbourOffset1[2] = _numberOfCells[0] + 2;
  tripletNeighbourOffset2[2] = _numberOfCells[0] + 3;
  tripletIndexOffset[3] = 1;
  tripletNeighbourOffset1[3] = _numberOfCells[0] + 2;
  tripletNeighbourOffset2[3] = _numberOfCells[0] + 3;

  tripletIndexOffset[4] = 0;
  tripletNeighbourOffset1[4] = 1;
  tripletNeighbourOffset2[4] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2);
  tripletIndexOffset[5] = 0;
  tripletNeighbourOffset1[5] = 1;
  tripletNeighbourOffset2[5] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2) + 1;
  tripletIndexOffset[6] = 0;
  tripletNeighbourOffset1[6] = 1;
  tripletNeighbourOffset2[6] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2) + (_numberOfCells[0] + 2);
  tripletIndexOffset[7] = 0;
  tripletNeighbourOffset1[7] = 1;
  tripletNeighbourOffset2[7] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2) + (_numberOfCells[0] + 2) + 1;
  tripletIndexOffset[8] = 0;
  tripletNeighbourOffset1[8] = _numberOfCells[0] + 2;
  tripletNeighbourOffset2[8] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2);
  tripletIndexOffset[9] = 0;
  tripletNeighbourOffset1[9] = _numberOfCells[0] + 2;
  tripletNeighbourOffset2[9] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2) + 1;
  tripletIndexOffset[10] = 0;
  tripletNeighbourOffset1[10] = _numberOfCells[0] + 2;
  tripletNeighbourOffset2[10] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2) + (_numberOfCells[0] + 2);
  tripletIndexOffset[11] = 0;
  tripletNeighbourOffset1[11] = _numberOfCells[0] + 2;
  tripletNeighbourOffset2[11] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2) + (_numberOfCells[0] + 2) + 1;
  tripletIndexOffset[12] = 0;
  tripletNeighbourOffset1[12] = _numberOfCells[0] + 3;
  tripletNeighbourOffset2[12] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2);
  tripletIndexOffset[13] = 0;
  tripletNeighbourOffset1[13] = _numberOfCells[0] + 3;
  tripletNeighbourOffset2[13] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2) + 1;
  tripletIndexOffset[14] = 0;
  tripletNeighbourOffset1[14] = _numberOfCells[0] + 3;
  tripletNeighbourOffset2[14] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2) + (_numberOfCells[0] + 2);
  tripletIndexOffset[15] = 0;
  tripletNeighbourOffset1[15] = _numberOfCells[0] + 3;
  tripletNeighbourOffset2[15] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2) + (_numberOfCells[0] + 2) + 1;
  tripletIndexOffset[16] = 1;
  tripletNeighbourOffset1[16] = _numberOfCells[0] + 2;
  tripletNeighbourOffset2[16] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2);
  tripletIndexOffset[17] = 1;
  tripletNeighbourOffset1[17] = _numberOfCells[0] + 2;
  tripletNeighbourOffset2[17] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2) + 1;
  tripletIndexOffset[18] = 1;
  tripletNeighbourOffset1[18] = _numberOfCells[0] + 2;
  tripletNeighbourOffset2[18] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2) + (_numberOfCells[0] + 2);
  tripletIndexOffset[19] = 1;
  tripletNeighbourOffset1[19] = _numberOfCells[0] + 2;
  tripletNeighbourOffset2[19] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2) + (_numberOfCells[0] + 2) + 1;
  tripletIndexOffset[20] = 1;
  tripletNeighbourOffset1[20] = _numberOfCells[0] + 3;
  tripletNeighbourOffset2[20] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2);
  tripletIndexOffset[21] = 1;
  tripletNeighbourOffset1[21] = _numberOfCells[0] + 3;
  tripletNeighbourOffset2[21] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2) + (_numberOfCells[0] + 2);
  tripletIndexOffset[22] = _numberOfCells[0] + 2;
  tripletNeighbourOffset1[22] = _numberOfCells[0] + 3;
  tripletNeighbourOffset2[22] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2);
  tripletIndexOffset[23] = _numberOfCells[0] + 2;
  tripletNeighbourOffset1[23] = _numberOfCells[0] + 3;
  tripletNeighbourOffset2[23] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2) + 1;

  tripletIndexOffset[24] = 0;
  tripletNeighbourOffset1[24] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2);
  tripletNeighbourOffset2[24] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2) + 1;
  tripletIndexOffset[25] = 1;
  tripletNeighbourOffset1[25] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2);
  tripletNeighbourOffset2[25] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2) + 1;
  tripletIndexOffset[26] = _numberOfCells[0] + 2;
  tripletNeighbourOffset1[26] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2);
  tripletNeighbourOffset2[26] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2) + 1;
  tripletIndexOffset[27] = _numberOfCells[0] + 3;
  tripletNeighbourOffset1[27] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2);
  tripletNeighbourOffset2[27] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2) + 1;
  tripletIndexOffset[28] = 0;
  tripletNeighbourOffset1[28] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2);
  tripletNeighbourOffset2[28] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2) + (_numberOfCells[0] + 2);
  tripletIndexOffset[29] = 1;
  tripletNeighbourOffset1[29] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2);
  tripletNeighbourOffset2[29] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2) + (_numberOfCells[0] + 2);
  tripletIndexOffset[30] = _numberOfCells[0] + 2;
  tripletNeighbourOffset1[30] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2);
  tripletNeighbourOffset2[30] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2) + (_numberOfCells[0] + 2);
  tripletIndexOffset[31] = _numberOfCells[0] + 3;
  tripletNeighbourOffset1[31] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2);
  tripletNeighbourOffset2[31] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2) + (_numberOfCells[0] + 2);
  tripletIndexOffset[32] = 0;
  tripletNeighbourOffset1[32] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2);
  tripletNeighbourOffset2[32] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2) + (_numberOfCells[0] + 2) + 1;
  tripletIndexOffset[33] = 1;
  tripletNeighbourOffset1[33] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2);
  tripletNeighbourOffset2[33] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2) + (_numberOfCells[0] + 2) + 1;
  tripletIndexOffset[34] = _numberOfCells[0] + 2;
  tripletNeighbourOffset1[34] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2);
  tripletNeighbourOffset2[34] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2) + (_numberOfCells[0] + 2) + 1;
  tripletIndexOffset[35] = _numberOfCells[0] + 3;
  tripletNeighbourOffset1[35] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2);
  tripletNeighbourOffset2[35] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2) + (_numberOfCells[0] + 2) + 1;
  tripletIndexOffset[36] = 0;
  tripletNeighbourOffset1[36] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2) + 1;
  tripletNeighbourOffset2[36] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2) + (_numberOfCells[0] + 2);
  tripletIndexOffset[37] = 1;
  tripletNeighbourOffset1[37] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2) + 1;
  tripletNeighbourOffset2[37] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2) + (_numberOfCells[0] + 2);
  tripletIndexOffset[38] = _numberOfCells[0] + 2;
  tripletNeighbourOffset1[38] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2) + 1;
  tripletNeighbourOffset2[38] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2) + (_numberOfCells[0] + 2);
  tripletIndexOffset[39] = _numberOfCells[0] + 3;
  tripletNeighbourOffset1[39] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2) + 1;
  tripletNeighbourOffset2[39] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2) + (_numberOfCells[0] + 2);
  tripletIndexOffset[40] = 0;
  tripletNeighbourOffset1[40] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2) + 1;
  tripletNeighbourOffset2[40] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2) + (_numberOfCells[0] + 2) + 1;
  tripletIndexOffset[41] = _numberOfCells[0] + 2;
  tripletNeighbourOffset1[41] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2) + 1;
  tripletNeighbourOffset2[41] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2) + (_numberOfCells[0] + 2) + 1;
  tripletIndexOffset[42] = 0;
  tripletNeighbourOffset1[42] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2) + (_numberOfCells[0] + 2);
  tripletNeighbourOffset2[42] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2) + (_numberOfCells[0] + 2) + 1;
  tripletIndexOffset[43] = 1;
  tripletNeighbourOffset1[43] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2) + (_numberOfCells[0] + 2);
  tripletNeighbourOffset2[43] = (_numberOfCells[0] + 2) * (_numberOfCells[1] + 2) + (_numberOfCells[0] + 2) + 1;
#endif

#if (MD_ERROR == MD_YES)
  for (unsigned int d = 0; d < MD_DIM; d++) {
    if (cellRange[d] > 2 * _indexOffset[d] + _numberOfCells[d] - 1) {
      std::cout << "ERROR simplemd::services::LinkedCellService::iterateCellPairs(): "
                   "defined Range does not fit into local sub-domain!"
                << std::endl;
      exit(EXIT_FAILURE);
    }
  }
#endif

  // start iteration();
  a.beginCellIteration();

#if (MD_OPENMP == MD_YES) // TODO
  if (useOpenMP) {
    const tarch::la::Vector<MD_DIM, unsigned int> size(simplemd::services::LinkedCellService::getInstance().getLocalNumberOfCells() +
                                                       2 * simplemd::services::LinkedCellService::getInstance().getLocalIndexOfFirstCell());

// iterate over the domain in a red-black manner
#if (MD_DIM > 2)
    for (unsigned int z = 0; z < 2; z++) {
#endif
#if (MD_DIM > 1)
      for (unsigned int y = 0; y < 2; y++) {
#endif
        for (unsigned int x = 0; x < 2; x++) {
          // determine range/ length of blocks for red-black traversal.
          // For odd block sizes, we need to do some more work in the
          // x/y/z==0-traversals. The second x/y/z==1-traversals are reduced by
          // the normal integer-rounding in this case.
          const tarch::la::Vector<MD_DIM, unsigned int> lengthVector((cellRange[0] + (cellRange[0] % 2) * (x == 0)) / 2
#if (MD_DIM > 1)
                                                                     ,
                                                                     (cellRange[1] + (cellRange[1] % 2) * (y == 0)) / 2
#endif
#if (MD_DIM > 2)
                                                                     ,
                                                                     (cellRange[2] + (cellRange[2] % 2) * (z == 0)) / 2
#endif
          );
          const int length = lengthVector[0]
#if (MD_DIM > 1)
                             * lengthVector[1]
#endif
#if (MD_DIM > 2)
                             * lengthVector[2]
#endif
              ;

// parallelise loop for all cells that are to be traversed in this way
#pragma omp parallel for
          for (int j = 0; j < length; j++) {
            // compute index of the current cell
            unsigned int index = 0;
#if (MD_DIM > 1)
            int helpIndex1 = j;
            int helpIndex2 = 0;
#endif
            unsigned int coordsCell1Buffer;
            unsigned int coordsCell2Buffer;

#if (MD_DIM > 2)
            // determine plane within traversed block
            helpIndex2 = helpIndex1 / (lengthVector[0] * lengthVector[1]);
            // save rest of index in helpIndex1
            helpIndex1 = helpIndex1 - helpIndex2 * (lengthVector[0] * lengthVector[1]);
            // compute contribution to index
            index += (lowerLeftFrontCell[2] + 2 * helpIndex2 + z) * size[0] * size[1];
#endif
#if (MD_DIM > 1)
            // determine plane within traversed block
            helpIndex2 = helpIndex1 / lengthVector[0];
            // save rest of index in helpIndex1
            helpIndex1 = helpIndex1 - helpIndex2 * lengthVector[0];
            // compute contribution to index
            index += (lowerLeftFrontCell[1] + 2 * helpIndex2 + y) * size[0];
            // compute contribution for last dimension
            index += (lowerLeftFrontCell[0] + 2 * helpIndex1 + x);
#else
        index = lowerLeftFrontCell[0] + 2 * j + x;
#endif
#if (MD_DEBUG == MD_YES)
            std::cout << "Handle cell " << index << std::endl;
#endif

            a.handleCell(_cells[index], index);
            // handle pairs (lower,left,back-oriented cells)
            for (unsigned int i = 0; i < MD_LINKED_CELL_NEIGHBOURS / 2; i++) {
#if (MD_DEBUG == MD_YES)
              std::cout << "iterateCellPairs: Pair index " << index + indexOffset[i] << "," << index + neighbourOffset[i] << std::endl;
#endif
              coordsCell1Buffer = index + indexOffset[i];
              coordsCell2Buffer = index + neighbourOffset[i];
              a.handleCellPair(_cells[coordsCell1Buffer], _cells[coordsCell2Buffer], coordsCell1Buffer, coordsCell2Buffer);
            }
          } // j
        }   // x
#if (MD_DIM > 1)
      } // y
#endif
#if (MD_DIM > 2)
    } // z
#endif
    // now: no open mp
  } else {
#endif

    tarch::la::Vector<MD_DIM, unsigned int> coords(0);
    unsigned int coordsCell1Buffer(0);
    unsigned int coordsCell2Buffer(0);
    unsigned int coordsCell3Buffer(0);
    unsigned int index;

// loop over domain
#if (MD_DIM > 2)
    for (coords[2] = lowerLeftFrontCell[2]; coords[2] < lowerLeftFrontCell[2] + cellRange[2]; coords[2]++) {
#endif
#if (MD_DIM > 1)
      for (coords[1] = lowerLeftFrontCell[1]; coords[1] < lowerLeftFrontCell[1] + cellRange[1]; coords[1]++) {
#endif
        for (coords[0] = lowerLeftFrontCell[0]; coords[0] < lowerLeftFrontCell[0] + cellRange[0]; coords[0]++) {
          // handle cell itself
          index = getLocalIndexFromLocalVector(coords);
#if (MD_DEBUG == MD_YES)
          std::cout << "iterateCellTriplets: Single index " << index << std::endl;
#endif

          a.handleCell(_cells[index], index);
          // handle pairs (lower,left,back-oriented cells)
          for (unsigned int i = 0; i < MD_LINKED_CELL_NEIGHBOURS / 2; i++) {
#if (MD_DEBUG == MD_YES)
            std::cout << "iterateCellTriplets: Pair index " << index + indexOffset[i] << "," << index + neighbourOffset[i] << std::endl;
#endif

            coordsCell1Buffer = index + indexOffset[i];
            coordsCell2Buffer = index + neighbourOffset[i];
            a.handleCellPair(_cells[coordsCell1Buffer], _cells[coordsCell2Buffer], coordsCell1Buffer, coordsCell2Buffer);

#if (MD_DIM == 2)
#if (MD_DEBUG == MD_YES)
            std::cout << "iterateCellTriplets: Triplet index " << index + indexOffset[i] << ","
                      << index + neighbourOffset[i] << "," << index + tripletOffset[i] << std::endl;
#endif 

            coordsCell3Buffer = index + tripletOffset[i];
            a.handleCellTriplet(_cells[coordsCell1Buffer], _cells[coordsCell2Buffer], _cells[coordsCell3Buffer],
                                coordsCell1Buffer, coordsCell2Buffer, coordsCell3Buffer);
          }
#elif (MD_DIM == 3)
          }

          for (unsigned int i = 0; i < 44; i++) {
#if (MD_DEBUG == MD_YES)
            std::cout << "iterateCellTriplets: Triplet index " << index + tripletIndexOffset[i] << ","
                      << index + tripletNeighbourOffset1[i] << "," << index + tripletNeighbourOffset2[i] << std::endl;
#endif

            coordsCell1Buffer = index + tripletIndexOffset[i];
            coordsCell2Buffer = index + tripletNeighbourOffset1[i];
            coordsCell3Buffer = index + tripletNeighbourOffset2[i];
            a.handleCellTriplet(_cells[coordsCell1Buffer], _cells[coordsCell2Buffer], _cells[coordsCell3Buffer],
                                coordsCell1Buffer, coordsCell2Buffer, coordsCell3Buffer);
          }
#else
          }
#endif
        } // coords(0)
#if (MD_DIM > 1)
      }
#endif
#if (MD_DIM > 2)
    }
#endif

#if (MD_OPENMP == MD_YES)
  }
#endif

  // end iteration();
  a.endCellIteration();
}

template <class A> void simplemd::services::LinkedCellService::iterateCellTriplets(A& a, const bool& useOpenMP) const {
  const tarch::la::Vector<MD_DIM, unsigned int> tripletIterationStart(0);
  const tarch::la::Vector<MD_DIM, unsigned int> tripletIterationLength(getLocalNumberOfCells() + getLocalIndexOfFirstCell());
  iterateCellTriplets(a, tripletIterationStart, tripletIterationLength, useOpenMP);
}

#endif // _MOLECULARDYNAMICS_SERVICES_LINKEDCELLSERVICE_H_
