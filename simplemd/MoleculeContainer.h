#ifndef _MOLECULARDYNAMICS_MOLECULARCONTAINER_H_
#define _MOLECULARDYNAMICS_MOLECULARCONTAINER_H_

#include <Kokkos_Core.hpp>

#include <simplemd/MolecularDynamicsDefinitions.h>
#include <simplemd/Molecule.h>
#include <simplemd/LinkedCell.h>
#include <simplemd/services/ParallelTopologyService.h>
#include <simplemd/services/MolecularPropertiesService.h>

namespace simplemd {
namespace services {
// forward declarations to remove circular dependencies
class ParallelTopologyService;
} // namespace services

class MoleculeContainer;
} // namespace simplemd

/**
 * @brief Contains molecules, creates linked cells and manages its own memory.
 *
 * The MoleculeContainer class serves as the sole container responsible for all molecule-related work.
 * Molecules are stored in a 2D Kokkos View, where each row corresponds to a linked cell (ghost included).
 * While creating a container, the maximum possible capacity (in molecules) of any cell is defined by cellCapacity.
 * Thus, each linked cell has a numMolecules() <= capacity.
 * All data with an index between numMolecules and capacity is treated as garbage data.
 */
class simplemd::MoleculeContainer {
public:
  /**
   * @brief Construct a new MoleculeContainer object
   *
   * Uses parallelTopologyService to store relevant information related to spatial layout,
   * so that 3D coordinated can be converted to 3D and 1D linked cell indices.
   *
   * @param parallelTopologyService Used to extract and store local number of cells, domain offset etc.
   * @param cellCapacity The maximum capacity of any cell. Cannot be changed, must provide ample room at compiletime.
   */
  MoleculeContainer(simplemd::services::ParallelTopologyService& parallelTopologyService, int cellCapacity);

  ~MoleculeContainer();
  /**
   * @brief Inserts a molecule into a specific linked cell.
   *
   * @param cellIdx The one-dimensional index of the linked cell to insert the molecule into (ghost included).
   * @param molecule The molecule to be inserted.
   */
  void insert(int cellIdx, simplemd::Molecule& molecule);

  /**
   * @brief Inserts a molecule into a specific linked cell.
   *
   * @param cellIdx The vector index of the linked cell to insert the molecule into (ghost included).
   * @param molecule The molecule to be inserted.
   */
  void insert(tarch::la::Vector<MD_DIM, unsigned int> cellIdx, simplemd::Molecule& molecule) { insert(vectorIndexToLinear(cellIdx), molecule); }

  /**
   * @brief Inserts a molecule into the container.
   *
   * The function calculates the appropriate linked cell to insert the molecule into.
   *
   * @param molecule The molecule to be inserted.
   */
  void insert(simplemd::Molecule& molecule);

  /**
   * @brief Removes a molecule, given its ID and linked cell index.
   *
   * This is done by swapping the last molecule of the linked cell into position moleculeIdx,
   * and decrementing the number of molecules in this cell.
   *
   * @param cellIdx 1D linked cell index of the molecule (ghost included).
   * @param moleculeIdx Index of the molecule within the linked cell.
   */
  void remove(int cellIdx, int moleculeIdx);

  /**
   * @brief Clears a linked cell.
   *
   * This function simply sets the number of molecules in the linked cell as 0.
   * All data in the linked cell is treated as garbage data, and insertion starts from index 0
   * within the linked cell again, overwriting previous data.
   *
   * @param cellIdx 1D index of the linked cell to clear (ghost included).
   */
  void clearLinkedCell(int cellIdx);

  /**
   * @brief Removes all outgoing molecules from a linked cell and moves them to the appropriate destination.
   *
   * @param cellIdx 1D index of the linked cell to sort (ghost included).
   */
  void sort(int cellIdx);

  /**
   * @brief Puts all molecules into their appropriate linked cells.
   *
   * Since a molecule can only travel to its immediate neighbours per timestep, sorting is done with a red-black
   * traversal of the whole domain (including ghost cells). This ensures that there are no concurrency issues since
   * no two linked cells should try to write to the same linked cell. Concurrency is implemented in quarter shells.
   */
  void sort();

  /**
   * @brief Get the Molecule at linked cell index i (ghost included), position j
   *
   * @param i
   * @param j
   * @return simplemd::Molecule&
   */
  simplemd::Molecule& getMoleculeAt(int i, int j) const;

  /**
   * @brief Returns the linked cell at 1D index idx (ghost included)
   *
   * @param idx
   * @return simplemd::LinkedCell
   */
  simplemd::LinkedCell operator[](unsigned int idx) const;

  /**
   * @brief Returns the linked cell at 3D index idx (ghost included)
   *
   * @param idx
   * @return simplemd::LinkedCell
   */
  simplemd::LinkedCell operator[](tarch::la::Vector<MD_DIM, unsigned int> cellIdx) const;

  /**
   * @brief Get the total number of cells in the container
   *
   * @return int
   */
  int getLocalNumberOfCellsScalarWithGhost() const;

  /**
   * @brief Returns the number of molecules in all cells
   *
   * @return const size_t
   */
  const size_t getLocalNumberOfMoleculesWithGhost() const;

  /**
   * @brief returns the index of the first (non-ghost) cell along each dimension
   *
   * @return tarch::la::Vector
   */
  const tarch::la::Vector<MD_DIM, unsigned int>& getLocalIndexOfFirstCell() const;

  /**
   * @brief returns the number of (non-ghost) cells along each dimension
   *
   * @return tarch::la::Vector
   */
  const tarch::la::Vector<MD_DIM, unsigned int> getLocalNumberOfCells() const;

  /**
   * @brief returns the local cell index vector for the local cell index cellIndex
   *
   * @return tarch::la::Vector
   */
  tarch::la::Vector<MD_DIM, unsigned int> getLocalCellIndexVector(const unsigned int cellIndex) const;

  /**
   * @brief can be used to apply a molecule-mapping which is iterated over all molecules of this process
   * uses static member in mapping class (A::IsParallel) to determine whether the parallel or serial iterator will be called
   */
  template <class A> void iterateMolecules(A& a);

  /**
   * @brief iterates over all cells in the range defined by the lower left front
   * corner cell lowerLeftFrontCell and the size of the domain cellRange.
   * cellRange defines a number of cells in each spatial direction that the
   * class A shall be applied to. lowerLeftFrontCell needs to be given in local
   * coordinates.
   *
   * uses static member in mapping class (A::IsParallel) to determine whether the parallel or serial iterator will be called.
   */
  template <class A>
  void iterateCells(A& a, const tarch::la::Vector<MD_DIM, unsigned int>& lowerLeftFrontCell, const tarch::la::Vector<MD_DIM, unsigned int>& cellRange);

  /**
   * @brief iterates over all cells in the inner part (i.e. does not consider the
   * ghost layer)
   *
   * uses static member in mapping class (A::IsParallel) to determine whether the parallel or serial iterator will be called
   */
  template <class A> void iterateCells(A& a);

  /**
   * @brief iterates over all cell pairs for the cells in the inner region of each
   * local process
   *
   * uses static member in mapping class (A::IsParallel) to determine whether the parallel or serial iterator will be called
   */
  template <class A> void iterateCellPairs(A& a);

  /**
   * @brief iterates over all cell pairs cell1 and cell2 with cell1 in the range
   * described by lowerLeftFrontCell and cellRange; cell2 does not need to lie
   * within the range (example: iterate only over lowerLeftFrontCell=(1,1,1) and
   * cellRange=(1,1,1). Then, we will consider amongst others the pair
   * (0,0,0),(1,1,1)).
   *
   * uses static member in mapping class (A::IsParallel) to determine whether the parallel or serial iterator will be called
   */
  template <class A>
  void iterateCellPairs(A& a, const tarch::la::Vector<MD_DIM, unsigned int>& lowerLeftFrontCell, const tarch::la::Vector<MD_DIM, unsigned int>& cellRange);

private:
  /**
   * @brief applies molecule mapping without any node-level parallelisation
   */
  template <class A> void iterateMoleculesSerial(A& a);

  /**
   * @brief applies molecule mapping while parallelising using Kokkos
   */
  template <class A> void iterateMoleculesParallel(A& a);

  /**
   * @brief iterates over cells in parallel using Kokkos
   */
  template <class A>
  void iterateCellsParallel(A& a, const tarch::la::Vector<MD_DIM, unsigned int>& lowerLeftFrontCell, const tarch::la::Vector<MD_DIM, unsigned int>& cellRange);

  /**
   * @brief iterates over all cells in the inner part (i.e. does not consider the
   * ghost layer) in parallel using Kokkos
   */
  template <class A> void iterateCellsParallel(A& a);

  /**
   * @brief iterates over all cell pairs for the cells in the inner region of each
   * local process in parallel using Kokkos
   */
  template <class A> void iterateCellPairsParallel(A& a);

  /**
   * @brief iterates over all cell pairs cell1 and cell2 with cell1 in the range
   * described by lowerLeftFrontCell and cellRange in parallel using Kokkos; cell2 does not need to lie
   * within the range (example: iterate only over lowerLeftFrontCell=(1,1,1) and
   * cellRange=(1,1,1). Then, we will consider amongst others the pair
   * (0,0,0),(1,1,1)).
   */
  template <class A>
  void iterateCellPairsParallel(A& a, const tarch::la::Vector<MD_DIM, unsigned int>& lowerLeftFrontCell,
                                const tarch::la::Vector<MD_DIM, unsigned int>& cellRange);

  /**
   * @brief iterates over cells without parallelization
   */
  template <class A>
  void iterateCellsSerial(A& a, const tarch::la::Vector<MD_DIM, unsigned int>& lowerLeftFrontCell, const tarch::la::Vector<MD_DIM, unsigned int>& cellRange);

  /**
   * @brief iterates over all cells in the inner part (i.e. does not consider the
   * ghost layer) without parallelisation */
  template <class A> void iterateCellsSerial(A& a);

  /**
   * @brief iterates over all cell pairs for the cells in the inner region of each
   * local process without parallelisation */
  template <class A> void iterateCellPairsSerial(A& a);

  /**
   * @brief iterates over all cell pairs cell1 and cell2 with cell1 in the range
   * described by lowerLeftFrontCell and cellRange without parallelisation; cell2 does not need to lie
   * within the range (example: iterate only over lowerLeftFrontCell=(1,1,1) and
   * cellRange=(1,1,1). Then, we will consider amongst others the pair
   * (0,0,0),(1,1,1)).
   */
  template <class A>
  void iterateCellPairsSerial(A& a, const tarch::la::Vector<MD_DIM, unsigned int>& lowerLeftFrontCell,
                              const tarch::la::Vector<MD_DIM, unsigned int>& cellRange);

  /**
   * @brief Converts a global 3D spatial coordinate to a local 1D linked cell index.
   *
   * This can be used to find the linked cell that a molecule belongs to.
   *
   * @param position 3D spatial coordinate.
   * @return unsigned int
   */
  unsigned int positionToCellIndex(const tarch::la::Vector<MD_DIM, double>& position) const;

  /**
   * @brief Converts a 3D local linked cell index into a 1D local linked cell index.
   *
   * @param vectorIndex 3D local index of the linked cell.
   * @return const unsigned int
   */
  const unsigned int vectorIndexToLinear(const tarch::la::Vector<MD_DIM, unsigned int>& vectorIndex) const;

  /**
   * @brief returns true if the local cell index cellIndex describes a linked cell
   * within the ghost layer
   *
   * @param cellIndex The linear cell index
   * @return bool
   */
  bool isGhostCell(const size_t cellIndex) const;

  /** number of cells per direction in the local domain */
  const tarch::la::Vector<MD_DIM, unsigned int> _numCells;
  /** The number of ghost cells around the local domain along each axis on each side*/
  const tarch::la::Vector<MD_DIM, unsigned int> _ghostCellLayerThickness;
  const tarch::la::Vector<MD_DIM, unsigned int> _numLocalCellsNoGhost;

  /** maximum number of particles a cell (a row of the view) can contain
   * if this is exceeded when writing to cell, the simulation behaviour is undefined
   */
  int _cellCapacity;

#if (MD_ERROR == MD_YES)
  inline void checkOperationWouldExceedCapacity(int sizePostOp) const;

  /** domain size */
  const tarch::la::Vector<MD_DIM, double> _domainSize;
#endif

  /** global domain offset */
  const tarch::la::Vector<MD_DIM, double> _domainOffset;

  /** mesh width of the linked cells */
  const tarch::la::Vector<MD_DIM, double> _meshWidth;

  /** global index of the first cell of this domain */
  const tarch::la::Vector<MD_DIM, unsigned int> _globalIndexOfFirstCell;

  /** local index of the first cell within this domain */
  const tarch::la::Vector<MD_DIM, unsigned int> _localIndexOfFirstCell;

  Kokkos::View<simplemd::Molecule**, Kokkos::LayoutRight, Kokkos::SharedSpace> _moleculeData;
  Kokkos::View<size_t*, Kokkos::LayoutRight, Kokkos::SharedSpace> _linkedCellNumMolecules;
  Kokkos::View<bool*, Kokkos::LayoutRight, Kokkos::SharedSpace> _linkedCellIsGhostCell;
};

template <class A> void simplemd::MoleculeContainer::iterateMolecules(A& a) {
  if constexpr (A::IsParallel) {
    iterateMoleculesParallel(a);
  } else {
    iterateMoleculesSerial(a);
  }
}

template <class A> void simplemd::MoleculeContainer::iterateMoleculesSerial(A& a) {
  a.beginMoleculeIteration();
  for (unsigned int i = 0; i < _linkedCellNumMolecules.size(); i++) {
    for (unsigned int j = 0; j < _linkedCellNumMolecules(i); j++) {
#if (MD_DEBUG == MD_YES)
      std::cout << "Handle molecule " << j << " in cell #" << i << std::endl;
#endif
      a.handleMolecule(getMoleculeAt(i, j));
    }
  }
  a.endMoleculeIteration();
}

template <class A> void simplemd::MoleculeContainer::iterateMoleculesParallel(A& a) {
#if (MD_DEBUG == MD_YES)
  iterateMoleculesSerial(a);
#else
  a.beginMoleculeIteration();
  Kokkos::parallel_for(
      _linkedCellNumMolecules.size(), KOKKOS_LAMBDA(const unsigned int i) {
        for (unsigned int j = 0; j < _linkedCellNumMolecules(i); j++) {
          a.handleMolecule(getMoleculeAt(i, j));
        }
      });
  Kokkos::fence(); // Ensure results are available on the host
  a.endMoleculeIteration();
#endif
}

template <class A>
void simplemd::MoleculeContainer::iterateCellsSerial(A& a, const tarch::la::Vector<MD_DIM, unsigned int>& lowerLeftFrontCell,
                                                     const tarch::la::Vector<MD_DIM, unsigned int>& cellRange) {
  unsigned int index = 0;
#if (MD_ERROR == MD_YES)
  for (unsigned int d = 0; d < MD_DIM; d++) {
    if (cellRange[d] == 0) {
      std::cout << "ERROR simplemd::MoleculeContainer::iterateCells: "
                   "cellRange("
                << d << ")==0!" << std::endl;
      exit(EXIT_FAILURE);
    }
    if (lowerLeftFrontCell[d] + cellRange[d] > 2 * _ghostCellLayerThickness[d] + _numLocalCellsNoGhost[d]) {
      std::cout << "ERROR simplemd::MoleculeContainer::iterateCells(): "
                   "defined Range does not fit into local sub-domain!"
                << std::endl;
      exit(EXIT_FAILURE);
    }
  }
#endif

  // start iteration();
  a.beginCellIteration();
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
        index = vectorIndexToLinear(coords);
        simplemd::LinkedCell cell = (*this)[index];
        a.handleCell(cell);
      }
#if (MD_DIM > 1)
    }
#endif
#if (MD_DIM > 2)
  }
#endif

  // end iteration();
  a.endCellIteration();
}

template <class A> void simplemd::MoleculeContainer::iterateCellsSerial(A& a) { iterateCellsSerial(a, _ghostCellLayerThickness, _numLocalCellsNoGhost); }

template <class A>
void simplemd::MoleculeContainer::iterateCellPairsSerial(A& a, const tarch::la::Vector<MD_DIM, unsigned int>& lowerLeftFrontCell,
                                                         const tarch::la::Vector<MD_DIM, unsigned int>& cellRange) {
  tarch::la::Vector<MD_LINKED_CELL_NEIGHBOURS / 2, unsigned int> neighbourOffset;
  tarch::la::Vector<MD_LINKED_CELL_NEIGHBOURS / 2, unsigned int> indexOffset;
#if (MD_DIM == 1)
  neighbourOffset[0] = 1;
  indexOffset[0] = 0;
#elif (MD_DIM == 2)
  indexOffset[0] = 0;
  neighbourOffset[0] = 1;
  indexOffset[1] = 0;
  neighbourOffset[1] = _numLocalCellsNoGhost[0] + 2;
  indexOffset[2] = 0;
  neighbourOffset[2] = _numLocalCellsNoGhost[0] + 3;
  indexOffset[3] = 1;
  neighbourOffset[3] = _numLocalCellsNoGhost[0] + 2;
#elif (MD_DIM == 3)
  indexOffset[0] = 0;
  neighbourOffset[0] = 1;
  indexOffset[1] = 0;
  neighbourOffset[1] = _numLocalCellsNoGhost[0] + 2;
  indexOffset[2] = 0;
  neighbourOffset[2] = _numLocalCellsNoGhost[0] + 3;
  indexOffset[3] = 0;
  neighbourOffset[3] = (_numLocalCellsNoGhost[0] + 2) * (_numLocalCellsNoGhost[1] + 2);
  indexOffset[4] = 0;
  neighbourOffset[4] = (_numLocalCellsNoGhost[0] + 2) * (_numLocalCellsNoGhost[1] + 2) + 1;
  indexOffset[5] = 0;
  neighbourOffset[5] = (_numLocalCellsNoGhost[0] + 2) * (_numLocalCellsNoGhost[1] + 2) + (_numLocalCellsNoGhost[0] + 2);
  indexOffset[6] = 0;
  neighbourOffset[6] = (_numLocalCellsNoGhost[0] + 2) * (_numLocalCellsNoGhost[1] + 2) + (_numLocalCellsNoGhost[0] + 2) + 1;

  indexOffset[7] = 1;
  neighbourOffset[7] = _numLocalCellsNoGhost[0] + 2;
  indexOffset[8] = 1;
  neighbourOffset[8] = (_numLocalCellsNoGhost[0] + 2) * (_numLocalCellsNoGhost[1] + 2);
  indexOffset[9] = 1;
  neighbourOffset[9] = (_numLocalCellsNoGhost[0] + 2) * (_numLocalCellsNoGhost[1] + 2) + (_numLocalCellsNoGhost[0] + 2);

  indexOffset[10] = _numLocalCellsNoGhost[0] + 2;
  neighbourOffset[10] = (_numLocalCellsNoGhost[0] + 2) * (_numLocalCellsNoGhost[1] + 2);
  indexOffset[11] = _numLocalCellsNoGhost[0] + 2;
  neighbourOffset[11] = (_numLocalCellsNoGhost[0] + 2) * (_numLocalCellsNoGhost[1] + 2) + 1;

  indexOffset[12] = (_numLocalCellsNoGhost[0] + 2) + 1;
  neighbourOffset[12] = (_numLocalCellsNoGhost[0] + 2) * (_numLocalCellsNoGhost[1] + 2);
#endif

#if (MD_ERROR == MD_YES)
  for (unsigned int d = 0; d < MD_DIM; d++) {
    if (cellRange[d] > 2 * _ghostCellLayerThickness[d] + _numLocalCellsNoGhost[d] - 1) {
      std::cout << "ERROR simplemd::MoleculeContainer::iterateCellPairs(): "
                   "defined Range does not fit into local sub-domain!"
                << std::endl;
      exit(EXIT_FAILURE);
    }
  }
#endif

  // start iteration();
  a.beginCellIteration();
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
        index = vectorIndexToLinear(coords);
#if (MD_DEBUG == MD_YES)
        std::cout << "iterateCellPairs: Single index " << index << std::endl;
#endif
        simplemd::LinkedCell cell = (*this)[index];
        a.handleCell(cell);
        // handle pairs (lower,left,back-oriented cells)
        for (unsigned int i = 0; i < MD_LINKED_CELL_NEIGHBOURS / 2; i++) {
#if (MD_DEBUG == MD_YES)
          std::cout << "iterateCellPairs: Pair index " << index + indexOffset[i] << "," << index + neighbourOffset[i] << std::endl;
#endif

          coordsCell1Buffer = index + indexOffset[i];
          coordsCell2Buffer = index + neighbourOffset[i];
          simplemd::LinkedCell cell1 = (*this)[coordsCell1Buffer];
          simplemd::LinkedCell cell2 = (*this)[coordsCell2Buffer];
          a.handleCellPair(cell1, cell2, coordsCell1Buffer, coordsCell2Buffer);
        }
      } // coords(0)
#if (MD_DIM > 1)
    }
#endif
#if (MD_DIM > 2)
  }
#endif

  // end iteration();
  a.endCellIteration();
}

template <class A> void simplemd::MoleculeContainer::iterateCellPairsSerial(A& a) {
  const tarch::la::Vector<MD_DIM, unsigned int> pairIterationStart(0);
  const tarch::la::Vector<MD_DIM, unsigned int> pairIterationLength(getLocalNumberOfCells() + getLocalIndexOfFirstCell());
  iterateCellPairsSerial(a, pairIterationStart, pairIterationLength);
}

template <class A>
void simplemd::MoleculeContainer::iterateCellsParallel(A& a, const tarch::la::Vector<MD_DIM, unsigned int>& lowerLeftFrontCell,
                                                       const tarch::la::Vector<MD_DIM, unsigned int>& cellRange) {
#if (MD_ERROR == MD_YES)
  for (unsigned int d = 0; d < MD_DIM; d++) {
    if (cellRange[d] == 0) {
      std::cout << "ERROR simplemd::MoleculeContainer::iterateCells: "
                   "cellRange("
                << d << ")==0!" << std::endl;
      exit(EXIT_FAILURE);
    }
    if (lowerLeftFrontCell[d] + cellRange[d] > 2 * _ghostCellLayerThickness[d] + _numLocalCellsNoGhost[d]) {
      std::cout << "ERROR simplemd::MoleculeContainer::iterateCells(): "
                   "defined Range does not fit into local sub-domain!"
                << std::endl;
      exit(EXIT_FAILURE);
    }
  }
#endif

  // start iteration();
  a.beginCellIteration();
  /**
   * The size<> Vector stores the number of cells, plus ghost layer.
   * If there are (1,1,1) ghost cells per dimension, getLocalIndexOfFirstCell will return (1,1,1)
   * Thus this is multiplied by 2 to account for ghost cells in both locations (begin, end) per axis
   * and then added to local number of cells
   */
  const tarch::la::Vector<MD_DIM, unsigned int> size(getLocalNumberOfCells() + 2u * getLocalIndexOfFirstCell());
  const int length = cellRange[0]
#if (MD_DIM > 1)
                     * cellRange[1]
#endif
#if (MD_DIM > 2)
                     * cellRange[2]
#endif
      ;
  // loop over domain, but with a single loop
  MoleculeContainer& container = (*this);
  Kokkos::parallel_for(
      length, KOKKOS_LAMBDA(const unsigned int i) {
// compute index of the current cell
#if (MD_DIM > 1)
        int helpIndex1 = i;
        int helpIndex2 = 0;
#endif
        unsigned int index = 0;

#if (MD_DIM > 2)
        // determine plane within traversed block
        helpIndex2 = helpIndex1 / (cellRange[0] * cellRange[1]);
        // save rest of index in helpIndex1
        helpIndex1 = helpIndex1 - helpIndex2 * (cellRange[0] * cellRange[1]);
        // compute contribution to index
        index += (lowerLeftFrontCell[2] + helpIndex2) * size[0] * size[1];
#endif
#if (MD_DIM > 1)
        // determine plane within traversed block
        helpIndex2 = helpIndex1 / cellRange[0];
        // save rest of index in helpIndex1
        helpIndex1 = helpIndex1 - helpIndex2 * cellRange[0];
        // compute contribution to index
        index += (lowerLeftFrontCell[1] + helpIndex2) * size[0];
        // compute contribution for last dimension
        index += (lowerLeftFrontCell[0] + helpIndex1);
#else
      index = lowerLeftFrontCell[0] + i;
#endif
#if (MD_DEBUG == MD_YES)
        std::cout << "Handle cell " << index << std::endl;
#endif

        // handle cell
        a.handleCell(container[index]);
      });          // Kokkos::parallel_for
  Kokkos::fence(); // Ensure results are available on the host
  // end iteration();
  a.endCellIteration();
}

template <class A> void simplemd::MoleculeContainer::iterateCellsParallel(A& a) { iterateCellsParallel(a, _ghostCellLayerThickness, _numLocalCellsNoGhost); }

template <class A>
void simplemd::MoleculeContainer::iterateCellPairsParallel(A& a, const tarch::la::Vector<MD_DIM, unsigned int>& lowerLeftFrontCell,
                                                           const tarch::la::Vector<MD_DIM, unsigned int>& cellRange) {
  tarch::la::Vector<MD_LINKED_CELL_NEIGHBOURS / 2, unsigned int> neighbourOffset;
  tarch::la::Vector<MD_LINKED_CELL_NEIGHBOURS / 2, unsigned int> indexOffset;
#if (MD_DIM == 1)
  neighbourOffset[0] = 1;
  indexOffset[0] = 0;
#elif (MD_DIM == 2)
  indexOffset[0] = 0;
  neighbourOffset[0] = 1;
  indexOffset[1] = 0;
  neighbourOffset[1] = _numLocalCellsNoGhost[0] + 2;
  indexOffset[2] = 0;
  neighbourOffset[2] = _numLocalCellsNoGhost[0] + 3;
  indexOffset[3] = 1;
  neighbourOffset[3] = _numLocalCellsNoGhost[0] + 2;
#elif (MD_DIM == 3)
  indexOffset[0] = 0;
  neighbourOffset[0] = 1;
  indexOffset[1] = 0;
  neighbourOffset[1] = _numLocalCellsNoGhost[0] + 2;
  indexOffset[2] = 0;
  neighbourOffset[2] = _numLocalCellsNoGhost[0] + 3;
  indexOffset[3] = 0;
  neighbourOffset[3] = (_numLocalCellsNoGhost[0] + 2) * (_numLocalCellsNoGhost[1] + 2);
  indexOffset[4] = 0;
  neighbourOffset[4] = (_numLocalCellsNoGhost[0] + 2) * (_numLocalCellsNoGhost[1] + 2) + 1;
  indexOffset[5] = 0;
  neighbourOffset[5] = (_numLocalCellsNoGhost[0] + 2) * (_numLocalCellsNoGhost[1] + 2) + (_numLocalCellsNoGhost[0] + 2);
  indexOffset[6] = 0;
  neighbourOffset[6] = (_numLocalCellsNoGhost[0] + 2) * (_numLocalCellsNoGhost[1] + 2) + (_numLocalCellsNoGhost[0] + 2) + 1;

  indexOffset[7] = 1;
  neighbourOffset[7] = _numLocalCellsNoGhost[0] + 2;
  indexOffset[8] = 1;
  neighbourOffset[8] = (_numLocalCellsNoGhost[0] + 2) * (_numLocalCellsNoGhost[1] + 2);
  indexOffset[9] = 1;
  neighbourOffset[9] = (_numLocalCellsNoGhost[0] + 2) * (_numLocalCellsNoGhost[1] + 2) + (_numLocalCellsNoGhost[0] + 2);

  indexOffset[10] = _numLocalCellsNoGhost[0] + 2;
  neighbourOffset[10] = (_numLocalCellsNoGhost[0] + 2) * (_numLocalCellsNoGhost[1] + 2);
  indexOffset[11] = _numLocalCellsNoGhost[0] + 2;
  neighbourOffset[11] = (_numLocalCellsNoGhost[0] + 2) * (_numLocalCellsNoGhost[1] + 2) + 1;

  indexOffset[12] = (_numLocalCellsNoGhost[0] + 2) + 1;
  neighbourOffset[12] = (_numLocalCellsNoGhost[0] + 2) * (_numLocalCellsNoGhost[1] + 2);
#endif

#if (MD_ERROR == MD_YES)
  for (unsigned int d = 0; d < MD_DIM; d++) {
    if (cellRange[d] > 2 * _ghostCellLayerThickness[d] + _numLocalCellsNoGhost[d] - 1) {
      std::cout << "ERROR simplemd::MoleculeContainer::iterateCellPairs(): "
                   "defined Range does not fit into local sub-domain!"
                << std::endl;
      exit(EXIT_FAILURE);
    }
  }
#endif

  // start iteration();
  a.beginCellIteration();
  /**
   * The size<> Vector stores the number of cells, plus ghost layer.
   * If there are (1,1,1) ghost cells per dimension, getLocalIndexOfFirstCell will return (1,1,1)
   * Thus this is multiplied by 2 to account for ghost cells in both locations (begin, end) per axis
   * and then added to local number of cells
   */
  const tarch::la::Vector<MD_DIM, unsigned int> size(getLocalNumberOfCells() + 2u * getLocalIndexOfFirstCell());

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
        MoleculeContainer& container = (*this);
        Kokkos::parallel_for(
            length, KOKKOS_LAMBDA(const unsigned int j) {
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

              simplemd::LinkedCell cell = operator[](index);
              a.handleCell(cell);
              // handle pairs (lower,left,back-oriented cells)
              for (unsigned int i = 0; i < MD_LINKED_CELL_NEIGHBOURS / 2; i++) {
#if (MD_DEBUG == MD_YES)
                std::cout << "iterateCellPairs: Pair index " << index + indexOffset[i] << "," << index + neighbourOffset[i] << std::endl;
#endif
                coordsCell1Buffer = index + indexOffset[i];
                coordsCell2Buffer = index + neighbourOffset[i];
                auto cell1 = container[coordsCell1Buffer];
                auto cell2 = container[coordsCell2Buffer];
                a.handleCellPair(cell1, cell2, coordsCell1Buffer, coordsCell2Buffer);
              }
            });          // j, Kokkos::parallel_for
        Kokkos::fence(); // Ensure results are available on the host
      } // x
#if (MD_DIM > 1)
    } // y
#endif
#if (MD_DIM > 2)
  } // z
#endif
  // end iteration();
  a.endCellIteration();
}

template <class A> void simplemd::MoleculeContainer::iterateCellPairsParallel(A& a) {
  const tarch::la::Vector<MD_DIM, unsigned int> pairIterationStart(0);
  const tarch::la::Vector<MD_DIM, unsigned int> pairIterationLength(getLocalNumberOfCells() + getLocalIndexOfFirstCell());
  iterateCellPairsParallel(a, pairIterationStart, pairIterationLength);
}

template <class A>
void simplemd::MoleculeContainer::iterateCellPairs(A& a, const tarch::la::Vector<MD_DIM, unsigned int>& lowerLeftFrontCell,
                                                   const tarch::la::Vector<MD_DIM, unsigned int>& cellRange) {
  if constexpr (A::IsParallel) {
    iterateCellPairsParallel(a, lowerLeftFrontCell, cellRange);
  } else {
    iterateCellPairsSerial(a, lowerLeftFrontCell, cellRange);
  }
}

template <class A> void simplemd::MoleculeContainer::iterateCellPairs(A& a) {
  const tarch::la::Vector<MD_DIM, unsigned int> pairIterationStart(0);
  const tarch::la::Vector<MD_DIM, unsigned int> pairIterationLength(getLocalNumberOfCells() + getLocalIndexOfFirstCell());
  iterateCellPairs(a, pairIterationStart, pairIterationLength);
}

template <class A>
void simplemd::MoleculeContainer::iterateCells(A& a, const tarch::la::Vector<MD_DIM, unsigned int>& lowerLeftFrontCell,
                                               const tarch::la::Vector<MD_DIM, unsigned int>& cellRange) {
  if constexpr (A::IsParallel) {
    iterateCellsParallel(a, lowerLeftFrontCell, cellRange);
  } else {
    iterateCellsSerial(a, lowerLeftFrontCell, cellRange);
  }
}

template <class A> void simplemd::MoleculeContainer::iterateCells(A& a) { iterateCells(a, _ghostCellLayerThickness, _numLocalCellsNoGhost); }

#endif // _MOLECULARDYNAMICS_MOLECULARCONTAINER_H_
