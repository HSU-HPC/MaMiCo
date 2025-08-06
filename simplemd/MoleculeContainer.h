#ifndef _MOLECULARDYNAMICS_MOLECULARCONTAINER_H_
#define _MOLECULARDYNAMICS_MOLECULARCONTAINER_H_

#include <Kokkos_Core.hpp>

#include <simplemd/MolecularDynamicsDefinitions.h>
#include <simplemd/Molecule.h>
#include <simplemd/LinkedCell.h>
#include <simplemd/services/ParallelTopologyService.h>

namespace simplemd {
class MoleculeContainer;
}

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
  MoleculeContainer(simplemd::services::ParallelTopologyService parallelTopologyService, int cellCapacity);

  /**
   * @brief Inserts a molecule into a specific linked cell.
   *
   * @param cellIdx The one-dimensional index of the linked cell to insert the molecule into (ghost included).
   * @param molecule The molecule to be inserted.
   */
  void insert(int cellIdx, simplemd::Molecule& molecule);

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
   * @brief Returns the linked cell at index idx (ghost included)
   *
   * @param idx
   * @return simplemd::LinkedCell
   */
  simplemd::LinkedCell operator[](unsigned int idx);

  /**
   * @brief Get the total number of cells in the container
   *
   * @return int
   */
  int getNumCells() const;

private:
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
   * @brief Returns the number of molecules in all cells
   * 
   * @return const size_t
   */
  const size_t getNumberMolecules () const;

  /** number of cells per direction in the local domain */
  tarch::la::Vector<MD_DIM, unsigned int> _numCells;

  /** maximum number of particles a cell (a row of the view) can contain
   * if this is exceeded when writing to cell, the simulation is stopped
   */
  int _cellCapacity;

#if (MD_ERROR == MD_YES)
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
};
#endif // _MOLECULARDYNAMICS_MOLECULARCONTAINER_H_
