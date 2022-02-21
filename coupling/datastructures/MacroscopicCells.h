// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_DATASTRUCTURES_MACROSCOPICCELLS_H_
#define _MOLECULARDYNAMICS_COUPLING_DATASTRUCTURES_MACROSCOPICCELLS_H_

#include "coupling/CouplingMDDefinitions.h"
#include "coupling/IndexConversion.h"
#include "coupling/datastructures/MacroscopicCell.h"
#include "coupling/interface/MDSolverInterface.h"

namespace coupling {
namespace datastructures {
template <class LinkedCell, unsigned int dim> class MacroscopicCells;
}
} // namespace coupling

/**
 *	@brief provides access to the macroscopic cells.
 *	@tparam LinkedCell linked cells that build up the
 *MacroscopicCellWithLinkedCells
 *	@tparam dim Number of dimensions; it can be 1, 2 or 3
 *  @author Philipp Neumann
 */
template <class LinkedCell, unsigned int dim>
class coupling::datastructures::MacroscopicCells {
public:
  /** Constructor: initialises the macroscopic cell
   *	@param numberLinkedCellsPerMacroscopicCell
   * 	@param indexConversion
   * 	@param mdSolverInterface
   */
  MacroscopicCells(
      tarch::la::Vector<dim, unsigned int> numberLinkedCellsPerMacroscopicCell,
      const coupling::IndexConversion<dim> &indexConversion,
      coupling::interface::MDSolverInterface<LinkedCell, dim>
          *mdSolverInterface);
  /** Destructor */
  ~MacroscopicCells();

  /** returns the pointer to the macroscopic cells with access to linked cell
   * structur. */
  coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell, dim> *
  getMacroscopicCellsWithLinkedCells();
  /** returns vector-of-pointers to macroscopic cells without access to linked
   * cells. We use this structure for data exchange between macroscopic and MD
   * solver.
   */
  const std::vector<coupling::datastructures::MacroscopicCell<dim> *> &
  getMacroscopicCells() const;

  /** apply the function apply(MacroscopicCell&,const unsigned int&) of a
   * generic class A to all local non-ghost macroscopic cells. This
   * functionality is used from the MacroscopicCellService to apply various
   * functionalities (such as momentum transfer) to the macroscopic cells which
   * cover the MD domain.
   */
  template <class A>
  void applyToLocalNonGhostMacroscopicCellsWithLinkedCells(A &a);
  /** apply the function apply(MacroscopicCell&,const unsigned int&) of a
   *generic class A to all local ghost macroscopic cells.
   *	@tparam A
   * 	@param a */
  template <class A>
  void applyToLocalGhostMacroscopicCellsWithLinkedCells(A &a);
  /** apply the function apply(MacroscopicCell&,const unsigned int&) of a
   * generic class A to all local macroscopic cells. We can have the same
   * traversal by applying both applyToLocalGhostMacroscopicCellsWithLinkedCells
   * and applyToLocalNonGhostMacroscopicCellsWithLinkedCells. However, this
   * implementation is more efficient.
   */
  template <class A> void applyToAllLocalMacroscopicCellsWithLinkedCells(A &a);
  /** apply the function apply(MacroscopicCell&,const unsigned int&) of a
   *generic class A to the first layer of GLOBAL non-ghost cells. If a cell is a
   *local ghost cell (due to parallelization, not due to being located outside
   *MD domain!), it will not be handled by this method. If a cell is a local
   *non-ghost cell and is not located at the boundary of the MD domain, it will
   *not be handled by this method. We use this traversal, e.g., for applying
   *boundary forces to molecules close to the outer boundary.
   *	@tparam A
   * 	@param a
   */
  template <class A>
  void applyToFirstLayerOfGlobalNonGhostCellsWithLinkedCells(A &a);

  /** \todo Helene!!
   */
  template <class A>
  void
  applyXLayersOfGlobalNonGhostCellsWithLinkedCells(A &a,
                                                   unsigned int layers2Use);

private:
  /** initialises the macroscopic cells: creates the buffer for the cells and
   * embeds linked cells into the macroscopic cells.
   * 	@param numberLinkedCellsPerMacroscopicCell
   * 	@param indexConversion
   * 	@param mdSolverInterface
   */
  coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell, dim> *
  initMacroscopicCellsWithLinkedCells(
      tarch::la::Vector<dim, unsigned int> numberLinkedCellsPerMacroscopicCell,
      const coupling::IndexConversion<dim> &indexConversion,
      coupling::interface::MDSolverInterface<LinkedCell, dim>
          *mdSolverInterface) const;
  /** initialises the macroscopic cells (without linked cells). This method
   * needs to be used in the constructor AFTER initialising the
   * _macroscopicCellsWithLinkedCells.
   * 	@param indexConversion
   */
  std::vector<coupling::datastructures::MacroscopicCell<dim> *>
  initMacroscopicCells(
      const coupling::IndexConversion<dim> &indexConversion) const;

  /** holds the macroscopic cells with linked cells. */
  coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell, dim>
      *_macroscopicCellsWithLinkedCells;
  /** holds pointers to all macroscopic cells with linked cells, but without
   * access to linked cells. This is used for interfacing to send-recv
   * operations.
   */
  std::vector<coupling::datastructures::MacroscopicCell<dim> *>
      _macroscopicCells;
  /** needed for index conversion. */
  const coupling::IndexConversion<dim> &_indexConversion;
};
#include "MacroscopicCells.cpph"

#endif // _MOLECULARDYNAMICS_COUPLING_DATASTRUCTURES_MACROSCOPICCELLS_H_
