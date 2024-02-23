// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_DATASTRUCTURES_COUPLINGCELLS_H_
#define _MOLECULARDYNAMICS_COUPLING_DATASTRUCTURES_COUPLINGCELLS_H_

#include "coupling/CouplingMDDefinitions.h"
#include "coupling/datastructures/CouplingCell.h"
#include "coupling/interface/MDSolverInterface.h"

namespace coupling {
namespace datastructures {
template <class LinkedCell, unsigned int dim> class CouplingCells;
}
} // namespace coupling

/**
 *	@brief provides access to the coupling cells.
 *	@tparam LinkedCell linked cells that build up the
 *CouplingCellWithLinkedCells
 *	@tparam dim Number of dimensions; it can be 1, 2 or 3
 *  @author Philipp Neumann
 */
template <class LinkedCell, unsigned int dim> class coupling::datastructures::CouplingCells {
public:
  /** Constructor: initialises the coupling cell
   *	@param numberLinkedCellsPerCouplingCell
   * 	@param mdSolverInterface
   */
  CouplingCells(tarch::la::Vector<dim, unsigned int> numberLinkedCellsPerCouplingCell,
                coupling::interface::MDSolverInterface<LinkedCell, dim>* mdSolverInterface);
  /** Destructor */
  ~CouplingCells();

  /** returns the pointer to the coupling cells with access to linked cell
   * structur. */
  coupling::datastructures::CouplingCellWithLinkedCells<LinkedCell, dim>* getCouplingCellsWithLinkedCells();
  /** returns vector-of-pointers to coupling cells without access to linked
   * cells. We use this structure for data exchange between macroscopic and MD
   * solver.
   */
  const std::vector<coupling::datastructures::CouplingCell<dim>*>& getCouplingCells() const;

  /** apply the function apply(CouplingCell&,const unsigned int&) of a
   * generic class A to all local non-ghost coupling cells. This
   * functionality is used from the CouplingCellService to apply various
   * functionalities (such as momentum transfer) to the coupling cells which
   * cover the MD domain.
   */
  template <class A> void applyToLocalNonGhostCouplingCellsWithLinkedCells(A& a);
  /** apply the function apply(CouplingCell&,const unsigned int&) of a
   *generic class A to all local ghost coupling cells.
   *	@tparam A
   * 	@param a */
  template <class A> void applyToLocalGhostCouplingCellsWithLinkedCells(A& a);
  /** apply the function apply(CouplingCell&,const unsigned int&) of a
   * generic class A to all local coupling cells. We can have the same
   * traversal by applying both applyToLocalGhostCouplingCellsWithLinkedCells
   * and applyToLocalNonGhostCouplingCellsWithLinkedCells. However, this
   * implementation is more efficient.
   */
  template <class A> void applyToAllLocalCouplingCellsWithLinkedCells(A& a);
  /** apply the function apply(CouplingCell&,const unsigned int&) of a
   *generic class A to the first layer of GLOBAL non-ghost cells. If a cell is a
   *local ghost cell (due to parallelization, not due to being located outside
   *MD domain!), it will not be handled by this method. If a cell is a local
   *non-ghost cell and is not located at the boundary of the MD domain, it will
   *not be handled by this method. We use this traversal, e.g., for applying
   *boundary forces to molecules close to the outer boundary.
   *	@tparam A
   * 	@param a
   */
  template <class A> void applyToFirstLayerOfGlobalNonGhostCellsWithLinkedCells(A& a);

  template <class A> void applyXLayersOfGlobalNonGhostCellsWithLinkedCells(A& a, unsigned int layers2Use);

private:
  /** initialises the coupling cells: creates the buffer for the cells and
   * embeds linked cells into the coupling cells.
   * 	@param numberLinkedCellsPerCouplingCell
   * 	@param mdSolverInterface
   */
  coupling::datastructures::CouplingCellWithLinkedCells<LinkedCell, dim>*
  initCouplingCellsWithLinkedCells(tarch::la::Vector<dim, unsigned int> numberLinkedCellsPerCouplingCell,
                                   coupling::interface::MDSolverInterface<LinkedCell, dim>* mdSolverInterface) const;
  /** initialises the coupling cells (without linked cells). This method
   * needs to be used in the constructor AFTER initialising the
   * _couplingCellsWithLinkedCells.
   */
  std::vector<coupling::datastructures::CouplingCell<dim>*> initCouplingCells() const;

  /** holds the coupling cells with linked cells. */
  coupling::datastructures::CouplingCellWithLinkedCells<LinkedCell, dim>* _couplingCellsWithLinkedCells;
  /** holds pointers to all coupling cells with linked cells, but without
   * access to linked cells. This is used for interfacing to send-recv
   * operations.
   */
  std::vector<coupling::datastructures::CouplingCell<dim>*> _couplingCells;
};
#include "CouplingCells.cpph"

#endif // _MOLECULARDYNAMICS_COUPLING_DATASTRUCTURES_COUPLINGCELLS_H_
