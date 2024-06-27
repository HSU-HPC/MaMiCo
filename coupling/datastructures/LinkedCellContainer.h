#ifndef _MOLECULARDYNAMICS_COUPLING_DATASTRUCTURES_COUPLINGCELLSWITHLINKEDCELLS_H_
#define _MOLECULARDYNAMICS_COUPLING_DATASTRUCTURES_COUPLINGCELLSWITHLINKEDCELLS_H_

#include "CellContainer.h"
#include "coupling/datastructures/CouplingCellWithLinkedCells.h"
#include "coupling/interface/MDSolverInterface.h"

namespace coupling {
namespace datastructures {
template <class LinkedCell, unsigned int dim> class LinkedCellContainer;
} // namespace datastructures
} // namespace coupling

/**
 *	@brief provides access to the coupling cells with linked cells.
 *	@tparam LinkedCell linked cells that build up the
 *CouplingCellWithLinkedCells
 *	@tparam dim Number of dimensions; it can be 1, 2 or 3
 */
template <class LinkedCell, unsigned int dim>
class coupling::datastructures::LinkedCellContainer : public coupling::datastructures::CellContainer<I02, dim> {
public:
  /** Constructor: initialises the coupling cell with linked cells
   *	@param numberLinkedCellsPerCouplingCell
   * 	@param mdSolverInterface
   */
  LinkedCellContainer(tarch::la::Vector<dim, unsigned int> numberLinkedCellsPerCouplingCell,
                               coupling::interface::MDSolverInterface<LinkedCell, dim>* mdSolverInterface);
  /** Destructor */
  ~LinkedCellContainer();

  /** returns the pointer to the coupling cells with access to linked cell
   * structur. */
  coupling::datastructures::CouplingCellWithLinkedCells<LinkedCell, dim>* getLinkedCellContainer();

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

  /** \todo Helene!!
   */
  template <class A> void applyXLayersOfGlobalNonGhostCellsWithLinkedCells(A& a, unsigned int layers2Use);

private:
  /** initialises the coupling cells: creates the buffer for the cells and
   * embeds linked cells into the coupling cells.
   * 	@param numberLinkedCellsPerCouplingCell
   * 	@param mdSolverInterface
   */
  coupling::datastructures::CouplingCellWithLinkedCells<LinkedCell, dim>*
  initLinkedCellContainer(tarch::la::Vector<dim, unsigned int> numberLinkedCellsPerCouplingCell,
                                   coupling::interface::MDSolverInterface<LinkedCell, dim>* mdSolverInterface) const;

  /** initialises the coupling cells (without linked cells).
   */
  std::vector<coupling::datastructures::CouplingCell<dim>*> initCouplingCells() const;
  /** holds the coupling cells with linked cells. */
  coupling::datastructures::CouplingCellWithLinkedCells<LinkedCell, dim>* _couplingCellsWithLinkedCells;
};

#include "LinkedCellContainer.cpph"

#endif // _MOLECULARDYNAMICS_COUPLING_DATASTRUCTURES_COUPLINGCELLSWITHLINKEDCELLS_H_
