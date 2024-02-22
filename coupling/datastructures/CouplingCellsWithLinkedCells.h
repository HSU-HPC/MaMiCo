#ifndef _MOLECULARDYNAMICS_COUPLING_DATASTRUCTURES_COUPLINGCELLSWITHLINKEDCELLS_H_
#define _MOLECULARDYNAMICS_COUPLING_DATASTRUCTURES_COUPLINGCELLSWITHLINKEDCELLS_H_

#include "CouplingCells.h"

namespace coupling {
namespace datastructures {
template <class LinkedCell, unsigned int dim> class CouplingCellsWithLinkedCells;
} // namespace datastructures
} // namespace coupling

/**
 *	@brief provides access to the coupling cells with linked cells.
 *	@tparam LinkedCell linked cells that build up the
 *CouplingCellWithLinkedCells
 *	@tparam dim Number of dimensions; it can be 1, 2 or 3
 */
template <class LinkedCell, unsigned int dim>
class coupling::datastructures::CouplingCellsWithLinkedCells : public coupling::datastructures::CouplingCells<dim> {
public:
  /** Constructor: initialises the coupling cell with linked cells
   *	@param numberLinkedCellsPerCouplingCell
   * 	@param indexConversion
   * 	@param mdSolverInterface
   */
  CouplingCellsWithLinkedCells(tarch::la::Vector<dim, unsigned int> numberLinkedCellsPerCouplingCell, const coupling::IndexConversion<dim>& indexConversion,
                               coupling::interface::MDSolverInterface<LinkedCell, dim>* mdSolverInterface);
  /** Destructor */
  ~CouplingCellsWithLinkedCells();

  /** returns the pointer to the coupling cells with access to linked cell
   * structur. */
  coupling::datastructures::CouplingCellWithLinkedCells<LinkedCell, dim>* getCouplingCellsWithLinkedCells();

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
   * 	@param indexConversion
   * 	@param mdSolverInterface
   */
  coupling::datastructures::CouplingCellWithLinkedCells<LinkedCell, dim>*
  initCouplingCellsWithLinkedCells(tarch::la::Vector<dim, unsigned int> numberLinkedCellsPerCouplingCell, const coupling::IndexConversion<dim>& indexConversion,
                                   coupling::interface::MDSolverInterface<LinkedCell, dim>* mdSolverInterface) const;

  /** initialises the coupling cells (without linked cells).
   * 	@param indexConversion
   */
  std::vector<coupling::datastructures::CouplingCell<dim>*> initCouplingCells(const coupling::IndexConversion<dim>& indexConversion) const;
  /** holds the coupling cells with linked cells. */
  coupling::datastructures::CouplingCellWithLinkedCells<LinkedCell, dim>* _couplingCellsWithLinkedCells;
};

#include "CouplingCellsWithLinkedCells.cpph"

#endif // _MOLECULARDYNAMICS_COUPLING_DATASTRUCTURES_COUPLINGCELLSWITHLINKEDCELLS_H_
