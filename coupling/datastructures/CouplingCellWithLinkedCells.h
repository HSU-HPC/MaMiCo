#ifndef _MOLECULARDYNAMICS_COUPLING_DATASTRUCTURES_COUPLINGCELLWITHLINKEDCELLS_H_
#define _MOLECULARDYNAMICS_COUPLING_DATASTRUCTURES_COUPLINGCELLWITHLINKEDCELLS_H_

#include "CouplingCell.h"

namespace coupling {
namespace datastructures {
template <class LinkedCell, unsigned int dim> class CouplingCellWithLinkedCells;
} // namespace datastructures
} // namespace coupling

/** describes a quadratic/ cubic coupling cell filled with fluid. It is built
 *up by a certain number of linked cells (from the MD algorithm). The linked
 *  cells need to exactly fill this cell; no overlap/ non-fitting boundaries
 *  shall be allowed.
 *  We can use the CouplingCell-structure to evaluate macroscopic quantities
 *  over a certain MD volume and then map macroscopic conserved quantities
 *  between macro- and microscopic simulations.
 *	@brief defines the cell type with cell-averaged quantities. Derived from
 *the class coupling::datastructures::CouplingCell
 *	@tparam LinkedCell linked cells that build up the
 *CouplingCellWithLinkedCells
 *	@tparam dim Number of dimensions; it can be 1, 2 or 3
 *  @author Philipp Neumann
 */
template <class LinkedCell, unsigned int dim> class coupling::datastructures::CouplingCellWithLinkedCells : public coupling::datastructures::CouplingCell<dim> {
public:
  /** Constructor: initialises the coupling cell based on the assumption of
   * having blockSize linked cells;
   *  @param blockSize represents the number of linked cells in all spatial
   * directions.
   */
  CouplingCellWithLinkedCells(tarch::la::Vector<dim, unsigned int> blockSize)
      : coupling::datastructures::CouplingCell<dim>(), _numberCells(getNumberCells(blockSize)), _linkedCells(NULL) {

    _linkedCells = new LinkedCell*[_numberCells];
    if (_linkedCells == NULL) {
      std::cout << "ERROR coupling::datastructures::CouplingCellWithLinkedCells: "
                   "_linkedCells == NULL"
                << std::endl;
      exit(EXIT_FAILURE);
    }
    // set each pointer to a NULL pointer
    for (unsigned int i = 0; i < _numberCells; i++) {
      _linkedCells[i] = NULL;
    }
  }
  /** Destructor */
  virtual ~CouplingCellWithLinkedCells() {
    if (_linkedCells != NULL) {
      for (unsigned int i = 0; i < _numberCells; i++) {
        _linkedCells[i] = NULL;
      }
      delete[] _linkedCells;
      _linkedCells = NULL;
    }
  }

  /** adds a linked cell to the coupling cell and puts it at position index.
   * We refer to the lexicographic ordering of the linked cells here.
   * @param cell the linked cell that should be inserted into the coupling
   * cell
   * @param index specifies the position, at which cell should be inserted
   */
  void addLinkedCell(LinkedCell& cell, const size_t index) { _linkedCells[index] = &cell; }

  /** This template fuction applies class A to all linked cells of this
   *coupling cell. The syntax is exactly the same as for regular cell
   *mappings so that a general cell mapping can also be directly applied to
   *single coupling cells only.
   *	@tparam A
   * 	@param a
   */
  template <class A> void iterateCells(A& a) {
    a.beginCellIteration();
    for (unsigned int i = 0; i < _numberCells; i++) {
      a.handleCell(*(_linkedCells[i]));
    }
    a.endCellIteration();
  }

  /** This template function applies class A to all linked cells of this
   *coupling cell. The syntax is exactly the same as for regular cell
   *mappings so that a general cell mapping can also be directly applied to
   *single coupling cells only. This method is const, i.e. id does not modify
   *anything except for the object a. This allows for further optimisations.
   *	@tparam A
   * 	@param a
   */
  template <class A> void iterateConstCells(A& a) const {
    a.beginCellIteration();
    for (unsigned int i = 0; i < _numberCells; i++) {
      a.handleCell(*(_linkedCells[i]));
    }
    a.endCellIteration();
  }

private:
  /** computes and returns the number of cells specified by the product of the
   * entries of blockSize in all dimensions
   *  @param blockSize represents the number of linked cells in all spatial
   * directions.
   */
  unsigned int getNumberCells(tarch::la::Vector<dim, unsigned int> blockSize) const {
    unsigned int num = 1;
    for (unsigned int d = 0; d < dim; d++) {
      num = num * blockSize[d];
    }
    return num;
  }

  /** total number of linked cells contained in this coupling cell */
  const unsigned int _numberCells;

  /** holds pointers to all linked cells that describe the microscopic dynamics
   * within the coupling cell */
  LinkedCell** _linkedCells;
};

#endif // _MOLECULARDYNAMICS_COUPLING_DATASTRUCTURES_COUPLINGCELLWITHLINKEDCELLS_H_
