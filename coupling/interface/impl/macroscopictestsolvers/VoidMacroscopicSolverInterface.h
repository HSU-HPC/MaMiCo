// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_INTERFACE_VOIDMACROSCOPICSOLVERINTERFACE_H_
#define _MOLECULARDYNAMICS_COUPLING_INTERFACE_VOIDMACROSCOPICSOLVERINTERFACE_H_

#include "coupling/interface/impl/macroscopictestsolvers/TestMacroscopicSolverInterface.h"

namespace coupling {
namespace interface {
template <unsigned int dim> class VoidMacroscopicSolverInterface;
}
} // namespace coupling

/** no operations in send/recv process.
 *  @author Philipp Neumann
 */
template <unsigned int dim> class coupling::interface::VoidMacroscopicSolverInterface : public coupling::interface::TestMacroscopicSolverInterface<dim> {
public:
  VoidMacroscopicSolverInterface() : TestMacroscopicSolverInterface<dim>() {}
  virtual ~VoidMacroscopicSolverInterface() {}

  /** returns true if the cell at position globalCellIndex shall be received
   * from the MD solver. This function does not send the cell, but only steers
   * the process.
   */
  virtual bool receiveMacroscopicQuantityFromMDSolver(tarch::la::Vector<dim, unsigned int> globalCellIndex) { return false; }

  /** returns true if the cell at position globalCellIndex shall be sent to the
   * MD solver. This function does not send the cell, but only steers the
   * process.
   */
  virtual bool sendMacroscopicQuantityToMDSolver(tarch::la::Vector<dim, unsigned int> globalCellIndex) { return false; }

  /** returns the ranks on which the macroscopic solver holds/requires data of
   * the macroscopic cell at index globalCellIndex.
   */
  virtual std::vector<unsigned int> getRanks(tarch::la::Vector<dim, unsigned int> globalCellIndex) {
    std::vector<unsigned int> result;
    result.push_back(0);
    return result;
  }

  virtual std::vector<coupling::datastructures::MacroscopicCell<dim>*> getMacroscopicCells4Sending() {
    std::vector<coupling::datastructures::MacroscopicCell<dim>*> result;
    return result;
  }
  virtual unsigned int* getMacroscopicCellIndices4Sending() { return NULL; }

  virtual std::vector<coupling::datastructures::MacroscopicCell<dim>*> getMacroscopicCells4Receiving() {
    std::vector<coupling::datastructures::MacroscopicCell<dim>*> result;
    return result;
  }
  virtual unsigned int* getMacroscopicCellIndices4Receiving() { return NULL; }
};

#endif // _MOLECULARDYNAMICS_COUPLING_INTERFACE_VOIDMACROSCOPICSOLVERINTERFACE_H_
