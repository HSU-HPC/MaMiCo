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

  /** no cells are excluded from the md2macro region; this is a no-op dummy solver.
   *  Replaces the legacy receive-/sendMacroscopicQuantity*MDSolver() pair, cf. MacroscopicSolverInterface.
   */
  unsigned int getOuterRegion() override { return 0; }

  /** returns the ranks on which the macroscopic solver holds/requires data of the coupling cell
   *  at index idx. All cells are located on rank 0.
   */
  std::vector<unsigned int> getRanks(I01 idx) override {
    std::vector<unsigned int> result;
    result.push_back(0);
    return result;
  }

  virtual std::vector<coupling::datastructures::CouplingCell<dim>*> getCouplingCells4Sending() {
    std::vector<coupling::datastructures::CouplingCell<dim>*> result;
    return result;
  }
  virtual unsigned int* getCouplingCellIndices4Sending() { return NULL; }

  virtual std::vector<coupling::datastructures::CouplingCell<dim>*> getCouplingCells4Receiving() {
    std::vector<coupling::datastructures::CouplingCell<dim>*> result;
    return result;
  }
  virtual unsigned int* getCouplingCellIndices4Receiving() { return NULL; }
};

#endif // _MOLECULARDYNAMICS_COUPLING_INTERFACE_VOIDMACROSCOPICSOLVERINTERFACE_H_
