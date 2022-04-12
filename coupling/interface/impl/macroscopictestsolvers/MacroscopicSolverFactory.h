// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_INTERFACE_MACROSCOPICSOLVERFACTORY_H_
#define _MOLECULARDYNAMICS_COUPLING_INTERFACE_MACROSCOPICSOLVERFACTORY_H_

#include "coupling/interface/impl/macroscopictestsolvers/UniformFlowSolverInterface.h"
#include "coupling/interface/impl/macroscopictestsolvers/VoidMacroscopicSolver.h"
#include "coupling/interface/impl/macroscopictestsolvers/VoidMacroscopicSolverInterface.h"

namespace coupling {
namespace interface {
class MacroscopicSolverFactory;
}
} // namespace coupling

class coupling::interface::MacroscopicSolverFactory {
public:
  enum SolverType { VOID_MACROSCOPIC_SOLVER = 0, UNIFORM_FLOW_SOLVER = 1 };

  static coupling::interface::MacroscopicSolverFactory& getInstance() {
    static coupling::interface::MacroscopicSolverFactory singleton;
    return singleton;
  }

  coupling::interface::TestMacroscopicSolver* getTestMacroscopicSolver(SolverType type) {
    switch (type) {
    case VOID_MACROSCOPIC_SOLVER:
      return new coupling::interface::VoidMacroscopicSolver();
      break;
    case UNIFORM_FLOW_SOLVER:
      return new coupling::interface::VoidMacroscopicSolver();
      break;
    default:
      break;
    }

    return NULL;
  }

  template <unsigned int dim>
  coupling::interface::TestMacroscopicSolverInterface<dim>*
  getTestMacroscopicSolverInterface(SolverType type, tarch::la::Vector<dim, unsigned int> globalNumberMacroscopicCells, double massPerMacroscopicCell,
                                    tarch::la::Vector<dim, double> momentumPerMacroscopicCell) {
    switch (type) {
    case VOID_MACROSCOPIC_SOLVER:
      return new coupling::interface::VoidMacroscopicSolverInterface<dim>();
      break;
    case UNIFORM_FLOW_SOLVER:
      return new coupling::interface::UniformFlowSolverInterface<dim>(globalNumberMacroscopicCells, massPerMacroscopicCell, momentumPerMacroscopicCell);
    default:
      break;
    }

    return NULL;
  }

  private : MacroscopicSolverFactory() {
  }
  ~MacroscopicSolverFactory() {}
};

#endif // _MOLECULARDYNAMICS_COUPLING_INTERFACE_MACROSCOPICSOLVERFACTORY_H_
