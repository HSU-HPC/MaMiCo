// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_INTERFACE_TESTMACROSCOPICSOLVERINTERFACE_H_
#define _MOLECULARDYNAMICS_COUPLING_INTERFACE_TESTMACROSCOPICSOLVERINTERFACE_H_

#include "coupling/datastructures/CouplingCell.h"
#include "coupling/interface/MacroscopicSolverInterface.h"

namespace coupling {
namespace interface {
template <unsigned int dim> class TestMacroscopicSolverInterface;
}
} // namespace coupling

template <unsigned int dim> class coupling::interface::TestMacroscopicSolverInterface : public coupling::interface::MacroscopicSolverInterface<dim> {
public:
  TestMacroscopicSolverInterface() : MacroscopicSolverInterface<dim>() {}
  virtual ~TestMacroscopicSolverInterface() {}

  virtual std::vector<coupling::datastructures::CouplingCell<dim>*> getCouplingCells4Sending() = 0;
  virtual unsigned int* getCouplingCellIndices4Sending() = 0;

  virtual std::vector<coupling::datastructures::CouplingCell<dim>*> getCouplingCells4Receiving() = 0;
  virtual unsigned int* getCouplingCellIndices4Receiving() = 0;
};
#endif // _MOLECULARDYNAMICS_COUPLING_INTERFACE_TESTMACROSCOPICSOLVERINTERFACE_H_
