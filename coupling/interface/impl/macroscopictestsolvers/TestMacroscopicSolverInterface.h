// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_INTERFACE_TESTMACROSCOPICSOLVERINTERFACE_H_
#define _MOLECULARDYNAMICS_COUPLING_INTERFACE_TESTMACROSCOPICSOLVERINTERFACE_H_

#include "coupling/interface/MacroscopicSolverInterface.h"
#include "coupling/datastructures/MacroscopicCell.h"


namespace coupling {
  namespace interface {
    template<unsigned int dim>
    class TestMacroscopicSolverInterface;
  }
}


template<unsigned int dim>
class coupling::interface::TestMacroscopicSolverInterface: public coupling::interface::MacroscopicSolverInterface<dim>{
  public:
    TestMacroscopicSolverInterface(): MacroscopicSolverInterface<dim>(){}
    virtual ~TestMacroscopicSolverInterface(){}

    virtual std::vector<coupling::datastructures::MacroscopicCell<dim>* > getMacroscopicCells4Sending() = 0;
    virtual unsigned int * getMacroscopicCellIndices4Sending() = 0;

    virtual std::vector<coupling::datastructures::MacroscopicCell<dim>* > getMacroscopicCells4Receiving() = 0;
    virtual unsigned int * getMacroscopicCellIndices4Receiving() = 0;
};
#endif // _MOLECULARDYNAMICS_COUPLING_INTERFACE_TESTMACROSCOPICSOLVERINTERFACE_H_
