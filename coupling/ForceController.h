// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_FORCECONTROLLER_H_
#define _MOLECULARDYNAMICS_COUPLING_FORCECONTROLLER_H_

#include "coupling/interface/MDSolverInterface.h"
#include "coupling/datastructures/MacroscopicCell.h"
#include "coupling/cell-mappings/LimitForce.h"

namespace coupling {
  template<class LinkedCell,unsigned int dim>
  class ForceController;
}


/** controller for forces acting at open MD boundaries. There is an interface method
 *  applyBoundaryForce which triggers potential boundary forcing in each macroscopic
 *  cell that is located at the very outer MD boundary (first layer of non-ghost macroscopic cells).
 *  @author Philipp Neumann
 */
template<class LinkedCell,unsigned int dim>
class coupling::ForceController {
  public:
    ForceController(
      coupling::interface::MDSolverInterface<LinkedCell,dim> * const mdSolverInterface
    ): _mdSolverInterface(mdSolverInterface),_limitForce(coupling::cellmappings::LimitForce(mdSolverInterface)){}
    virtual ~ForceController(){}

    virtual void limitForce(coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell,dim>& cell, const unsigned int &currentLocalMacroscopicCellIndex){
      cell.iterateCells(_limitForce);
    };

  protected:
    coupling::interface::MDSolverInterface<LinkedCell,dim> * const _mdSolverInterface;
    coupling::cellmappings::LimitForce<LinkedCell,dim> _limitForce;
};
#endif // _MOLECULARDYNAMICS_COUPLING_FORCECONTROLLER_H_
