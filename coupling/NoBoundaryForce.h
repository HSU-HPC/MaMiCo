// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_NOBOUNDARYFORCE_H_
#define _MOLECULARDYNAMICS_COUPLING_NOBOUNDARYFORCE_H_

#include "coupling/BoundaryForceController.h"

namespace coupling {
  template<class LinkedCell,unsigned int dim>
  class NoBoundaryForce;
}

/** dummy implementation, applying no boundary force.
 *  @author Philipp Neumann
 */
template<class LinkedCell,unsigned int dim>
class coupling::NoBoundaryForce: public coupling::BoundaryForceController<LinkedCell,dim> {
  public:
    NoBoundaryForce(
      coupling::interface::MDSolverInterface<LinkedCell,dim>* mdSolverInterface
    ): coupling::BoundaryForceController<LinkedCell,dim>(mdSolverInterface){}
    virtual ~NoBoundaryForce(){}

    virtual void applyBoundaryForce(coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell,dim>& cell, const unsigned int &currentLocalMacroscopicCellIndex){}
};
#endif // _MOLECULARDYNAMICS_COUPLING_NOBOUNDARYFORCE_H_

