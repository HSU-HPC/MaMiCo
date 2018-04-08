// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_ZHOUBOUNDARYFORCECONTROLLER_H_
#define _MOLECULARDYNAMICS_COUPLING_ZHOUBOUNDARYFORCECONTROLLER_H_

#include "coupling/cell-mappings/ZhouBoundaryForce.h"
#include "coupling/BoundaryForceController.h"

namespace coupling {
  template<class LinkedCell,unsigned int dim>
  class ZhouBoundaryForceController;
}


/** applies the interpolating boundary force from Zhou et al. in every boundary cell. For details on the forcing,
 *  check out the descriptions in cell-mappings/ZhouBoundaryForce.
 *  @author Philipp Neumann
 */
template<class LinkedCell,unsigned int dim>
class coupling::ZhouBoundaryForceController: public coupling::BoundaryForceController<LinkedCell,dim> {
  public:
    ZhouBoundaryForceController(
      const double& density,
      const double& temperature,
      const tarch::la::Vector<2*dim,bool>& boundary,
      coupling::interface::MDSolverInterface<LinkedCell,dim> * const mdSolverInterface
    ): coupling::BoundaryForceController<LinkedCell,dim>(mdSolverInterface), _density(density), _temperature(temperature),_boundary(boundary){}
    virtual ~ZhouBoundaryForceController(){}

    virtual void applyBoundaryForce(coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell,dim>& cell, const unsigned int &currentLocalMacroscopicCellIndex){
      const tarch::la::Vector<dim,double> domainOffset(coupling::BoundaryForceController<LinkedCell,dim>::_mdSolverInterface->getGlobalMDDomainOffset());
      const tarch::la::Vector<dim,double> domainSize(coupling::BoundaryForceController<LinkedCell,dim>::_mdSolverInterface->getGlobalMDDomainSize());
      const double epsilon(coupling::BoundaryForceController<LinkedCell,dim>::_mdSolverInterface->getMoleculeEpsilon());
      const double sigma(coupling::BoundaryForceController<LinkedCell,dim>::_mdSolverInterface->getMoleculeSigma());
      coupling::cellmappings::ZhouBoundaryForce<LinkedCell,dim> zhouBoundaryForce(
        _density,_temperature,epsilon,sigma,_boundary,domainOffset,domainSize,
        coupling::BoundaryForceController<LinkedCell,dim>::_mdSolverInterface
      );
      cell.iterateCells(zhouBoundaryForce);
    }

  private:
    const double _density;
    const double _temperature;
    const tarch::la::Vector<2*dim,bool> _boundary;
};
#endif // _MOLECULARDYNAMICS_COUPLING_ZHOUBOUNDARYFORCECONTROLLER_H_

