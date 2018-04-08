// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_TRANSFERSTRATEGIES_DIRECTTRANSFERSTRATEGY_H_
#define _MOLECULARDYNAMICS_COUPLING_TRANSFERSTRATEGIES_DIRECTTRANSFERSTRATEGY_H_

#include "coupling/transferstrategies/TransferStrategy.h"
#include "coupling/cell-mappings/ComputeMassMapping.h"
#include "coupling/cell-mappings/ComputeMomentumMapping.h"
#include "tarch/la/Vector.h"

namespace coupling {
  namespace transferstrategies {
    template<class LinkedCell, unsigned int dim>
    class DirectTransferStrategy;
  }
}


/** transfers and introduces mass and momentum directly into MD and to macroscopic solver.
 *  So, for example, if mass M is coming from the macroscopic solver, M is to be inserted into MD.
 *
 *  @author Philipp Neumann
 */
template<class LinkedCell,unsigned int dim>
class coupling::transferstrategies::DirectTransferStrategy:
public coupling::transferstrategies::TransferStrategy<LinkedCell,dim> {
  public:
    DirectTransferStrategy(
      coupling::interface::MDSolverInterface<LinkedCell,dim> * const mdSolverInterface,
      const coupling::IndexConversion<dim> &indexConversion
    ): coupling::transferstrategies::TransferStrategy<LinkedCell,dim>(mdSolverInterface,indexConversion),
       _massMapping(mdSolverInterface), _momentumMapping(mdSolverInterface){}
    virtual ~DirectTransferStrategy(){}

    virtual void processInnerMacroscopicCellBeforeReceivingMacroscopicSolverData(
      coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell,dim> &cell, const unsigned int &index
    ){
      // reset quantities
      const tarch::la::Vector<dim,double> zero(0.0);
      cell.setMicroscopicMass(0.0);
      cell.setMicroscopicMomentum(zero);
    }

    virtual void processOuterMacroscopicCellBeforeReceivingMacroscopicSolverData(
      coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell,dim> &cell, const unsigned int &index
    ){
      // reset quantities
      const tarch::la::Vector<dim,double> zero(0.0);
      cell.setMicroscopicMass(0.0);
      cell.setMicroscopicMomentum(zero);
    }

    virtual void processInnerMacroscopicCellBeforeSendingMDSolverData(
      coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell,dim> &cell, const unsigned int &index
    ){
      cell.iterateConstCells(_massMapping);
      cell.iterateConstCells(_momentumMapping);
      cell.setMacroscopicMass(_massMapping.getMass());
      cell.setMacroscopicMomentum(_momentumMapping.getMomentum());
    }

  private:
    coupling::cellmappings::ComputeMassMapping<LinkedCell,dim> _massMapping;
    coupling::cellmappings::ComputeMomentumMapping<LinkedCell,dim> _momentumMapping;
};
#endif // _MOLECULARDYNAMICS_COUPLING_TRANSFERSTRATEGIES_DIRECTTRANSFERSTRATEGY_H_
