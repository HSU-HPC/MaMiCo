// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_TRANSFERSTRATEGIES_TRANSFERSTRATEGY4SCHWARZCOUPLING_H_
#define _MOLECULARDYNAMICS_COUPLING_TRANSFERSTRATEGIES_TRANSFERSTRATEGY4SCHWARZCOUPLING_H_

#include "coupling/transferstrategies/TransferStrategy.h"
#include "coupling/cell-mappings/ComputeMassMapping.h"
#include "coupling/cell-mappings/ComputeMomentumMapping.h"

namespace coupling {
  namespace transferstrategies {
    template<class LinkedCell,unsigned int dim>
    class TransferStrategy4SchwarzCoupling;
  }
}


/** transfer strategy for Schwarz coupling algorithm, adopted from Dupuis et al.
 *  We currently sample over the last 20% of the coupling interval, i.e. of the numberMDsteps time steps in MD.
 *  The other 80% are used for equilibration.
 *  @author Philipp Neumann
 */
template<class LinkedCell,unsigned int dim>
class coupling::transferstrategies::TransferStrategy4SchwarzCoupling:
public coupling::transferstrategies::TransferStrategy<LinkedCell,dim> {
public:

  TransferStrategy4SchwarzCoupling(
    coupling::interface::MDSolverInterface<LinkedCell,dim> * const mdSolverInterface,
    const coupling::IndexConversion<dim> &indexConversion,
    unsigned int numberMDSteps
  ): coupling::transferstrategies::TransferStrategy<LinkedCell,dim>(mdSolverInterface,indexConversion),
     _massMapping(mdSolverInterface),
     _momentumMapping(mdSolverInterface),
     _timestepCounter(0),
     _sampleCounter(0),
     _numberMDSteps(numberMDSteps),
     _sampleEveryTimestep(1)
  {}
  virtual ~TransferStrategy4SchwarzCoupling(){}

  virtual void beginProcessInnerMacroscopicCellsBeforeReceivingMacroscopicSolverData();
  virtual void processInnerMacroscopicCellAfterReceivingMacroscopicSolverData(
    coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell,dim> &cell, const unsigned int &index
  );
  virtual void processOuterMacroscopicCellAfterReceivingMacroscopicSolverData(
    coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell,dim> &cell, const unsigned int &index
  );
  virtual void processInnerMacroscopicCellBeforeSendingMDSolverData(
    coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell,dim> &cell, const unsigned int &index
  );
  virtual void processOuterMacroscopicCellBeforeSendingMDSolverData(
    coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell,dim> &cell, const unsigned int &index
  );
  virtual void beginProcessInnerMacroscopicCellsAfterMDTimestep();
  virtual void processInnerMacroscopicCellAfterMDTimestep(
    coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell,dim> &cell, const unsigned int &index
  );

  protected:
    virtual bool sample() const;

  private:
    coupling::cellmappings::ComputeMassMapping<LinkedCell,dim> _massMapping;
    coupling::cellmappings::ComputeMomentumMapping<LinkedCell,dim> _momentumMapping;
    unsigned int _timestepCounter;
    unsigned int _sampleCounter;
    const unsigned int _numberMDSteps;
    const int _sampleEveryTimestep;
};
#include "coupling/transferstrategies/TransferStrategy4SchwarzCoupling.cpph"

#endif // _MOLECULARDYNAMICS_COUPLING_TRANSFERSTRATEGIES_TRANSFERSTRATEGY4SCHWARZCOUPLING_H_
