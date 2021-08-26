// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_TRANSFERSTRATEGIES_AVERAGINGTRANSFERSTRATEGY_H_
#define _MOLECULARDYNAMICS_COUPLING_TRANSFERSTRATEGIES_AVERAGINGTRANSFERSTRATEGY_H_

#include <map>
#include <list>
#include "coupling/transferstrategies/TransferStrategy.h"
#include "coupling/cell-mappings/ComputeMassMapping.h"
#include "coupling/cell-mappings/ComputeMomentumMapping.h"

namespace coupling {
  namespace transferstrategies {
    template<class LinkedCell,unsigned int dim>
    class AveragingTransferStrategy;
  }
}


/** this class is used for pure averaging operations on the macroscopic cells.
 *  This can be used e.g. to measure errors in averaging over time, to estimate number
 *  of samples etc.
 *  @author Philipp Neumann
 */
template<class LinkedCell,unsigned int dim>
class coupling::transferstrategies::AveragingTransferStrategy:
public coupling::transferstrategies::TransferStrategy<LinkedCell,dim> {
public:

  AveragingTransferStrategy(
    coupling::interface::MDSolverInterface<LinkedCell,dim> * const mdSolverInterface,
    const coupling::IndexConversion<dim> &indexConversion
  ): coupling::transferstrategies::TransferStrategy<LinkedCell,dim>(mdSolverInterface,indexConversion),
     _massMapping(mdSolverInterface),
     _momentumMapping(mdSolverInterface),
     _sampleCounter(0),
     _rank(indexConversion.getThisRank())
  {
  }
  virtual ~AveragingTransferStrategy(){}

  virtual void beginProcessInnerMacroscopicCellsBeforeReceivingMacroscopicSolverData(){
    // reset sample counter for each coupling cycle
    _sampleCounter = 0;
  }
  virtual void processInnerMacroscopicCellBeforeReceivingMacroscopicSolverData(
    coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell,dim> &cell, const unsigned int &index
  ){
    // reset buffers for sampling mass and momentum in each inner macroscopic cell
    cell.setMacroscopicMass(0.0);
    cell.setMacroscopicMomentum(tarch::la::Vector<dim,double>(0.0));
  }

  virtual void beginProcessInnerMacroscopicCellsAfterMDTimestep(){
    // output information of last sampling...
    if (_rank==0){
      std::cout << "Global quantities of sampling no. " << _sampleCounter << " on rank 0: mass=" << _avgMass << ", momentum=" << _avgMomentum << std::endl;
    }
    // reset avg. mass and momentum...
    _avgMass = 0.0;
    _avgMomentum = tarch::la::Vector<dim,double>(0.0);
    // and increment sample counter
    _sampleCounter++;
  }
  virtual void processInnerMacroscopicCellAfterMDTimestep(
    coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell,dim> &cell, const unsigned int &index
  ){
    // compute total mass/momentum from previous samples
    const double oldMass                            = (_sampleCounter-1)*cell.getMacroscopicMass();
    const tarch::la::Vector<dim,double> oldMomentum = ((double)(_sampleCounter-1))*cell.getMacroscopicMomentum();

    // compute new averaged mass and momentum
    cell.iterateConstCells(_massMapping);
    cell.iterateConstCells(_momentumMapping);
    const double mass                            = (1.0/_sampleCounter)*(oldMass+_massMapping.getMass());
    const tarch::la::Vector<dim,double> momentum = (1.0/_sampleCounter)*(oldMomentum+_momentumMapping.getMomentum());
    _avgMass += mass;
    _avgMomentum = _avgMomentum + momentum;
    // set mass and momentum in buffers
    cell.setMacroscopicMass(mass);
    cell.setMacroscopicMomentum(momentum);
  }

  private:
    // compute mass and momentum in a cell
    coupling::cellmappings::ComputeMassMapping<LinkedCell,dim> _massMapping;
    coupling::cellmappings::ComputeMomentumMapping<LinkedCell,dim> _momentumMapping;
    // required for counting the samples
    unsigned int _sampleCounter;

    // this rank -> only output info for rank 0
    const unsigned int _rank;
    // avg. velocity
    tarch::la::Vector<dim,double> _avgMomentum;
    // avg. mass
    double _avgMass;
};
#endif // _MOLECULARDYNAMICS_COUPLING_TRANSFERSTRATEGIES_AVERAGINGTRANSFERSTRATEGY_H_
