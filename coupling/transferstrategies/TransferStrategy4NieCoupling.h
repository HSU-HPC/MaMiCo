// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_TRANSFERSTRATEGIES_TRANSFERSTRATEGY4NIECOUPLING_H_
#define _MOLECULARDYNAMICS_COUPLING_TRANSFERSTRATEGIES_TRANSFERSTRATEGY4NIECOUPLING_H_

#include "coupling/transferstrategies/TransferStrategy.h"
#include "coupling/cell-mappings/ComputeMassMapping.h"
#include "coupling/cell-mappings/ComputeMomentumMapping.h"

namespace coupling {
  namespace transferstrategies {
    template<class LinkedCell,unsigned int dim>
    class TransferStrategy4NieCoupling;
  }
}


/** transfer strategy for Nie coupling algorithm, adopted from:
 *  X.B. Nie, S.Y. Chen, W.N. E, M.O. Robbins
 *  A continuum and molecular dynamics hybrid method for micro- and nano-fluid flow
 *  J. Fluid. Mech. 500: 55-64, 2004
 *  We basically sample in every MD time step.
 *  @author Philipp Neumann
 */
template<class LinkedCell,unsigned int dim>
class coupling::transferstrategies::TransferStrategy4NieCoupling: public coupling::transferstrategies::TransferStrategy<LinkedCell,dim> {
  public:

    TransferStrategy4NieCoupling(
      coupling::interface::MDSolverInterface<LinkedCell,dim> * const mdSolverInterface,
      const coupling::IndexConversion<dim> &indexConversion,
      unsigned int numberMDSteps,double shiftTimestep,
      tarch::la::Vector<2*dim,bool> massFluxBoundary
    );
    virtual ~TransferStrategy4NieCoupling();

    /** stores the old cont.-velocity field solution and resets time step counter */
    virtual void beginProcessInnerMacroscopicCellsBeforeReceivingMacroscopicSolverData();
    /** store old microscopic mass in excess-mass buffer and reset microscopic mass buffer */
    virtual void processInnerMacroscopicCellBeforeReceivingMacroscopicSolverData(
      coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell,dim> &cell, const unsigned int &index
    );
    /** converts momentum into velocity values; stores velocity values in new cont.-velocity field solution and sets correct velocity value for first MD time step in
     *  microscopic momentum buffer. */
    virtual void processInnerMacroscopicCellAfterReceivingMacroscopicSolverData(
      coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell,dim> &cell, const unsigned int &index
    );
    virtual void endProcessInnerMacroscopicCellsAfterReceivingMacroscopicSolverData()override;

    virtual void beginProcessInnerMacroscopicCellsBeforeSendingMDSolverData();
    /** divides accumulated mass and momentum values by time step counter.
     */
    virtual void processInnerMacroscopicCellBeforeSendingMDSolverData(
      coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell,dim> &cell, const unsigned int &index
    );
    virtual void endProcessInnerMacroscopicCellsBeforeSendingMDSolverData();
    /** increments time step counter */
    virtual void beginProcessInnerMacroscopicCellsAfterMDTimestep();
    /** computes current velocity (linear time interpolation) in this cell and accumulates mass/momentum for sampling */
    virtual void processInnerMacroscopicCellAfterMDTimestep(
      coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell,dim> &cell, const unsigned int &index
    );

  private:
    /** returns the local number of macroscopic cells incl. ghost layers */
    unsigned int getLocalNumberMacroscopicCells(const coupling::IndexConversion<dim> &indexConversion) const;

    /** computes the mass flux in the outermost inner macroscopic cells. For all other cells, 0.0 is returned. */
    double computeMassFlux(const double& mass, const tarch::la::Vector<dim,double>& velocity, const unsigned int index) ;

    coupling::cellmappings::ComputeMassMapping<LinkedCell,dim> _massMapping;
    coupling::cellmappings::ComputeMomentumMapping<LinkedCell,dim> _momentumMapping;
    tarch::la::Vector<dim,double> *_oldSolution; // solution of velocity field at the beginning of coupling cycle (=end of old coupling cycle)
    tarch::la::Vector<dim,double> *_newSolution; // solution of velocity field at the end of coupling cycle
    const unsigned int _numberMDSteps;    // number of MD time steps per coupling cycle
    const double _shiftTimestep;          // offset by which evaluation of time interval is shifted compared to continuum solver. A value of zero implies interpolation between new and old solution,
                                          // corresponding to a time interval t->t+dt.
                                          // A value of one implies extrapolation between the times t+dt (newSolution) and t+2dt using oldSolution and newSolution. All inter-/extrapolations are linear.
    const unsigned int _numberLocalCells; // local number of cells (incl. ghost layers)
    unsigned int _timestepCounter;        // time step counter within a coupling cycle (should run from 0 to _numberMDSteps)
    double *_excessMass;                  // mass that was transferred in an earlier coupled step
    const tarch::la::Vector<2*dim,bool> _massFluxBoundary; // true in each entry if west/east, south/north, bottom/top boundary is a mass flux boundary
    //double _totalMass; int _cellCount;
};
#include "coupling/transferstrategies/TransferStrategy4NieCoupling.cpph"

#endif // _MOLECULARDYNAMICS_COUPLING_TRANSFERSTRATEGIES_TRANSFERSTRATEGY4NIECOUPLING_H_
