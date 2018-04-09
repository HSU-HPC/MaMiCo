// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_TRANSFERSTRATEGIES_STRATEGY_H_
#define _MOLECULARDYNAMICS_COUPLING_TRANSFERSTRATEGIES_STRATEGY_H_

#include "coupling/datastructures/MacroscopicCell.h"
#include "coupling/interface/MDSolverInterface.h"
#include "coupling/IndexConversion.h"

namespace coupling {
  namespace transferstrategies {
    template<class LinkedCell, unsigned int dim>
    class TransferStrategy;
  }
}


/** interface for transfer strategy, that is for macroscopic cell operations that are carried out before/after send/recv-operations between
 *  the two solvers, or in each MD step for sampling purposes.
 *  @author Philipp Neumann
 */
template<class LinkedCell,unsigned int dim>
class coupling::transferstrategies::TransferStrategy {
  public:
    TransferStrategy(
      coupling::interface::MDSolverInterface<LinkedCell,dim> * const mdSolverInterface,
      const coupling::IndexConversion<dim> &indexConversion
    ): _mdSolverInterface(mdSolverInterface), _indexConversion(indexConversion){}
    virtual ~TransferStrategy(){}

    /** is called before macroscopic cell information is transferred from the macroscopic solver to the macroscopic
     *  cell. Here, you might e.g. reset all macroscopic cell values to zero.
     *  Afterwards, all values from the send/receive-buffer are added to the mass/momentum buffers of the macroscopic
     *  cells.
     *  This method is only applied to macroscopic cells that cover parts of the MD domain; it is not applied
     *  in the outer macroscopic cells.
     *  The methods begin...() and end..() are called write before/after the traversal of the local cells.
     */
    virtual void processInnerMacroscopicCellBeforeReceivingMacroscopicSolverData(
      coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell,dim> &cell, const unsigned int &index
    ){}
    virtual void beginProcessInnerMacroscopicCellsBeforeReceivingMacroscopicSolverData(){}
    virtual void endProcessInnerMacroscopicCellsBeforeReceivingMacroscopicSolverData(){}

    /** is called before macroscopic cell information is transferred from the macroscopic solver to the macroscopic
     *  cell. Here, you might e.g. reset all macroscopic cell values to zero.
     *  Afterwards, all values from the send/receive-buffer are added to the mass/momentum buffers of the macroscopic
     *  cells.
     *  This method is only applied to outer macroscopic cells, that is cells that are located outside the MD domain.
     */
    virtual void processOuterMacroscopicCellBeforeReceivingMacroscopicSolverData(
      coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell,dim> &cell, const unsigned int &index
    ){}
    virtual void beginProcessOuterMacroscopicCellsBeforeReceivingMacroscopicSolverData(){}
    virtual void endProcessOuterMacroscopicCellsBeforeReceivingMacroscopicSolverData(){}

    /** is called for every macroscopic cell right after the microscopicMass and -Momentum-values
     *  have been filled in with information from the macroscopic solver.
     *  These values might be manipulated within this function. Example: Mass and momentum have been
     *  transferred from the macroscopic solver to MD, but only the difference between MD and macroscopic
     *  solver shall be introduced to MD. Then, this function determines this difference and stores
     *  the result again in microscopicMass and -Momentum.
     *  This method is only applied to macroscopic cells that cover parts of the MD domain; it is not applied
     *  in the outer macroscopic cells.
     */
    virtual void processInnerMacroscopicCellAfterReceivingMacroscopicSolverData(
      coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell,dim> &cell, const unsigned int &index
    ){}
    virtual void beginProcessInnerMacroscopicCellsAfterReceivingMacroscopicSolverData(){}
    virtual void endProcessInnerMacroscopicCellsAfterReceivingMacroscopicSolverData(){}
    virtual void processOuterMacroscopicCellAfterReceivingMacroscopicSolverData(
      coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell,dim> &cell, const unsigned int &index
    ){}
    virtual void beginProcessOuterMacroscopicCellsAfterReceivingMacroscopicSolverData(){}
    virtual void endProcessOuterMacroscopicCellsAfterReceivingMacroscopicSolverData(){}

    /** is called for every macroscopic cell before sending the macroscopicMass and -Momentum
     *  data to the macroscopic solver and before noise reduction invocation. 
     *  Example: Compute mass and momentum and store the results in
     *  macroscopicMass and -Momentum. The total mass and momentum from the MD system will then be
     *  sent to the macroscopic solver.
     *  This method is only applied to macroscopic cells that cover parts of the MD domain; it is not applied
     *  in the outer macroscopic cells.
     */
    virtual void processInnerMacroscopicCellBeforeSendingMDSolverData(
      coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell,dim> &cell, const unsigned int &index
    ){}
    virtual void beginProcessInnerMacroscopicCellsBeforeSendingMDSolverData(){}
    virtual void endProcessInnerMacroscopicCellsBeforeSendingMDSolverData(){}

    /** is called for every macroscopic cell before sending the macroscopicMass and -Momentum
     *  data to the macroscopic solver and before noise reduction invocation. 
     *  Example: Compute mass and momentum and store the results in
     *  macroscopicMass and -Momentum. The total mass and momentum from the MD system will then be
     *  sent to the macroscopic solver.
     *  This method is only applied to outer macroscopic cells, that is cells that are located outside the MD domain.
     */
    virtual void processOuterMacroscopicCellBeforeSendingMDSolverData(
      coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell,dim> &cell, const unsigned int &index
    ){}
    virtual void beginProcessOuterMacroscopicCellsBeforeSendingMDSolverData(){}
    virtual void endProcessOuterMacroscopicCellsBeforeSendingMDSolverData(){}

    /** required to collect cell data during an MD simulation. For example, if we need time-averaged
     *  data within a macroscopic cell, we can compute mass and momentum in each timestep and add
     *  it to the microscopicMass and -Momentum buffers.
     *  This method is only applied to macroscopic cells that cover parts of the MD domain; it is not applied
     *  in the outer macroscopic cells.
     */
    virtual void processInnerMacroscopicCellAfterMDTimestep(
      coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell,dim> &cell, const unsigned int &index
    ){}
    virtual void beginProcessInnerMacroscopicCellsAfterMDTimestep(){}
    virtual void endProcessInnerMacroscopicCellsAfterMDTimestep(){}

  protected:
    coupling::interface::MDSolverInterface<LinkedCell,dim> * const _mdSolverInterface;
    const coupling::IndexConversion<dim> &_indexConversion;
};
#endif // _MOLECULARDYNAMICS_COUPLING_TRANSFERSTRATEGIES_STRATEGY_H_
