// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_TRANSFERSTRATEGIES_STRATEGY_H_
#define _MOLECULARDYNAMICS_COUPLING_TRANSFERSTRATEGIES_STRATEGY_H_

#include "coupling/IndexConversion.h"
#include "coupling/datastructures/CouplingCell.h"
#include "coupling/interface/MDSolverInterface.h"

namespace coupling {
namespace transferstrategies {
template <class LinkedCell, unsigned int dim> class TransferStrategy;
}
} // namespace coupling

/** interface for transfer strategy, that is for macroscopic cell operations
 * that are carried out before/after send/recv-operations between the two
 * solvers, or in each MD step for sampling purposes.
 *  @author Philipp Neumann
 *  @tparam LinkedCell the LinkedCell class is given by the implementation of
 * linked cells in the molecular dynamics simulation
 *  @tparam dim refers to the spacial dimension of the simulation, can be 1, 2,
 * or 3  */
template <class LinkedCell, unsigned int dim> class coupling::transferstrategies::TransferStrategy {
public:
  /** @brief a simple destructor
   *  @param mdSolverInterface interface to md solver
   *  @param indexConversion an instance of the indexConversion for the current
   * simulation  */
  TransferStrategy(coupling::interface::MDSolverInterface<LinkedCell, dim>* const mdSolverInterface, const coupling::IndexConversion<dim>& indexConversion)
      : _mdSolverInterface(mdSolverInterface), _indexConversion(indexConversion) {}

  /** @brief a dummy destructor */
  virtual ~TransferStrategy() {}

  /** Here, you might e.g. reset all macroscopic cell values to zero.
   *  This method is only applied to macroscopic cells that cover parts of the
   * MD domain; it is not applied in the outer macroscopic cells.
   *  @brief is called on the inner macroscopic cells before the data from the
   * macro solver is applied
   *  @param cell the macroscopic cell to apply
   *  @param index the index of the macroscopic cell */
  virtual void processInnerCouplingCellBeforeReceivingMacroscopicSolverData(coupling::datastructures::CouplingCellWithLinkedCells<LinkedCell, dim>& cell,
                                                                               const unsigned int& index) {}

  /** This method is only applied to macroscopic cells that cover parts of the
   * MD domain; it is not applied in the outer macroscopic cells. Is called
   * before the
   * processInnerCouplingCellBeforeReceivingMacroscopicSolverData() method
   *  @brief the method is called before the cells are processed, e.g. general
   * values might be set here */
  virtual void beginProcessInnerCouplingCellsBeforeReceivingMacroscopicSolverData() {}

  /** This method is only applied to macroscopic cells that cover parts of the
   * MD domain; it is not applied in the outer macroscopic cells. Is called
   * after the processInnerCouplingCellBeforeReceivingMacroscopicSolverData()
   * method
   *  @brief the method is called after the cells are processed, e.g. some
   * general evaluation might happen like sum/counter */
  virtual void endProcessInnerCouplingCellsBeforeReceivingMacroscopicSolverData() {}

  /** Here, you might e.g. reset all macroscopic cell values to zero.
   *  This method is only applied to outer macroscopic cells that are not part
   * of the inner MD domain; it is not applied to the cells in the ghost layer
   *  @brief is called on the outer macroscopic cells before the data from the
   * macro solver is applied
   *  @param cell the macroscopic cell to apply
   *  @param index the index of the macroscopic cell */
  virtual void processOuterCouplingCellBeforeReceivingMacroscopicSolverData(coupling::datastructures::CouplingCellWithLinkedCells<LinkedCell, dim>& cell,
                                                                               const unsigned int& index) {}

  /** This method is only applied to macroscopic cells that are not part the
   * inner MD domain; it is not applied to the ghost layer Is called before the
   * processOuterCouplingCellBeforeReceivingMacroscopicSolverData() method
   *  @brief the method is called before the outer cells are processed, e.g.
   * general values might be set here */
  virtual void beginProcessOuterCouplingCellsBeforeReceivingMacroscopicSolverData() {}

  /** This method is only applied to macroscopic cells that are not part the
   * inner MD domain; it is not applied to the ghost layer Is called after the
   * processOuterCouplingCellBeforeReceivingMacroscopicSolverData() method
   *  @brief the method is called after the outer cells are processed, e.g.
   * general values might be set here */
  virtual void endProcessOuterCouplingCellsBeforeReceivingMacroscopicSolverData() {}

  /** These values might be manipulated within this function. Example: Mass and
   * momentum have been transferred from the macroscopic solver to MD, but only
   * the difference between MD and macroscopic solver shall be introduced to MD.
   * Then, this function determines this difference and stores the result again
   * in microscopicMass and -Momentum. This method is only applied to
   * macroscopic cells that cover parts of the MD domain; it is not applied in
   * the outer macroscopic cells.
   *  @brieff is called for every macroscopic cell after the microscopicMass and
   * -Momentum have been filled in with data from the macroscopic solver.
   *  @param cell the macroscopic cell to apply
   *  @param index the index of the macroscopic cell */
  virtual void processInnerCouplingCellAfterReceivingMacroscopicSolverData(coupling::datastructures::CouplingCellWithLinkedCells<LinkedCell, dim>& cell,
                                                                              const unsigned int& index) {}

  /** This method is only applied to macroscopic cells that cover parts of the
   * MD domain; it is not applied the outer macroscopic cells. Is called before
   * the processInnerCouplingCellAfterReceivingMacroscopicSolverData() method
   *  @brief the method is called before the inner cells are processed, e.g.
   * general values might be set here */
  virtual void beginProcessInnerCouplingCellsAfterReceivingMacroscopicSolverData() {}

  /** This method is only applied to macroscopic cells that cover parts of the
   * MD domain; it is not applied the outer macroscopic cells. Is called after
   * the processInnerCouplingCellAfterReceivingMacroscopicSolverData() method
   *  @brief the method is called after the inner cells are processed, e.g.
   * general values might be set here */
  virtual void endProcessInnerCouplingCellsAfterReceivingMacroscopicSolverData() {}

  /** These values might be manipulated within this function. Example: Mass and
   * momentum have been transferred from the macroscopic solver to MD, but only
   * the difference between MD and macroscopic solver shall be introduced to MD.
   * Then, this function determines this difference and stores the result again
   * in microscopicMass and -Momentum. This method is only applied to
   * macroscopic cells in the outer parts of the MD domain; it is not applied in
   * the ghost layer.
   *  @brief is called for every macroscopic cell after the microscopicMass and
   * -Momentum have been filled in with data from the macroscopic solver.
   *  @param cell the macroscopic cell to apply
   *  @param index the index of the macroscopic cell */
  virtual void processOuterCouplingCellAfterReceivingMacroscopicSolverData(coupling::datastructures::CouplingCellWithLinkedCells<LinkedCell, dim>& cell,
                                                                              const unsigned int& index) {}

  /** This method is only applied to macroscopic cells are part of the outer MD
   * domain; it is not applied the inner macroscopic cells. Is called before the
   * processOuterCouplingCellAfterReceivingMacroscopicSolverData() method
   *  @brief the method is called before the outer cells are processed, e.g.
   * general values might be set here */
  virtual void beginProcessOuterCouplingCellsAfterReceivingMacroscopicSolverData() {}

  /** This method is only applied to macroscopic cells are part of the outer MD
   * domain; it is not applied the inner macroscopic cells. Is called after the
   * processOuterCouplingCellAfterReceivingMacroscopicSolverData() method
   *  @brief the method is called after the outer cells are processed, e.g.
   * general values might be set here */
  virtual void endProcessOuterCouplingCellsAfterReceivingMacroscopicSolverData() {}

  /** Example: Compute mass and momentum and store the results in
   *  macroscopicMass and -Momentum. The total mass and momentum from the MD
   * system will then be sent to the macroscopic solver. This method is only
   * applied to macroscopic cells that cover parts of the MD domain; it is not
   * applied in the outer macroscopic cells.
   *  @brief is called for every macroscopic cell before sending the
   * macroscopicMass and -Momentum data to the macroscopic solver and before
   * noise reduction invocation.
   *  @param cell the macroscopic cell to apply
   *  @param index the index of the macroscopic cell */
  virtual void processInnerCouplingCellBeforeSendingMDSolverData(coupling::datastructures::CouplingCellWithLinkedCells<LinkedCell, dim>& cell,
                                                                    const unsigned int& index) {}

  /** This method is only applied to macroscopic cells that cover parts of the
   * MD domain; it is not applied the outer macroscopic cells. Is called before
   * the processInnerCouplingCellBeforeSendingMDSolverData() method
   *  @brief the method is called before the inner cell data is send, e.g.
   * general values might be set here */
  virtual void beginProcessInnerCouplingCellsBeforeSendingMDSolverData() {}

  /** This method is only applied to macroscopic cells that cover parts of the
   * MD domain; it is not applied the outer macroscopic cells. Is called after
   * the processInnerCouplingCellBeforeSendingMDSolverData() method
   *  @brief the method is called after the inner cell data is send, e.g.
   * general values might be set here */
  virtual void endProcessInnerCouplingCellsBeforeSendingMDSolverData() {}

  /** is called for every macroscopic cell before sending the macroscopicMass
   * and -Momentum data to the macroscopic solver and before noise reduction
   * invocation. Example: Compute mass and momentum and store the results in
   *  macroscopicMass and -Momentum. The total mass and momentum from the MD
   * system will then be sent to the macroscopic solver. This method is only
   * applied to outer macroscopic cells, that is cells that are located outside
   * the MD domain.
   *  @brief is called for outer macroscopic cell before sending the
   * macroscopicMass and -Momentum data to the macroscopic solver and before
   * noise reduction invocation.
   *  @param cell the macroscopic cell to apply
   *  @param index the index of the macroscopic cell */
  virtual void processOuterCouplingCellBeforeSendingMDSolverData(coupling::datastructures::CouplingCellWithLinkedCells<LinkedCell, dim>& cell,
                                                                    const unsigned int& index) {}

  /** This method is only applied to macroscopic cells are part of the outer MD
   * domain; it is not applied the inner macroscopic cells. Is called before the
   * processOuterCouplingCellBeforeSendingMDSolverData() method
   *  @brief the method is called before the outer cells are processed, e.g.
   * general values might be set here */
  virtual void beginProcessOuterCouplingCellsBeforeSendingMDSolverData() {}

  /** This method is only applied to macroscopic cells are part of the outer MD
   * domain; it is not applied the inner macroscopic cells. Is called after the
   * processOuterCouplingCellBeforeSendingMDSolverData() method
   *  @brief the method is called after the outer cells are processed, e.g.
   * general values might be set here */
  virtual void endProcessOuterCouplingCellsBeforeSendingMDSolverData() {}

  /** required to collect cell data during an MD simulation. For example, if we
   * need time-averaged data within a macroscopic cell, we can compute mass and
   * momentum in each timestep and add it to the microscopicMass and -Momentum
   * buffers. This method is only applied to macroscopic cells that cover parts
   * of the MD domain; it is not applied in the outer macroscopic cells.
   *  @param cell the macroscopic cell to apply
   *  @param index the index of the macroscopic cell */
  virtual void processInnerCouplingCellAfterMDTimestep(coupling::datastructures::CouplingCellWithLinkedCells<LinkedCell, dim>& cell,
                                                          const unsigned int& index) {}

  /** This method is only applied to macroscopic cells that cover parts of the
   * MD domain; it is not applied the outer macroscopic cells. Is called before
   * the processInnerCouplingCellAfterMDTimestep() method
   *  @brief the method is called before the inner cell data is send, e.g.
   * general values might be set here */
  virtual void beginProcessInnerCouplingCellsAfterMDTimestep() {}

  /** This method is only applied to macroscopic cells that cover parts of the
   * MD domain; it is not applied the outer macroscopic cells. Is called after
   * the processInnerCouplingCellAfterMDTimestep() method
   *  @brief the method is called after the inner cell data is send, e.g.
   * general values might be set here */
  virtual void endProcessInnerCouplingCellsAfterMDTimestep() {}

protected:
  /** interface for the md solver */
  coupling::interface::MDSolverInterface<LinkedCell, dim>* const _mdSolverInterface;
  /** an instance of the indexConversion current setup */
  const coupling::IndexConversion<dim>& _indexConversion;
};
#endif // _MOLECULARDYNAMICS_COUPLING_TRANSFERSTRATEGIES_STRATEGY_H_
