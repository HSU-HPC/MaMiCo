// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_TRANSFERSTRATEGIES_TRANSFERSTRATEGY4NIECOUPLING_H_
#define _MOLECULARDYNAMICS_COUPLING_TRANSFERSTRATEGIES_TRANSFERSTRATEGY4NIECOUPLING_H_

#include "coupling/cell-mappings/ComputeMassMapping.h"
#include "coupling/cell-mappings/ComputeMomentumMapping.h"
#include "coupling/transferstrategies/TransferStrategy.h"

namespace coupling {
namespace transferstrategies {
template <class LinkedCell, unsigned int dim> class TransferStrategy4NieCoupling;
}
} // namespace coupling

/** transfer strategy for Nie coupling algorithm, adopted from:
 *  X.B. Nie, S.Y. Chen, W.N. E, M.O. Robbins
 *  A continuum and molecular dynamics hybrid method for micro- and nano-fluid
 * flow J. Fluid. Mech. 500: 55-64, 2004 We basically sample in every MD time
 * step.
 *  @author Philipp Neumann
 *  @tparam LinkedCell the LinkedCell class is given by the implementation of
 * linked cells in the molecular dynamics simulation
 *  @tparam dim  refers to the spacial dimension of the simulation, can be 1, 2,
 * or 3 */
template <class LinkedCell, unsigned int dim>
class coupling::transferstrategies::TransferStrategy4NieCoupling : public coupling::transferstrategies::TransferStrategy<LinkedCell, dim> {
public:
  /** @brief a simple
   *  @param mdSolverInterface interface for the md solver
   *  @param shiftTimestep parameter for the inter- and extrapolation of the
   * time steps refers to the Nie algorithm
   *  @param massFluxBoundary indicates at which boundary mass shall be
   * transferred (true) or not (false); in their order the entries refer to:
   * west, east, south, north, bottom, top boundary  */
  TransferStrategy4NieCoupling(coupling::interface::MDSolverInterface<LinkedCell, dim>* const mdSolverInterface, unsigned int numberMDSteps,
                               double shiftTimestep, tarch::la::Vector<2 * dim, bool> massFluxBoundary);

  /** @brief a dummy destructor */
  virtual ~TransferStrategy4NieCoupling();

  /** @brief stores the old cont.-velocity field solution and resets time step
   * counter */
  void beginProcessInnerCouplingCellsBeforeReceivingMacroscopicSolverData() override;

  /** @brief store old microscopic mass in excess-mass buffer and reset
   * microscopic mass buffer
   *  @param cell the coupling cell to be processed
   *  @param index the index of the coupling cell */
  void processInnerCouplingCellBeforeReceivingMacroscopicSolverData(coupling::datastructures::CouplingCellWithLinkedCells<LinkedCell, dim>& cell,
                                                                    I02 index) override;

  /** stores velocity values in new cont.-velocity field solution and sets
   * correct velocity value for first MD time step in microscopic momentum
   * buffer.
   *  @brief converts momentum into velocity values;
   *  @param cell the coupling cell to be processed
   *  @param index the index of the coupling cell */
  void processInnerCouplingCellAfterReceivingMacroscopicSolverData(coupling::datastructures::CouplingCellWithLinkedCells<LinkedCell, dim>& cell,
                                                                   I02 index) override;

  // virtual void beginProcessInnerCouplingCellsBeforeSendingMDSolverData();

  /** @brief divides accumulated mass and momentum values by time step counter.
   */
  void processInnerCouplingCellBeforeSendingMDSolverData(coupling::datastructures::CouplingCellWithLinkedCells<LinkedCell, dim>& cell, I02 index) override;

  // virtual void endProcessInnerCouplingCellsBeforeSendingMDSolverData();

  /** @brief increments time step counter */
  void beginProcessInnerCouplingCellsAfterMDTimestep() override;

  /** @brief computes current velocity (linear time interpolation) in this cell
   * and accumulates mass/momentum for sampling
   *  @param cell the coupling cell to be processed
   *  @param index the index of the coupling cell */
  void processInnerCouplingCellAfterMDTimestep(coupling::datastructures::CouplingCellWithLinkedCells<LinkedCell, dim>& cell, I02 index) override;

private:
  /** computes the mass flux in the outermost inner coupling cells. For all
   * other cells, 0.0 is returned. */
  double computeMassFlux(const double& mass, const tarch::la::Vector<dim, double>& velocity, I01 index);
  /** class to compute the amount of mass in every single cell*/
  coupling::cellmappings::ComputeMassMapping<LinkedCell, dim> _massMapping;
  /** class to compute the amount of momentum in every single cell*/
  coupling::cellmappings::ComputeMomentumMapping<LinkedCell, dim> _momentumMapping;
  /** solution of velocity field at the beginning of coupling cycle (=end of old
   * coupling cycle) */
  tarch::la::Vector<dim, double>* _oldSolution;
  /** solution of velocity field at the end of coupling cycle */
  tarch::la::Vector<dim, double>* _newSolution;
  /** number of MD time steps per coupling cycle */
  const unsigned int _numberMDSteps;
  /** offset by which evaluation of time interval is shifted compared to
   * continuum solver. A value of zero implies interpolation between new and old
   * solution, corresponding to a time interval t->t+dt. A value of one implies
   * extrapolation between the times t+dt (newSolution) and t+2dt using
   * oldSolution and newSolution. All inter-/extrapolations are linear. */
  const double _shiftTimestep;
  /** local number of cells (incl. ghost layers) */
  const unsigned int _numberLocalCells;
  /** time step counter within a coupling cycle (should run from 0 to
   * _numberMDSteps) */
  unsigned int _timestepCounter;
  /** mass that was transferred in an earlier coupled step */
  double* _excessMass;
  /** true in each entry if west/east, south/north, bottom/top boundary is a
   * mass flux boundary */
  const tarch::la::Vector<2 * dim, bool> _massFluxBoundary;
};
#include "coupling/transferstrategies/TransferStrategy4NieCoupling.cpph"

#endif // _MOLECULARDYNAMICS_COUPLING_TRANSFERSTRATEGIES_TRANSFERSTRATEGY4NIECOUPLING_H_
