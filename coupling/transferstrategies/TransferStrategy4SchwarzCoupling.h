// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_TRANSFERSTRATEGIES_TRANSFERSTRATEGY4SCHWARZCOUPLING_H_
#define _MOLECULARDYNAMICS_COUPLING_TRANSFERSTRATEGIES_TRANSFERSTRATEGY4SCHWARZCOUPLING_H_

#include "coupling/cell-mappings/ComputeMassMapping.h"
#include "coupling/cell-mappings/ComputeMomentumMapping.h"
#include "coupling/transferstrategies/TransferStrategy.h"

namespace coupling {
namespace transferstrategies {
template <class LinkedCell, unsigned int dim> class TransferStrategy4SchwarzCoupling;
}
} // namespace coupling

/** transfer strategy for Schwarz coupling algorithm, adopted from Dupuis et al.
 *  We currently sample over the last 20% of the coupling interval, i.e. of the
 * numberMDsteps time steps in MD. The other 80% are used for equilibration.
 *  @author Philipp Neumann
 *  @tparam LinkedCell the LinkedCell class is given by the implementation of
 * linked cells in the molecular dynamics simulation
 *  @tparam dim  refers to the spacial dimension of the simulation, can be 1, 2,
 * or 3*/
template <class LinkedCell, unsigned int dim>
class coupling::transferstrategies::TransferStrategy4SchwarzCoupling : public coupling::transferstrategies::TransferStrategy<LinkedCell, dim> {
public:
  /** @brief a simple
   *  @param mdSolverInterface interface for the md solver
   *  @param numberMDSteps number of md time steps within one coupling time step
   */
  TransferStrategy4SchwarzCoupling(coupling::interface::MDSolverInterface<LinkedCell, dim>* const mdSolverInterface, unsigned int numberMDSteps)
      : coupling::transferstrategies::TransferStrategy<LinkedCell, dim>(mdSolverInterface), _massMapping(mdSolverInterface),
        _momentumMapping(mdSolverInterface), _timestepCounter(0), _sampleCounter(0), _numberMDSteps(numberMDSteps), _sampleEveryTimestep(1) {}

  /** @brief a dummy destructor*/
  virtual ~TransferStrategy4SchwarzCoupling() {}

  /** @brief the sample counter is reseted (0)*/
  void beginProcessInnerCouplingCellsBeforeReceivingMacroscopicSolverData() override;

  /** the momentum is converted to velocity (velocity = momentum/mass) and
   * stored in the microscopicMomentum the mass is not applied, therefore the
   * macroscopicMass is set to 0 the macroscopic quantities are reseted (0)
   *  @brief the data received from the macro solver is processed
   *  @param cell the coupling cell to be processed
   *  @param index the index of the coupling cell */
  void processInnerCouplingCellAfterReceivingMacroscopicSolverData(coupling::datastructures::CouplingCellWithLinkedCells<LinkedCell, dim>& cell,
                                                                   I02 index) override;

  /** @brief the momentum is converted to velocity (velocity = momentum/mass)
   * and stored in the microscopicMomentum the mass is not applied, therefore
   * the macroscopicMass is set to 0 the macroscopic quantities are reseted (0)
   *  @brief the data received from the macro solver is processed
   *  @param cell the coupling cell to be processed
   *  @param index the index of the coupling cell */
  void processOuterCouplingCellAfterReceivingMacroscopicSolverData(coupling::datastructures::CouplingCellWithLinkedCells<LinkedCell, dim>& cell,
                                                                   I02 index) override;

  /** @brief the data collected during the md time steps is averaged
   * (/numberMDSteps)
   *  @param cell the coupling cell to be processed
   *  @param index the index of the coupling cell */
  void processInnerCouplingCellBeforeSendingMDSolverData(coupling::datastructures::CouplingCellWithLinkedCells<LinkedCell, dim>& cell, I02 index) override;

  /** @brief the macroscopic quantities are reseted (0)
   *  @param cell the coupling cell to be processed
   *  @param index the index of the coupling cell */
  void processOuterCouplingCellBeforeSendingMDSolverData(coupling::datastructures::CouplingCellWithLinkedCells<LinkedCell, dim>& cell, I02 index) override;

  /** @brief the _timestepCounter is incremented and if sample()==true the
   * _sampleCounter too */
  void beginProcessInnerCouplingCellsAfterMDTimestep() override;

  /** @brief the mass and momentum are evaluated in the cell and stored in the
   * macroscopic quantities
   *  @param cell the coupling cell to be processed
   *  @param index the index of the coupling cell */
  void processInnerCouplingCellAfterMDTimestep(coupling::datastructures::CouplingCellWithLinkedCells<LinkedCell, dim>& cell, I02 index) override;

protected:
  /** @brief depending on the sampling interval (_sampleEveryTimestep) true is
   *  returned if sampling has to be done in the current time step */
  virtual bool sample() const;

private:
  /** computes the mass in every single cell */
  coupling::cellmappings::ComputeMassMapping<LinkedCell, dim> _massMapping;
  /** computes the momentum in every single cell */
  coupling::cellmappings::ComputeMomentumMapping<LinkedCell, dim> _momentumMapping;
  /** counter for the time steps */
  unsigned int _timestepCounter;
  /** counter for the samples */
  unsigned int _sampleCounter;
  /** number of md time steps within one coupling time step */
  const unsigned int _numberMDSteps;
  /** the interval of sample point; 1 means every md time step is sampled; 10
   * means every 10th is used */
  const int _sampleEveryTimestep;
};
#include "coupling/transferstrategies/TransferStrategy4SchwarzCoupling.cpph"

#endif // _MOLECULARDYNAMICS_COUPLING_TRANSFERSTRATEGIES_TRANSFERSTRATEGY4SCHWARZCOUPLING_H_
