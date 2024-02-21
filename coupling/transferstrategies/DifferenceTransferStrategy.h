// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_TRANSFERSTRATEGIES_DIFFERENCETRANSFERSTRATEGY_H_
#define _MOLECULARDYNAMICS_COUPLING_TRANSFERSTRATEGIES_DIFFERENCETRANSFERSTRATEGY_H_

#include "coupling/cell-mappings/ComputeMassMapping.h"
#include "coupling/cell-mappings/ComputeMomentumMapping.h"
#include "coupling/transferstrategies/TransferStrategy.h"
#include "tarch/la/Vector.h"

namespace coupling {
namespace transferstrategies {
template <class LinkedCell, unsigned int dim> class DifferenceTransferStrategy;
}
} // namespace coupling

/** from MD to macroscopic solver: mass and momentum are computed at the current
 * timestep and mapped to macroscopic solver. The values are averaged over a
 * certain time. from macroscopic solver to MD: The difference in mass and
 * momentum between both solvers is computed and introduced on MD side.
 *  @author Philipp Neumann
 *  @tparam LinkedCell the LinkedCell class is given by the implementation of
 * linked cells in the molecular dynamics simulation
 *  @tparam dim  refers to the spacial dimension of the simulation, can be 1, 2,
 * or 3*/
template <class LinkedCell, unsigned int dim>
class coupling::transferstrategies::DifferenceTransferStrategy : public coupling::transferstrategies::TransferStrategy<LinkedCell, dim> {
public:
  /** @brief a simple constructor
   *  @param mdSolverInterface interface for the md solver
   *  @param indexConversion instance of the indexConversion
   *  @param numberMDSteps number of md steps within one coupling time step */
  DifferenceTransferStrategy(coupling::interface::MDSolverInterface<LinkedCell, dim>* const mdSolverInterface,
                             const coupling::IndexConversion<dim>& indexConversion, unsigned int numberMDsteps)
      : coupling::transferstrategies::TransferStrategy<LinkedCell, dim>(mdSolverInterface, indexConversion), _numberMDsteps(numberMDsteps), _zero(0.0),
        _massMapping(mdSolverInterface), _momentumMapping(mdSolverInterface) {}

  /** @brief a dummy destructor */
  virtual ~DifferenceTransferStrategy() {}

  /** @brief
   *  @param cell macroscopic cell to process
   *  @param index index of the macroscopic cell */
  virtual void processInnerCouplingCellBeforeReceivingMacroscopicSolverData(coupling::datastructures::CouplingCellWithLinkedCells<LinkedCell, dim>& cell,
                                                                               const unsigned int& index) {
    // reset quantities
    cell.setMicroscopicMass(0.0);
    cell.setMicroscopicMomentum(_zero);
  }

  /** @brief the microscopicMass and -Momentum are reseted to zero
   *  @param cell macroscopic cell to process
   *  @param index index of the macroscopic cell */
  virtual void processOuterCouplingCellBeforeReceivingMacroscopicSolverData(coupling::datastructures::CouplingCellWithLinkedCells<LinkedCell, dim>& cell,
                                                                               const unsigned int& index) {
    // reset quantities
    cell.setMicroscopicMass(0.0);
    cell.setMicroscopicMomentum(_zero);
  }

  /** @brief difference between microscopic and macroscopic values is evaluated
   * and stored in the macroscopic quantity
   *  @param cell macroscopic cell to process
   *  @param index index of the macroscopic cell */
  virtual void processInnerCouplingCellAfterReceivingMacroscopicSolverData(coupling::datastructures::CouplingCellWithLinkedCells<LinkedCell, dim>& cell,
                                                                              const unsigned int& index) {
    // compute difference between macroscopic and microscopic mass and momentum
    // values; this value is set in the microscopic data buffer of the
    // macroscopic cell.
    const double diffMass = cell.getMicroscopicMass() - cell.getMacroscopicMass();
    const tarch::la::Vector<dim, double> diffMomentum = cell.getMicroscopicMomentum() - cell.getMacroscopicMomentum();
    cell.setMicroscopicMass(diffMass);
    cell.setMicroscopicMomentum(diffMomentum);
    // reset macroscopic-buffers
    cell.setMacroscopicMass(0.0);
    cell.setMacroscopicMomentum(_zero);
  }

  /** @brief the quantities (mass & momentum) are averaged (divided by the
   * amount of md time steps)
   *  @param cell macroscopic cell to process
   *  @param index index of the macroscopic cell */
  virtual void processInnerCouplingCellBeforeSendingMDSolverData(coupling::datastructures::CouplingCellWithLinkedCells<LinkedCell, dim>& cell,
                                                                    const unsigned int& index) {
    // average quantities
    const double mass = cell.getMacroscopicMass() / ((double)_numberMDsteps);
    const tarch::la::Vector<dim, double> momentum = cell.getMacroscopicMomentum() * (1.0 / ((double)_numberMDsteps));
    cell.setMacroscopicMass(mass);
    cell.setMacroscopicMomentum(momentum);
  }

  /** @brief compute current mass and momentum and add it to averaged buffer
   * value
   *  @param cell macroscopic cell to process
   *  @param index index of the macroscopic cell */
  virtual void processInnerCouplingCellAfterMDTimestep(coupling::datastructures::CouplingCellWithLinkedCells<LinkedCell, dim>& cell,
                                                          const unsigned int& index) {
    cell.iterateConstCells(_massMapping);
    const double mass = _massMapping.getMass();
    cell.iterateConstCells(_momentumMapping);
    const tarch::la::Vector<dim, double> momentum = _momentumMapping.getMomentum();
    cell.addMacroscopicMass(mass);
    cell.addMacroscopicMomentum(momentum);
  }

private:
  /** number of md time steps within one coupling time step  */
  const unsigned int _numberMDsteps;
  /** vector containing zeros, the dimension is the spacial dimension of the
   * simulation */
  const tarch::la::Vector<dim, double> _zero;
  /** class to compute the mass within every single cell */
  coupling::cellmappings::ComputeMassMapping<LinkedCell, dim> _massMapping;
  /** class to compute the momentum within every single cell */
  coupling::cellmappings::ComputeMomentumMapping<LinkedCell, dim> _momentumMapping;
};
#endif // _MOLECULARDYNAMICS_COUPLING_TRANSFERSTRATEGIES_DIFFERENCETRANSFERSTRATEGY_H_
