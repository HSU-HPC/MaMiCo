// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_TRANSFERSTRATEGIES_DIRECTTRANSFERSTRATEGY_H_
#define _MOLECULARDYNAMICS_COUPLING_TRANSFERSTRATEGIES_DIRECTTRANSFERSTRATEGY_H_

#include "coupling/cell-mappings/ComputeMassMapping.h"
#include "coupling/cell-mappings/ComputeMomentumMapping.h"
#include "coupling/transferstrategies/TransferStrategy.h"
#include "tarch/la/Vector.h"

namespace coupling {
namespace transferstrategies {
template <class LinkedCell, unsigned int dim> class DirectTransferStrategy;
}
} // namespace coupling

/** transfers and introduces mass and momentum directly into MD and to
 * macroscopic solver. So, for example, if mass M is coming from the macroscopic
 * solver, M is to be inserted into MD.
 *  @author Philipp Neumann
 *  @tparam LinkedCell the LinkedCell class is given by the implementation of
 * linked cells in the molecular dynamics simulation
 *  @tparam dim  refers to the spacial dimension of the simulation, can be 1, 2,
 * or 3 */
template <class LinkedCell, unsigned int dim>
class coupling::transferstrategies::DirectTransferStrategy
    : public coupling::transferstrategies::TransferStrategy<LinkedCell, dim> {
public:
  /** @brief a simple constructor
   *  @param mdSolverInterface interface for the md solver
   *  @param indexConversion an instance of the indexConversion */
  DirectTransferStrategy(coupling::interface::MDSolverInterface<LinkedCell, dim>
                             *const mdSolverInterface,
                         const coupling::IndexConversion<dim> &indexConversion)
      : coupling::transferstrategies::TransferStrategy<LinkedCell, dim>(
            mdSolverInterface, indexConversion),
        _massMapping(mdSolverInterface), _momentumMapping(mdSolverInterface) {}

  /** @brief a dummy destructor*/
  virtual ~DirectTransferStrategy() {}

  /** @brief the microscopicMass and -Momentum are set to 0
   *  @param cell macroscopic cell to process
   *  @param index index of the macroscopic cell */
  virtual void processInnerMacroscopicCellBeforeReceivingMacroscopicSolverData(
      coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell, dim>
          &cell,
      const unsigned int &index) {
    // reset quantities
    const tarch::la::Vector<dim, double> zero(0.0);
    cell.setMicroscopicMass(0.0);
    cell.setMicroscopicMomentum(zero);
  }

  /** @brief the microscopicMass and -Momentum are set to 0
   *  @param cell macroscopic cell to process
   *  @param index index of the macroscopic cell */
  virtual void processOuterMacroscopicCellBeforeReceivingMacroscopicSolverData(
      coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell, dim>
          &cell,
      const unsigned int &index) {
    // reset quantities
    const tarch::la::Vector<dim, double> zero(0.0);
    cell.setMicroscopicMass(0.0);
    cell.setMicroscopicMomentum(zero);
  }

  /** @brief the mass and momentum is evaluated for the cell and written to the
   * macroscopic quantities
   *  @param cell macroscopic cell to process
   *  @param index index of the macroscopic cell */
  virtual void processInnerMacroscopicCellBeforeSendingMDSolverData(
      coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell, dim>
          &cell,
      const unsigned int &index) {
    cell.iterateConstCells(_massMapping);
    cell.iterateConstCells(_momentumMapping);
    cell.setMacroscopicMass(_massMapping.getMass());
    cell.setMacroscopicMomentum(_momentumMapping.getMomentum());
  }

private:
  /** necessary to compute the mass in every single cell */
  coupling::cellmappings::ComputeMassMapping<LinkedCell, dim> _massMapping;
  /** necessary to compute the momentum in every single cell  */
  coupling::cellmappings::ComputeMomentumMapping<LinkedCell, dim>
      _momentumMapping;
};
#endif // _MOLECULARDYNAMICS_COUPLING_TRANSFERSTRATEGIES_DIRECTTRANSFERSTRATEGY_H_
