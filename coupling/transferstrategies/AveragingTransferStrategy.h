// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_TRANSFERSTRATEGIES_AVERAGINGTRANSFERSTRATEGY_H_
#define _MOLECULARDYNAMICS_COUPLING_TRANSFERSTRATEGIES_AVERAGINGTRANSFERSTRATEGY_H_

#include "coupling/cell-mappings/ComputeMassMapping.h"
#include "coupling/cell-mappings/ComputeMomentumMapping.h"
#include "coupling/transferstrategies/TransferStrategy.h"
#include <list>
#include <map>

namespace coupling {
/** transferstrategies define how the data will be applied, send, and received
 *  @brief the namespace is used for transferstrategies */
namespace transferstrategies {
template <class LinkedCell, unsigned int dim> class AveragingTransferStrategy;
}
} // namespace coupling

/** this class is used for pure averaging operations on the macroscopic cells.
 *  This can be used e.g. to measure errors in averaging over time, to estimate
 * number of samples etc.
 *  @author Philipp Neumann
 *  @tparam LinkedCell the LinkedCell class is given by the implementation of
 * linked cells in the molecular dynamics simulation
 *  @tparam dim  refers to the spacial dimension of the simulation, can be 1, 2,
 * or 3  */
template <class LinkedCell, unsigned int dim>
class coupling::transferstrategies::AveragingTransferStrategy
    : public coupling::transferstrategies::TransferStrategy<LinkedCell, dim> {
public:
  /** @brief a simple constructor
   *  @param mdSolverInterface interface for the md solver
   *  @param indexConversion instance of the indexConversion*/
  AveragingTransferStrategy(
      coupling::interface::MDSolverInterface<LinkedCell, dim>
          *const mdSolverInterface,
      const coupling::IndexConversion<dim> &indexConversion)
      : coupling::transferstrategies::TransferStrategy<LinkedCell, dim>(
            mdSolverInterface, indexConversion),
        _massMapping(mdSolverInterface), _momentumMapping(mdSolverInterface),
        _sampleCounter(0), _rank(indexConversion.getThisRank()) {}

  /** @brief a dummy destructor */
  virtual ~AveragingTransferStrategy() {}

  /** @brief reset the sample counter before processing any cell */
  virtual void
  beginProcessInnerMacroscopicCellsBeforeReceivingMacroscopicSolverData() {
    // reset sample counter for each coupling cycle
    _sampleCounter = 0;
  }

  /** @brief macroscopicMass and -Momentum are reset before the data from the
   * macro solver is transferred
   *  @param cell the macroscopic cell to process
   *  @param index the index of the macroscopic cell */
  virtual void processInnerMacroscopicCellBeforeReceivingMacroscopicSolverData(
      coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell, dim>
          &cell,
      const unsigned int &index) {
    // reset buffers for sampling mass and momentum in each inner macroscopic
    // cell
    cell.setMacroscopicMass(0.0);
    cell.setMacroscopicMomentum(tarch::la::Vector<dim, double>(0.0));
  }

  /** @brief values are reseted before the cells are processes and on rank=0
   * info is written to the stdstream */
  virtual void beginProcessInnerMacroscopicCellsAfterMDTimestep() {
    // output information of last sampling...
    if (_rank == 0) {
      std::cout << "Global quantities of sampling no. " << _sampleCounter
                << " on rank 0: mass=" << _avgMass
                << ", momentum=" << _avgMomentum << std::endl;
    }
    // reset avg. mass and momentum...
    _avgMass = 0.0;
    _avgMomentum = tarch::la::Vector<dim, double>(0.0);
    // and increment sample counter
    _sampleCounter++;
  }

  /** the macroscopicMass and -Momentum are averaged over all md time steps
   *  @brief the averaging operation is applied to the cell
   *  @param cell the macroscopic cell to process
   *  @param index the index of the macroscopic cell */
  virtual void processInnerMacroscopicCellAfterMDTimestep(
      coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell, dim>
          &cell,
      const unsigned int &index) {
    // compute total mass/momentum from previous samples
    const double oldMass = (_sampleCounter - 1) * cell.getMacroscopicMass();
    const tarch::la::Vector<dim, double> oldMomentum =
        ((double)(_sampleCounter - 1)) * cell.getMacroscopicMomentum();

    // compute new averaged mass and momentum
    cell.iterateConstCells(_massMapping);
    cell.iterateConstCells(_momentumMapping);
    const double mass =
        (1.0 / _sampleCounter) * (oldMass + _massMapping.getMass());
    const tarch::la::Vector<dim, double> momentum =
        (1.0 / _sampleCounter) * (oldMomentum + _momentumMapping.getMomentum());
    _avgMass += mass;
    _avgMomentum = _avgMomentum + momentum;
    // set mass and momentum in buffers
    cell.setMacroscopicMass(mass);
    cell.setMacroscopicMomentum(momentum);
  }

private:
  /** necessary to compute the mass in the cells*/
  coupling::cellmappings::ComputeMassMapping<LinkedCell, dim> _massMapping;
  /** necessary to compute the current momentum in the cell */
  coupling::cellmappings::ComputeMomentumMapping<LinkedCell, dim>
      _momentumMapping;
  /** counter for the samples*/
  unsigned int _sampleCounter;
  /** rank of the mpi process*/
  const unsigned int _rank;
  /** the momentum of every cell is summed in this variable and divided by the
   * _sampleCounter so it holds the average momentum in the end */
  tarch::la::Vector<dim, double> _avgMomentum;
  /** the mass of every cell is summed in this variable and divided by the
   * _sampleCounter so it holds the average mass in the end */
  double _avgMass;
};
#endif // _MOLECULARDYNAMICS_COUPLING_TRANSFERSTRATEGIES_AVERAGINGTRANSFERSTRATEGY_H_
