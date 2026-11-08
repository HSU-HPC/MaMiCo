// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_ADDITIVEMOMENTUMINSERTION_H_
#define _MOLECULARDYNAMICS_COUPLING_ADDITIVEMOMENTUMINSERTION_H_

#include "coupling/MomentumInsertion.h"
#include "coupling/cell-mappings/ComputeMassMapping.h"
#include "coupling/cell-mappings/SetMomentumMapping.h"
#include "coupling/datastructures/CouplingCell.h"
#include "tarch/la/Vector.h"

/** @brief everything necessary for coupling operations, is defined in here */
namespace coupling {
template <class LinkedCell, unsigned int dim> class AdditiveMomentumInsertion;
}

/** This class allows to add momentum to molecules. In each MD timestep, it
 * takes a respective fraction from the momentum buffer of a coupling cell
 * and adds this momentum to the molecules of the coupling cell.
 *  @brief used to manipulate the momentum/velocity of the molecules contained
 * in a coupling cell.
 *  @tparam LinkedCell the LinkedCell class is given by the implementation of
 * linked cells in the molecular dynamics simulation
 *  @tparam dim the integer dim refers to the spacial dimension of the
 * simulation, can be 1, 2, or 3
 *  @author Philipp Neumann
 */
template <class LinkedCell, unsigned int dim> class coupling::AdditiveMomentumInsertion : public coupling::MomentumInsertion<LinkedCell, dim> {
public:
  /** @brief A simple constructor
   *   @param mdSolverInterface The interface of the molecular dynamics solver
   * in application
   *   @param numberMDTimestepsPerCouplingCycle The number of molecular dynamics
   * timesteps within one coupling cycle */
  AdditiveMomentumInsertion(coupling::interface::MDSolverInterface<LinkedCell, dim>* const mdSolverInterface, unsigned int numberMDTimestepsPerCouplingCycle)
      : coupling::MomentumInsertion<LinkedCell, dim>(mdSolverInterface), _numberMDTimestepsPerCouplingCycle(numberMDTimestepsPerCouplingCycle) {}

  /** @brief A simple destructor*/
  ~AdditiveMomentumInsertion() {}

  /** returns 1 since momentum shall be inserted in every time step
   *  @brief returns the interval of time steps for momentum insertion
   *  @returns 1 */
  unsigned int getTimeIntervalPerMomentumInsertion() const override { return 1; }

  /** inserts the momentum of the coupling cell and distributes it over all
   * molecules. This method does not conserve the kinetic energy of the
   * respective coupling cell. To conserve the energy as well, see the
   * description of MomentumController on details how to do that.
   *  @brief inserts momentum to the cell
   *  @param cell coupling cell to insert momentum to
   *  @param idx coupling cell index of the cell */
  void insertMomentum(coupling::datastructures::CouplingCellWithLinkedCells<LinkedCell, dim>& cell, I02 idx) const override {
    const unsigned int timeIntervalMomentumInsertion = getTimeIntervalPerMomentumInsertion();
    // determine fraction of momentum that is to be inserted in this frame
    double fraction = 1.0 / ((_numberMDTimestepsPerCouplingCycle / timeIntervalMomentumInsertion) +
                             (_numberMDTimestepsPerCouplingCycle % timeIntervalMomentumInsertion != 0));
    tarch::la::Vector<dim, double> momentum(fraction * cell.getMicroscopicMomentum());
    tarch::la::Vector<dim, double> zeroMomentum(0.0);

    // if there is some momentum to be transferred, do so
    if (tarch::la::dot(momentum, momentum) != 0.0) {
      coupling::cellmappings::ComputeMassMapping<LinkedCell, dim> computeMassMapping(coupling::MomentumInsertion<LinkedCell, dim>::_mdSolverInterface);
      cell.iterateConstCells(computeMassMapping);
      unsigned int numberParticles = computeMassMapping.getNumberOfParticles();
      coupling::cellmappings::SetMomentumMapping<LinkedCell, dim> setMomentumMapping(zeroMomentum, momentum, numberParticles,
                                                                                     coupling::MomentumInsertion<LinkedCell, dim>::_mdSolverInterface);
      cell.iterateCells(setMomentumMapping);
    }
  }

  void setInnerImposition(bool enable) override {
    throw std::runtime_error(std::string("coupling::AdditiveMomentumInsertion::setInnerImposition not implemented"));
  }

private:
  const unsigned int _numberMDTimestepsPerCouplingCycle; ///< The number of molecular dynamics
                                                         ///< timesteps within one coupling
                                                         ///< cycle */
};
#endif // _MOLECULARDYNAMICS_COUPLING_ADDITIVEMOMENTUMINSERTION_H_
