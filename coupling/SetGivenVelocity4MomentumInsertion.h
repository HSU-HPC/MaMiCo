// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_SETGIVENVELOCITY4MOMENTUMINSERTION_H_
#define _MOLECULARDYNAMICS_COUPLING_SETGIVENVELOCITY4MOMENTUMINSERTION_H_

#include "coupling/MomentumInsertion.h"
#include "coupling/cell-mappings/ComputeMassMapping.h"
#include "coupling/cell-mappings/ComputeMomentumMapping.h"
#include "coupling/cell-mappings/SetMomentumMapping.h"
#include "coupling/interface/MDSolverInterface.h"

namespace coupling {
template <class LinkedCell, unsigned int dim> class SetGivenVelocity4MomentumInsertion;
}

/** interpretes the microscopicMomentum-buffer as velocity and sets this value
 * in the respective coupling cell.
 *  @author Philipp Neumann
 *  @tparam LinkedCell the LinkedCell class is given by the implementation of
 * linked cells in the molecular dynamics simulation
 *  @tparam dim  refers to the spacial dimension of the simulation, can be 1, 2,
 * or 3  */
template <class LinkedCell, unsigned int dim> class coupling::SetGivenVelocity4MomentumInsertion : public coupling::MomentumInsertion<LinkedCell, dim> {
public:
  /** @brief a simple constructor
   *  @param mdSolverInterface interface for the md solver */
  SetGivenVelocity4MomentumInsertion(coupling::interface::MDSolverInterface<LinkedCell, dim>* const mdSolverInterface)
      : coupling::MomentumInsertion<LinkedCell, dim>(mdSolverInterface) {}

  /** @brief a simple destructor */
  virtual ~SetGivenVelocity4MomentumInsertion() {}

  /** @brief returns 1, since momentum insertions will be applied in every md
   * timestep
   *  @returns the time step interval for momentum insertion */
  virtual unsigned int getTimeIntervalPerMomentumInsertion() const { return 1; }

  /** This method does not conserve the kinetic energy of the respective
   * coupling cell. To conserve the energy as well, see the description of
   * MomentumController on details how to do that.
   *  @brief updates the momentum based on the microscopic momentum
   *  @param cell coupling cell
   *  @param currentCouplingCellIndex index of the coupling cell
   */
  virtual void insertMomentum(coupling::datastructures::CouplingCellWithLinkedCells<LinkedCell, dim>& cell,
                              const unsigned int& currentCouplingCellIndex) const {
    coupling::cellmappings::ComputeMassMapping<LinkedCell, dim> massMapping(coupling::MomentumInsertion<LinkedCell, dim>::_mdSolverInterface);
    coupling::cellmappings::ComputeMomentumMapping<LinkedCell, dim> momentumMapping(coupling::MomentumInsertion<LinkedCell, dim>::_mdSolverInterface);

    // old and new momentum
    tarch::la::Vector<dim, double> oldMomentum(0.0);
    tarch::la::Vector<dim, double> newMomentum(0.0);
    unsigned int particleCounter = 0;

    cell.iterateConstCells(massMapping);
    cell.iterateConstCells(momentumMapping);
    // set old momentum
    oldMomentum = momentumMapping.getMomentum();
    // set new momentum (only number of particles is missing)
    newMomentum = massMapping.getMass() * cell.getMicroscopicMomentum();
    particleCounter = massMapping.getNumberOfParticles();

    // set new momentum (based on velocity stored in microscopic
    // momentum-buffer)
    coupling::cellmappings::SetMomentumMapping<LinkedCell, dim> setMomentumMapping(oldMomentum, newMomentum, particleCounter,
                                                                                   coupling::MomentumInsertion<LinkedCell, dim>::_mdSolverInterface);
    cell.iterateCells(setMomentumMapping);
  }
};
#endif // _MOLECULARDYNAMICS_COUPLING_SETGIVENVELOCITY4MOMENTUMINSERTION_H_
