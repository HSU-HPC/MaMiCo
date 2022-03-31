// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_ADDITIVEMOMENTUMINSERTION_H_
#define _MOLECULARDYNAMICS_COUPLING_ADDITIVEMOMENTUMINSERTION_H_

#include "coupling/MomentumInsertion.h"
#include "coupling/cell-mappings/ComputeMassMapping.h"
#include "coupling/cell-mappings/SetMomentumMapping.h"
#include "coupling/datastructures/MacroscopicCell.h"
#include "tarch/la/Vector.h"

namespace coupling {
template <class LinkedCell, unsigned int dim> class AdditiveMomentumInsertion;
}

/** used to manipulate the momentum/ velocity of the molecules contained in a
 * macroscopic cell.
 *  This class allows to add momentum to molecules. In each MD timestep, it
 * takes a respective fraction
 *  from the momentum buffer of a macroscopic cell and adds this momentum to the
 * molecules of the
 *  macroscopic cell.
 *
 *  @author Philipp Neumann
 */
template <class LinkedCell, unsigned int dim> class coupling::AdditiveMomentumInsertion : public coupling::MomentumInsertion<LinkedCell, dim> {
public:
  AdditiveMomentumInsertion(coupling::interface::MDSolverInterface<LinkedCell, dim> *const mdSolverInterface, unsigned int numberMDTimestepsPerCouplingCycle)
      : coupling::MomentumInsertion<LinkedCell, dim>(mdSolverInterface), _numberMDTimestepsPerCouplingCycle(numberMDTimestepsPerCouplingCycle) {}
  ~AdditiveMomentumInsertion() {}

  /** insert momentum each timestep */
  virtual unsigned int getTimeIntervalPerMomentumInsertion() const { return 1; }

  /** inserts a fraction 'fraction' from the momentum of the macroscopic cell
   * 'cell' and distributes
   *  it over all molecules.
   *  This method does not conserve the kinetic energy of the respective
   * macroscopic cell. To conserve
   *  the energy as well, see the description of MomentumController on details
   * how to do that.
   */
  virtual void insertMomentum(coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell, dim> &cell,
                              const unsigned int &currentMacroscopicCellIndex) const {
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

private:
  const unsigned int _numberMDTimestepsPerCouplingCycle;
};
#endif // _MOLECULARDYNAMICS_COUPLING_ADDITIVEMOMENTUMINSERTION_H_
