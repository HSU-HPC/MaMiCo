// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_MOMENTUMINSERTION_H_
#define _MOLECULARDYNAMICS_COUPLING_MOMENTUMINSERTION_H_

#include "coupling/datastructures/MacroscopicCell.h"
#include "tarch/la/Vector.h"

namespace coupling {
template <class LinkedCell, unsigned int dim> class MomentumInsertion;
}

/** used to manipulate the momentum/ velocity of the molecules contained in a
 * macroscopic cell.
 *
 *  @author Philipp Neumann
 */
template <class LinkedCell, unsigned int dim> class coupling::MomentumInsertion {
public:
  MomentumInsertion(coupling::interface::MDSolverInterface<LinkedCell, dim> *const mdSolverInterface) : _mdSolverInterface(mdSolverInterface) {}
  virtual ~MomentumInsertion() {}

  /** returns the number of MD steps between subsequent momentum insertions */
  virtual unsigned int getTimeIntervalPerMomentumInsertion() const = 0;

  /** inserts a fraction 'fraction' from the momentum of the macroscopic cell
   * 'cell' and distributes
   *  it over all molecules.
   *  This method does not conserve the kinetic energy of the respective
   * macroscopic cell. To conserve
   *  the energy as well, see the description of MomentumController on details
   * how to do that.
   */
  virtual void insertMomentum(coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell, dim> &cell,
                              const unsigned int &currentMacroscopicCellIndex) const = 0;

protected:
  coupling::interface::MDSolverInterface<LinkedCell, dim> *const _mdSolverInterface;
};

#endif // _MOLECULARDYNAMICS_COUPLING_MOMENTUMINSERTION_H_
