// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_NOMOMENTUMINSERTION_H_
#define _MOLECULARDYNAMICS_COUPLING_NOMOMENTUMINSERTION_H_

#include "coupling/MomentumInsertion.h"
namespace coupling {
template <class LinkedCell, unsigned int dim> class NoMomentumInsertion;
}

/** does not do anything with the momentum. Empty insertion mechanism.
 * @author Philipp Neumann
 */
template <class LinkedCell, unsigned int dim> class coupling::NoMomentumInsertion : public coupling::MomentumInsertion<LinkedCell, dim> {
public:
  NoMomentumInsertion(coupling::interface::MDSolverInterface<LinkedCell, dim> *const mdSolverInterface)
      : MomentumInsertion<LinkedCell, dim>(mdSolverInterface) {}
  virtual ~NoMomentumInsertion() {}

  /** returns the number of MD steps between subsequent momentum insertions */
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
                              const unsigned int &currentMacroscopicCellIndex) const {}
};

#endif // _MOLECULARDYNAMICS_COUPLING_NOMOMENTUMINSERTION_H_
