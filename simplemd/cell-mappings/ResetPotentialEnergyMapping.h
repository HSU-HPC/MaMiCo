// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_CELLMAPPINGS_RESETPOTENTIALENERGYMAPPING_H_
#define _MOLECULARDYNAMICS_CELLMAPPINGS_RESETPOTENTIALENERGYMAPPING_H_

#include "simplemd/LinkedCell.h"

namespace simplemd {
namespace cellmappings {
class ResetPotentialEnergyMapping;
}
} // namespace simplemd

/** resets the potential energy for all molecules within a certain linked cell.
 *  @author Philipp Neumann
 */
class simplemd::cellmappings::ResetPotentialEnergyMapping {
public:
  ResetPotentialEnergyMapping() : _zero(0.0) {}
  ~ResetPotentialEnergyMapping() {}

  void beginCellIteration() {}
  void endCellIteration() {}
  void handleCell(LinkedCell &cell, const unsigned int &cellIndex) {
    for (std::list<Molecule *>::const_iterator it = cell.begin(); it != cell.end(); it++) {
      (*it)->setPotentialEnergy(_zero);
    }
  }

private:
  const double _zero;
};

#endif // _MOLECULARDYNAMICS_CELLMAPPINGS_RESETPOTENTIALENERGYMAPPING_H_
