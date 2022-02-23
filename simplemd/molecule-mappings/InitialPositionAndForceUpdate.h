// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_MOLECULEMAPPINGS_INITIALPOSITIONANDFORCEUPDATE_H_
#define _MOLECULARDYNAMICS_MOLECULEMAPPINGS_INITIALPOSITIONANDFORCEUPDATE_H_

#include "simplemd/Molecule.h"

namespace simplemd {
namespace moleculemappings {
class InitialPositionAndForceUpdate;
}
} // namespace simplemd

/** does the initial update of forces and position.
 *
 *  @author Philipp Neumann
 */
class simplemd::moleculemappings::InitialPositionAndForceUpdate {
public:
  InitialPositionAndForceUpdate(const double &dt, const double &mass)
      : _dt(dt), _a(_dt / (2.0 * mass)), _zero(0.0) {}
  ~InitialPositionAndForceUpdate() {}

  void beginMoleculeIteration() {}
  void endMoleculeIteration() {}
  void handleMolecule(Molecule &molecule) {
    tarch::la::Vector<MD_DIM, double> &position = molecule.getPosition();
    position +=
        _dt * (molecule.getConstVelocity() + _a * molecule.getConstForce());
    molecule.setForceOld(molecule.getConstForce());
    molecule.setForce(_zero);
  }

private:
  const double _dt;
  const double _a;
  const tarch::la::Vector<MD_DIM, double> _zero;
};

#endif // _MOLECULARDYNAMICS_MOLECULEMAPPINGS_INITIALPOSITIONANDFORCEUPDATE_H_
