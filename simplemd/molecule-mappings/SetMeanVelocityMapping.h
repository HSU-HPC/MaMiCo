// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_MOLECULEMAPPINGS_SETMEANVELOCITY_MAPPING_H_
#define _MOLECULARDYNAMICS_MOLECULEMAPPINGS_SETMEANVELOCITY_MAPPING_H_

#include "simplemd/Molecule.h"

namespace simplemd {
namespace moleculemappings {
class SetMeanVelocityMapping;
}
} // namespace simplemd

/** sets a certain new mean velocity.
 *  oldVelocity corresponds to the old avg. velocity, newVelocity is the new
 * value.
 *
 *  @author Philipp Neumann
 */
class simplemd::moleculemappings::SetMeanVelocityMapping {
public:
  SetMeanVelocityMapping(const tarch::la::Vector<MD_DIM, double>& oldVelocity, const tarch::la::Vector<MD_DIM, double>& newVelocity)
      : _oldVelocity(oldVelocity), _newVelocity(newVelocity) {}
  ~SetMeanVelocityMapping() {}

  void beginMoleculeIteration() {}
  void endMoleculeIteration() {}
  void handleMolecule(Molecule& molecule) {
    tarch::la::Vector<MD_DIM, double>& velocity = molecule.getVelocity();
    velocity = velocity - _oldVelocity + _newVelocity;
  }

private:
  const tarch::la::Vector<MD_DIM, double> _oldVelocity;
  const tarch::la::Vector<MD_DIM, double> _newVelocity;
};

#endif // _MOLECULARDYNAMICS_MOLECULEMAPPINGS_SETMEANVELOCITY_MAPPING_H_
