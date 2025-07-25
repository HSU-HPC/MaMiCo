// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_MOLECULEMAPPINGS_COMPUTEMEANVELOCITY_MAPPING_H_
#define _MOLECULARDYNAMICS_MOLECULEMAPPINGS_COMPUTEMEANVELOCITY_MAPPING_H_

#include "simplemd/Molecule.h"

namespace simplemd {
namespace moleculemappings {
class ComputeMeanVelocityMapping;
}
} // namespace simplemd

/** computes mean velocity.
 *
 *  @author Philipp Neumann
 */
class simplemd::moleculemappings::ComputeMeanVelocityMapping {
public:
  ComputeMeanVelocityMapping() {}
  ~ComputeMeanVelocityMapping() {}

  void beginMoleculeIteration() {
    _particleCounter = 0;
    for (unsigned int d = 0; d < MD_DIM; d++) {
      _meanVelocity[d] = 0.0;
    }
  }
  void endMoleculeIteration() { _meanVelocity = (1.0 / ((double)_particleCounter)) * _meanVelocity; }
  void handleMolecule(Molecule& molecule) {
    const tarch::la::Vector<MD_DIM, double>& velocity = molecule.getVelocity();
    _meanVelocity += velocity;
    _particleCounter++;
  }

  tarch::la::Vector<MD_DIM, double> getMeanVelocity() const { return _meanVelocity; }

  static const bool IsParallel = false;

private:
  tarch::la::Vector<MD_DIM, double> _meanVelocity;
  unsigned int _particleCounter;
};

#endif // _MOLECULARDYNAMICS_MOLECULEMAPPINGS_COMPUTEMEANVELOCITY_MAPPING_H_
