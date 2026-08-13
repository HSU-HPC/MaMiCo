// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_MOLECULEMAPPINGS_COMPUTEMEANVELOCITY_MAPPING_H_
#define _MOLECULARDYNAMICS_MOLECULEMAPPINGS_COMPUTEMEANVELOCITY_MAPPING_H_

#include "simplemd/Molecule.h"
#include "tarch/utils/Utils.h"

namespace simplemd {
namespace moleculemappings {
class ComputeMeanVelocityMapping;
}
} // namespace simplemd

/** computes mean velocity. Uses fixed-point math so that result does not depend on iteration order
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
      _meanVelocity[d] = 0;
    }
  }

  void endMoleculeIteration() { _meanVelocity = _meanVelocity / _particleCounter; }

  void handleMolecule(Molecule& molecule) {
    DEFINE_DECIMAL_FP_LIMITS(3);
    tarch::la::Vector<MD_DIM, double> v{stepFP3 * molecule.getVelocity()};
    for (unsigned int d = 0; d < MD_DIM; d++) {
      _meanVelocity[d] += (long long)(v[d]);
    }
    _particleCounter++;
  }

  tarch::la::Vector<MD_DIM, double> getMeanVelocity() const {
    DEFINE_DECIMAL_FP_LIMITS(3);
    tarch::la::Vector<MD_DIM, double> res;
    for (unsigned int d = 0; d < MD_DIM; d++) {
      res[d] = _meanVelocity[d] * minFP3;
    }
    return res;
  }

  static const bool IsParallel = false;

private:
  tarch::la::Vector<MD_DIM, long long> _meanVelocity;
  long long _particleCounter;
};

#endif // _MOLECULARDYNAMICS_MOLECULEMAPPINGS_COMPUTEMEANVELOCITY_MAPPING_H_
