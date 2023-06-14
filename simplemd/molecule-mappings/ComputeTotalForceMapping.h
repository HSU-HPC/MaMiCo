// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_MOLECULEMAPPINGS_COMPUTETOTALFORCEMAPPING_H_
#define _MOLECULARDYNAMICS_MOLECULEMAPPINGS_COMPUTETOTALFORCEMAPPING_H_

#include "simplemd/Molecule.h"
#include <cmath>

namespace simplemd {
namespace moleculemappings {
class ComputeTotalForceMapping;
}
} // namespace simplemd

/** computes the total force of a molecule from the partial two- and three-body-forces in the adaptive resolution case.
 *
 *  @author Maximilian Mayr
 */
class simplemd::moleculemappings::ComputeTotalForceMapping {
public:
  ComputeTotalForceMapping(const double& interfaceStart, const double& interfaceLength, const unsigned int& dimension)
      : _interfaceStart(interfaceStart), _interfaceLength(interfaceLength), _dimension(dimension), _zero(0.0) {}
  ~ComputeTotalForceMapping() {}

  void beginMoleculeIteration() {}
  void endMoleculeIteration() {}
  void handleMolecule(Molecule& molecule) {
    double position = molecule.getPosition()[_dimension - 1];
    double weight = calculateWeight(position);
    tarch::la::Vector<MD_DIM, double> totalForce = weight * molecule.getConstThreeBodyForce() + (1.0 - weight) * molecule.getConstTwoBodyForce();
    molecule.setForce(totalForce);
    molecule.setTwoBodyForce(_zero);
    molecule.setThreeBodyForce(_zero);
  }

private:
  const double _interfaceStart;
  const double _interfaceLength;
  const unsigned int _dimension;
  const tarch::la::Vector<MD_DIM, double> _zero;

  double calculateWeight(const double& position) {
    if (position <= _interfaceStart) {
        return 0.0;
    }
    else if (position >= _interfaceStart + _interfaceLength) {
        return 1.0;
    }
    else {
        return std::sin(MD_PI / (2.0 * _interfaceLength) * (position - _interfaceStart)) *
               std::sin(MD_PI / (2.0 * _interfaceLength) * (position - _interfaceStart));
    }
  }
};

#endif // _MOLECULARDYNAMICS_MOLECULEMAPPINGS_COMPUTETOTALFORCEMAPPING_H_