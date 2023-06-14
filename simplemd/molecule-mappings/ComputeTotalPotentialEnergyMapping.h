// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_MOLECULEMAPPINGS_COMPUTETOTALPOTENTIALENERGYMAPPING_H_
#define _MOLECULARDYNAMICS_MOLECULEMAPPINGS_COMPUTETOTALPOTENTIALENERGYMAPPING_H_

#include "simplemd/Molecule.h"
#include <cmath>

namespace simplemd {
namespace moleculemappings {
class ComputeTotalPotentialEnergyMapping;
}
} // namespace simplemd

/** computes the total potential of a molecule from the partial two- and three-body-potentials in the adaptive resolution case.
 *
 *  @author Maximilian Mayr
 */
class simplemd::moleculemappings::ComputeTotalForceMapping {
public:
  ComputeTotalPotentialEnergyMapping(const double& domainSize, const double& interfaceStart, const double& interfaceLength, const unsigned int& dimension)
      : _domainSize(domainSize), _interfaceStart(interfaceStart), _interfaceLength(interfaceLength), _dimension(dimension) {}
  ~ComputeTotalPotentialEnergyMapping() {}

  void beginMoleculeIteration() {}
  void endMoleculeIteration() {}
  void handleMolecule(Molecule& molecule) {
    double position = molecule.getPosition()[dimension - 1];
    double weight = calculateWeight(position);
    double totalPotentialEnergy = weight * molecule.getConstThreeBodyPotentialEnergy() + (1.0 - weight) * molecule.getConstTwoBodyPotentialEnergy();
    molecule.setPotentialEnergy(totalPotentialEnergy);
    molecule.setTwoBodyPotentialEnergy(0.0);
    molecule.setThreeBodyPotentialEnergy(0.0);
  }

private:
  const double _domainSize;
  const double _interfaceStart;
  const double _interfaceLength;
  const unsigned int _dimension;

  double calculateWeight(const double& position) {
    if (position < _interfaceStart) {
        return 0.0;
    }
    else if (position > _interfaceStart + _interfaceLength) {
        return 1.0;
    }
    else {
        return std::sin(MD_PI / (2.0 * _interfaceLength) * (position - _interfaceStart)) *
               std::sin(MD_PI / (2.0 * _interfaceLength) * (position - _interfaceStart));
    }
  }
};

#endif // _MOLECULARDYNAMICS_MOLECULEMAPPINGS_COMPUTETOTALPOTENTIALENERGYMAPPING_H_