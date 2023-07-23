// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_CELLMAPPINGS_COMPUTETOTALPOTENTIALENERGYMAPPING_H_
#define _MOLECULARDYNAMICS_CELLMAPPINGS_COMPUTETOTALPOTENTIALENERGYMAPPING_H_

#include "simplemd/Molecule.h"
#include "simplemd/LinkedCell.h"
#include <cmath>

namespace simplemd {
namespace cellmappings {
class ComputeTotalPotentialEnergyMapping;
}
} // namespace simplemd

/** computes the total potential of a molecule from the partial two- and three-body-potentials in the adaptive resolution case.
 *
 *  @author Maximilian Mayr
 */
class simplemd::cellmappings::ComputeTotalPotentialEnergyMapping {
public:
  ComputeTotalPotentialEnergyMapping(const double& interfaceStart, const double& interfaceLength, const unsigned int& dimension)
      : _interfaceStart(interfaceStart), _interfaceLength(interfaceLength), _dimension(dimension) {}
  ~ComputeTotalPotentialEnergyMapping() {}

  void beginCellIteration() {}
  void endCellIteration() {}
  void handleCell(LinkedCell& cell, const unsigned int& cellIndex);
  void handleCellPair(LinkedCell& cell1, LinkedCell& cell2, const unsigned int& cellIndex1, const unsigned int& cellIndex2) {} 
  void handleCellTriplet(LinkedCell& cell1, LinkedCell& cell2, LinkedCell& cell3,
                         const unsigned int& cellIndex1, const unsigned int& cellIndex2, const unsigned int& cellIndex3) {}

private:
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