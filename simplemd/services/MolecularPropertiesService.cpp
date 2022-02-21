// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#include "simplemd/services/MolecularPropertiesService.h"

simplemd::services::MolecularPropertiesService::MolecularPropertiesService(
    const double &mass, const double &epsilon, const double &sigma,
    const double &cutOffRadius, const double &kB)
    : _properties(mass, epsilon, sigma, cutOffRadius, kB) {
#if (MD_DEBUG == MD_YES)
  std::cout << "init MolecularProperties, mass: " << _properties.getMass()
            << std::endl;
  std::cout << " epsilon: " << _properties.getEpsilon() << std::endl;
  std::cout << " sigma: " << _properties.getSigma() << " cut-off rad. "
            << std::endl;
  std::cout << _properties.getCutOffRadius() << " kB: " << _properties.getKB()
            << std::endl;
#endif
}

const simplemd::MolecularProperties &
simplemd::services::MolecularPropertiesService::getMolecularProperties() const {
  return _properties;
}

void simplemd::services::MolecularPropertiesService::shutdown() {}
