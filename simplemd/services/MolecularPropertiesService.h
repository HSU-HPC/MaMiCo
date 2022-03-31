// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_SERVICES_MOLECULARPROPERTIESSERVICE_H_
#define _MOLECULARDYNAMICS_SERVICES_MOLECULARPROPERTIESSERVICE_H_

#include <iostream>
#include <cstdlib>
#include "simplemd/MolecularProperties.h"
#include "simplemd/MolecularDynamicsDefinitions.h"

namespace simplemd {
namespace services { class MolecularPropertiesService; }
}

/** access to the molecular properties.
 *  @author Philipp Neumann
 */
class simplemd::services::MolecularPropertiesService {
public:
  /** initialises the properties */
  MolecularPropertiesService(const double &mass, const double &epsilon,
                             const double &sigma, const double &cutOffRadius,
                             const double &kB);

  /** returns the properties */
  const MolecularProperties &getMolecularProperties() const;

  void shutdown();

  ~MolecularPropertiesService() {}

private:
  const MolecularProperties _properties;
};
#endif // _MOLECULARDYNAMICS_SERVICES_MOLECULARPROPERTIESSERVICE_H_
