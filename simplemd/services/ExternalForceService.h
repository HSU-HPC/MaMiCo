// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_SERVICES_EXTERNALFORCESERVICE_H_
#define _MOLECULARDYNAMICS_SERVICES_EXTERNALFORCESERVICE_H_

#include "simplemd/MolecularDynamicsDefinitions.h"
#include "simplemd/Molecule.h"
#include "simplemd/configurations/ExternalForceConfiguration.h"
#include <cstdlib>
#include <iostream>
#include <vector>
#include <Kokkos_Core.hpp>

namespace simplemd {
namespace services {
class ExternalForceService;
}
} // namespace simplemd

/** adds external forces to the molecules
 *  @author Philipp Neumann
 */
class simplemd::services::ExternalForceService {
public:
  /** initialises the properties */
  ExternalForceService(const std::vector<simplemd::configurations::ExternalForceConfiguration>& configs);

  KOKKOS_FUNCTION void addExternalForce(tarch::la::Vector<MD_DIM, double>& force);

  ~ExternalForceService() {}

private:
  const std::vector<simplemd::configurations::ExternalForceConfiguration>& _configurations;
};
#endif // _MOLECULARDYNAMICS_SERVICES_EXTERNALFORCESERVICE_H_
