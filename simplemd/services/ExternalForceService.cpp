// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#include "simplemd/services/ExternalForceService.h"

simplemd::services::ExternalForceService::ExternalForceService(
  const std::vector<simplemd::configurations::ExternalForceConfiguration>& configs
): _configurations(configs){}


void simplemd::services::ExternalForceService::addExternalForce(tarch::la::Vector<MD_DIM,double>& force){
  for (unsigned int i = 0; i < _configurations.size(); i++){
    force += _configurations[i].getExternalForce();
  }
}
