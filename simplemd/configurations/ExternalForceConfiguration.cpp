// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#include "simplemd/configurations/ExternalForceConfiguration.h"

#include "tarch/configuration/ParseConfiguration.h"

const std::string simplemd::configurations::ExternalForceConfiguration::VALUE("value");

void simplemd::configurations::ExternalForceConfiguration::parseSubtag(tinyxml2::XMLElement* node) {
  _isValid = true;
  tarch::configuration::ParseConfiguration::readVectorMandatory<MD_DIM, double>(_externalForce, node, VALUE);
}
