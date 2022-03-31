// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#include "simplemd/configurations/MolecularDynamicsConfiguration.h"

void simplemd::configurations::MolecularDynamicsConfiguration::parseSubtag(tinyxml2::XMLElement *node) {
  const simplemd::configurations::ProfilePlotterConfiguration tagProfilePlotter;
  const simplemd::configurations::ExternalForceConfiguration externalForce;

  for (tinyxml2::XMLElement *child = node->FirstChildElement(); child != NULL; child = child->NextSiblingElement()) {
    std::string nodename(child->Value());
    if (nodename == _domainConfiguration.getTag()) {
      _domainConfiguration.parseSubtag(child);
    } else if (nodename == _moleculeConfiguration.getTag()) {
      _moleculeConfiguration.parseSubtag(child);
    } else if (nodename == _simulationConfiguration.getTag()) {
      _simulationConfiguration.parseSubtag(child);
    } else if (nodename == _vtkConfiguration.getTag()) {
      _vtkConfiguration.parseSubtag(child);
    } else if (nodename == _mpiConfiguration.getTag()) {
      _mpiConfiguration.parseSubtag(child);
    } else if (nodename == _rdfConfiguration.getTag()) {
      _rdfConfiguration.parseSubtag(child);
    } else if (nodename == _checkpointConfiguration.getTag()) {
      _checkpointConfiguration.parseSubtag(child);
    } else if (nodename == tagProfilePlotter.getTag()) {
      _profilePlotterConfigurations.push_back(simplemd::configurations::ProfilePlotterConfiguration());
      _profilePlotterConfigurations[_profilePlotterConfigurations.size() - 1].parseSubtag(child);
    } else if (nodename == externalForce.getTag()) {
      _externalForceConfigurations.push_back(simplemd::configurations::ExternalForceConfiguration());
      _externalForceConfigurations[_externalForceConfigurations.size() - 1].parseSubtag(child);
    } else {
      std::cout << "Unknown subtag " << nodename << "!" << std::endl;
      exit(EXIT_FAILURE);
    }
  }
}

std::string simplemd::configurations::MolecularDynamicsConfiguration::getTag() const { return "molecular-dynamics"; }

bool simplemd::configurations::MolecularDynamicsConfiguration::isValid() const {
  bool isValid = true;
  for (unsigned int i = 0; i < _profilePlotterConfigurations.size(); i++) {
    isValid = isValid && _profilePlotterConfigurations[i].isValid();
  }
  for (unsigned int i = 0; i < _externalForceConfigurations.size(); i++) {
    isValid = isValid && _externalForceConfigurations[i].isValid();
  }

  return isValid && _domainConfiguration.isValid() && _moleculeConfiguration.isValid() && _vtkConfiguration.isValid() && _simulationConfiguration.isValid() &&
         _mpiConfiguration.isValid() && _rdfConfiguration.isValid() && _checkpointConfiguration.isValid();
}
