// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#include "simplemd/configurations/MoleculeConfiguration.h"

#include "tarch/configuration/ParseConfiguration.h"

const std::string
    simplemd::configurations::MoleculeConfiguration::MEAN_VELOCITY(
        "mean-velocity");
const std::string
    simplemd::configurations::MoleculeConfiguration::TEMPERATURE("temperature");
const std::string simplemd::configurations::MoleculeConfiguration::MASS("mass");
const std::string
    simplemd::configurations::MoleculeConfiguration::EPSILON("epsilon");
const std::string
    simplemd::configurations::MoleculeConfiguration::SIGMA("sigma");

simplemd::configurations::MoleculeConfiguration::MoleculeConfiguration()
    : _meanVelocity(0.0), _temperature(0.0), _mass(0.0), _epsilon(0.0),
      _sigma(0.0), _isValid(true) {}

void simplemd::configurations::MoleculeConfiguration::parseSubtag(
    tinyxml2::XMLElement *node) {
  // read mean velocity
  tarch::configuration::ParseConfiguration::readVector<MD_DIM, double>(
      _meanVelocity, node, MEAN_VELOCITY);

  // read temperature
  tarch::configuration::ParseConfiguration::readDoubleMandatory(
      _temperature, node, TEMPERATURE);
  if (_temperature <= 0.0) {
    std::cout << TEMPERATURE << " is smaller than or equal zero!" << std::endl;
    _isValid = false;
    exit(EXIT_FAILURE);
  }

  // read mass
  tarch::configuration::ParseConfiguration::readDoubleMandatory(_mass, node,
                                                                MASS);
  if (_mass <= 0.0) {
    std::cout << MASS << " is smaller than or equal zero!" << std::endl;
    _isValid = false;
    exit(EXIT_FAILURE);
  }

  // read epsilon
  tarch::configuration::ParseConfiguration::readDoubleMandatory(_epsilon, node,
                                                                EPSILON);
  if (_epsilon <= 0.0) {
    std::cout << EPSILON << " is smaller than or equal zero!" << std::endl;
    _isValid = false;
    exit(EXIT_FAILURE);
  }

  // read sigma
  tarch::configuration::ParseConfiguration::readDoubleMandatory(_sigma, node,
                                                                SIGMA);
  if (_sigma <= 0.0) {
    std::cout << SIGMA << " is smaller than or equal zero!" << std::endl;
    _isValid = false;
    exit(EXIT_FAILURE);
  }

#if (MD_DEBUG == MD_YES)
  std::cout << "Temperature: " << _temperature << std::endl;
  std::cout << "Sigma:       " << _sigma << std::endl;
  std::cout << "Epsilon:     " << _epsilon << std::endl;
  std::cout << "Mass:        " << _mass << std::endl;
  std::cout << "Mean vel.:   " << _meanVelocity << std::endl;
#endif
}

std::string simplemd::configurations::MoleculeConfiguration::getTag() const {
  return "molecule-configuration";
}

bool simplemd::configurations::MoleculeConfiguration::isValid() const {
  return _isValid;
}
