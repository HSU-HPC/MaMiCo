// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#include "coupling/configurations/ParticleInsertionConfiguration.h"

const std::string coupling::configurations::ParticleInsertionConfiguration::INSERT_DELETE_MASS_EVERY_TIMESTEP("insert-every-timestep");
const std::string coupling::configurations::ParticleInsertionConfiguration::RSIGMA_COEFF("r_sigma");
const std::string coupling::configurations::ParticleInsertionConfiguration::MEAN_POTENTIAL_ENERGY_FACTOR("u_0-factor");
const std::string coupling::configurations::ParticleInsertionConfiguration::UOVERLAP_COEFF("u_ovlp");
const std::string coupling::configurations::ParticleInsertionConfiguration::STEPREF_COEFF("reference-stepsize");
const std::string coupling::configurations::ParticleInsertionConfiguration::ITER_MAX("maximum-number-of-iterations");
const std::string coupling::configurations::ParticleInsertionConfiguration::RESTART_MAX("maximum-number-of-restarts");
const std::string coupling::configurations::ParticleInsertionConfiguration::TOLERANCE("tolerance");
const std::string coupling::configurations::ParticleInsertionConfiguration::OFFSET_FROM_OUTER_BOUNDARY("offset-from-outer-boundary");

void coupling::configurations::ParticleInsertionConfiguration::parseSubtag(tinyxml2::XMLElement *node) {
  int buf = -1;

  // specify the type of insertion
  std::string thisType;
  tarch::configuration::ParseConfiguration::readStringMandatory(thisType, node, "type");
  // in this case, no insertion should be carried out, no more parameters need
  // to be parsed -> thus: return
  if (thisType == "none") {
    _particleInsertionType = NO_INSERTION;
    return;
  } else if (thisType == "usher") {
    _particleInsertionType = USHER;
  } else {
    std::cout << "ERROR "
                 "coupling::configurations::ParticleInsertionConfiguration::par"
                 "seSubtag(): 'type' not specified!"
              << std::endl;
    exit(EXIT_FAILURE);
  }

  buf = 1;
  tarch::configuration::ParseConfiguration::readIntOptional(buf, node, INSERT_DELETE_MASS_EVERY_TIMESTEP);
  if (buf <= 0) {
    std::cout << "ERROR coupling::configurations::ParticleInsertionConfiguration::";
    std::cout << "parseSubtag(): " << INSERT_DELETE_MASS_EVERY_TIMESTEP << " smaller than or equal zero!" << std::endl;
    exit(EXIT_FAILURE);
  } else {
    _insertDeleteMassEveryTimestep = (unsigned int)buf;
  }

  _rSigmaCoeff = 0.9;
  tarch::configuration::ParseConfiguration::readDoubleOptional(_rSigmaCoeff, node, RSIGMA_COEFF);
  if (_rSigmaCoeff <= 0.0) {
    std::cout << "ERROR coupling::configurations::ParticleInsertionConfiguration::";
    std::cout << "parseSubtag(): " << RSIGMA_COEFF << " smaller than or equal zero!" << std::endl;
    exit(EXIT_FAILURE);
  }

  _meanPotentialEnergyFactor = 1.0;
  tarch::configuration::ParseConfiguration::readDoubleOptional(_meanPotentialEnergyFactor, node, MEAN_POTENTIAL_ENERGY_FACTOR);
  if (_meanPotentialEnergyFactor <= 0.0) {
    std::cout << "ERROR coupling::configurations::ParticleInsertionConfiguration::";
    std::cout << "parseSubtag(): " << MEAN_POTENTIAL_ENERGY_FACTOR << " smaller than or equal zero!" << std::endl;
    exit(EXIT_FAILURE);
  }

  // default setting of last commit
  _tolerance = _meanPotentialEnergyFactor * 2.0;
  tarch::configuration::ParseConfiguration::readDoubleOptional(_tolerance, node, TOLERANCE);
  if (_tolerance <= 0.0) {
    std::cout << "ERROR "
                 "coupling::configurations::ParticleInsertionConfiguration::par"
                 "seSubtag(): "
              << TOLERANCE << " is smaller than or equal zero!" << std::endl;
    exit(EXIT_FAILURE);
  }

  _uOverlapCoeff = 10000.0;
  tarch::configuration::ParseConfiguration::readDoubleOptional(_uOverlapCoeff, node, UOVERLAP_COEFF);
  if (_uOverlapCoeff <= 0.0) {
    std::cout << "ERROR coupling::configurations::ParticleInsertionConfiguration::";
    std::cout << "parseSubtag(): " << UOVERLAP_COEFF << " smaller than or equal zero!" << std::endl;
    exit(EXIT_FAILURE);
  }

  // if not defined, choose optimal ref. step strategy (see usher-paper)
  _stepRefCoeff = -1.0;
  tarch::configuration::ParseConfiguration::readDoubleOptional(_stepRefCoeff, node, STEPREF_COEFF);
  if ((_stepRefCoeff != -1.0) && (_stepRefCoeff <= 0.0)) {
    std::cout << "ERROR coupling::configurations::ParticleInsertionConfiguration::";
    std::cout << "parseSubtag(): " << STEPREF_COEFF << " smaller than or equal zero!" << std::endl;
    exit(EXIT_FAILURE);
  }

  // if not defined, set offset for particle insertions from outer boundary to
  // 0.0; may yield instabilities for dense fluids
  _offsetFromOuterBoundary = 0.0;
  tarch::configuration::ParseConfiguration::readDoubleOptional(_offsetFromOuterBoundary, node, OFFSET_FROM_OUTER_BOUNDARY);
  if (_offsetFromOuterBoundary < 0.0) {
    std::cout << "ERROR coupling::configurations::ParticleInsertionConfiguration::";
    std::cout << "parseSubtag(): " << OFFSET_FROM_OUTER_BOUNDARY << " smaller than zero!" << std::endl;
    exit(EXIT_FAILURE);
  }

  tarch::configuration::ParseConfiguration::readIntMandatory(buf, node, ITER_MAX);
  if (buf <= 0) {
    std::cout << "ERROR coupling::configurations::ParticleInsertionConfiguration::";
    std::cout << "parseSubtag(): " << ITER_MAX << " smaller than or equal zero!" << std::endl;
    exit(EXIT_FAILURE);
  } else {
    _iterMax = (unsigned int)buf;
  }

  tarch::configuration::ParseConfiguration::readIntMandatory(buf, node, RESTART_MAX);
  if (buf <= 0) {
    std::cout << "ERROR coupling::configurations::ParticleInsertionConfiguration::";
    std::cout << "parseSubtag(): " << RESTART_MAX << " smaller or equal zero!" << std::endl;
    exit(EXIT_FAILURE);
  } else {
    _restartMax = (unsigned int)buf;
  }
}

std::string coupling::configurations::ParticleInsertionConfiguration::getTag() const { return "particle-insertion"; }

bool coupling::configurations::ParticleInsertionConfiguration::isValid() const { return true; }
