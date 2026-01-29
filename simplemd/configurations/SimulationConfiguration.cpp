// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#include "simplemd/configurations/SimulationConfiguration.h"
#include "tarch/configuration/ParseConfiguration.h"

const std::string simplemd::configurations::SimulationConfiguration::DT("dt");
const std::string simplemd::configurations::SimulationConfiguration::NUMBER_OF_TIMESTEPS("number-of-timesteps");
const std::string simplemd::configurations::SimulationConfiguration::COMPUTE_QUANTITIES_EVERY_TIMESTEP("compute-macroscopic-quantities-every-timestep");
const std::string simplemd::configurations::SimulationConfiguration::FIX_SEED("fix-seed");
const std::string
    simplemd::configurations::SimulationConfiguration::OVERLAP_COMMUNICATION_WITH_FORCE_COMPUTATION("overlap-communication-with-force-computation");

simplemd::configurations::SimulationConfiguration::SimulationConfiguration()
    : _dt(0.0), _numberOfTimesteps(0), _fixSeed(false), _overlapCommWithForceComputation(false), _isValid(true) {}

void simplemd::configurations::SimulationConfiguration::parseSubtag(tinyxml2::XMLElement* node) {
  int intBuf = -1;

  // get timestep size dt
  tarch::configuration::ParseConfiguration::readDoubleMandatory(_dt, node, DT);
  if (_dt <= 0.0) {
    std::cout << DT << " is smaller than or equal zero!" << std::endl;
    _isValid = false;
    exit(EXIT_FAILURE);
  }

  // get number of timesteps
  tarch::configuration::ParseConfiguration::readIntMandatory(intBuf, node, NUMBER_OF_TIMESTEPS);
  if (intBuf <= 0) {
    std::cout << NUMBER_OF_TIMESTEPS << " is smaller than or equal zero: " << intBuf << std::endl;
    _isValid = false;
    exit(EXIT_FAILURE);
  }
  _numberOfTimesteps = (unsigned int)(intBuf);

  // get quantity evaluation info
  tarch::configuration::ParseConfiguration::readIntMandatory(intBuf, node, COMPUTE_QUANTITIES_EVERY_TIMESTEP);
  if (intBuf < 0) {
    std::cout << REORGANISE_MEMORY_EVERY_TIMESTEP << " is smaller than zero: " << intBuf << std::endl;
    _isValid = false;
    exit(EXIT_FAILURE);
  }
  _computeMacroscopicQuantitiesEveryTimestep = (unsigned int)(intBuf);

  // read fix-seed parameter, if available; is kept at false otherwise
  _fixSeed = false;
  tarch::configuration::ParseConfiguration::readBoolOptional(_fixSeed, node, FIX_SEED);

  // read if in parallel simulation force and send operations shall be
  // overlapped
  _overlapCommWithForceComputation = false;
  tarch::configuration::ParseConfiguration::readBoolOptional(_overlapCommWithForceComputation, node, OVERLAP_COMMUNICATION_WITH_FORCE_COMPUTATION);

#if (MD_DEBUG == MD_YES)
  std::cout << "Timestep dt:   " << _dt << std::endl;
  std::cout << "No. timesteps: " << _numberOfTimesteps << std::endl;
#endif
}

std::string simplemd::configurations::SimulationConfiguration::getTag() const { return "simulation-configuration"; }

bool simplemd::configurations::SimulationConfiguration::isValid() const { return _isValid; }
