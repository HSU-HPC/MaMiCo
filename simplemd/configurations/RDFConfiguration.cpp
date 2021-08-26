// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#include "simplemd/configurations/RDFConfiguration.h"
#include "tarch/configuration/ParseConfiguration.h"

const std::string simplemd::configurations::RDFConfiguration::START_AT_TIMESTEP("start-at-timestep");
const std::string simplemd::configurations::RDFConfiguration::EVALUATE_EVERY_TIMESTEP("evaluate-every-timestep");
const std::string simplemd::configurations::RDFConfiguration::WRITE_EVERY_TIMESTEP("write-every-timestep");
const std::string simplemd::configurations::RDFConfiguration::NUMBER_OF_POINTS("number-of-points");

simplemd::configurations::RDFConfiguration::RDFConfiguration():
_startAtTimestep(0),_evaluateEveryTimestep(0),_writeEveryTimestep(0),_numberPoints(0),
_isDefined(false),_isValid(true){}



void simplemd::configurations::RDFConfiguration::
parseSubtag( tinyxml2::XMLElement* node ){
  // set is-defined
  _isDefined = true;
  int intBuf = -1;

  // parse write-every-timestep
  tarch::configuration::ParseConfiguration::readIntMandatory(intBuf,node,WRITE_EVERY_TIMESTEP);
  if (intBuf < 0){
    std::cout << WRITE_EVERY_TIMESTEP << " is smaller than zero: " << intBuf << std::endl;
    _isValid = false;
    exit(EXIT_FAILURE);
  }
  _writeEveryTimestep = (unsigned int) (intBuf);

  // parse evaluate-every-timestep
  tarch::configuration::ParseConfiguration::readIntMandatory(intBuf,node,EVALUATE_EVERY_TIMESTEP);
  if (intBuf < 0){
    std::cout << EVALUATE_EVERY_TIMESTEP << " is smaller than zero: " << intBuf << std::endl;
    _isValid = false;
    exit(EXIT_FAILURE);
  }
  _evaluateEveryTimestep = (unsigned int) (intBuf);

  // start-at-timestep
  tarch::configuration::ParseConfiguration::readIntMandatory(intBuf,node,START_AT_TIMESTEP);
  if (intBuf < 0){
    std::cout << START_AT_TIMESTEP << " is smaller than zero: " << intBuf << std::endl;
    _isValid = false;
    exit(EXIT_FAILURE);
  }
  _startAtTimestep = (unsigned int) (intBuf);

  tarch::configuration::ParseConfiguration::readIntMandatory(intBuf,node,NUMBER_OF_POINTS);
  if (intBuf < 1){
    std::cout << NUMBER_OF_POINTS << " is smaller than one: " << intBuf << std::endl;
    _isValid = false;
    exit(EXIT_FAILURE);
  }
  _numberPoints = (unsigned int) (intBuf);

  #if (MD_DEBUG==MD_YES)
  std::cout <<  "Start at step:       " << _startAtTimestep << std::endl;
  std::cout <<  "Write every step:    " << _writeEveryTimestep << std::endl;
  std::cout <<  "Evaluate every step: " << _evaluateEveryTimestep << std::endl;
  #endif
}


std::string simplemd::configurations::RDFConfiguration::getTag() const {
  return "rdf-configuration";
}


bool simplemd::configurations::RDFConfiguration::isValid() const {
  return _isValid;
}


