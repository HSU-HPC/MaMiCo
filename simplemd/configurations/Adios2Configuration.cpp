// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#include "simplemd/configurations/Adios2Configuration.h"
#include "tarch/configuration/ParseConfiguration.h"

const std::string simplemd::configurations::Adios2Configuration::FILENAME("directory-name");
const std::string simplemd::configurations::Adios2Configuration::WRITE_EVERY_TIMESTEP("write-every-timestep");

simplemd::configurations::Adios2Configuration::Adios2Configuration() : _filename(""), _writeEveryTimestep(0), _isValid(true) {}

void simplemd::configurations::Adios2Configuration::parseSubtag(tinyxml2::XMLElement* node) {
  int buffer = -1;
  // parse write-every-timestep
  tarch::configuration::ParseConfiguration::readIntMandatory(buffer, node, WRITE_EVERY_TIMESTEP);
  if (buffer < 0) {
    std::cout << WRITE_EVERY_TIMESTEP << " is smaller than zero: " << buffer << std::endl;
    _isValid = false;
    exit(EXIT_FAILURE);
  }
#ifndef BUILD_WITH_ADIOS2
  if (buffer != 0) {
    std::cout << " BUILD_WITH_ADIOS2 disabled but Adios2Writer enabled! " << std::endl;
    _isValid = false;
    exit(EXIT_FAILURE);
  }
#endif

  _writeEveryTimestep = (unsigned int)(buffer);

  // parse file stem
  tarch::configuration::ParseConfiguration::readStringMandatory(_filename, node, FILENAME);

#if (MD_DEBUG == MD_YES)
  std::cout << "Adios2-Directory-name:     " << _filename << std::endl;
  std::cout << "Write every step: " << _writeEveryTimestep << std::endl;
#endif
}

std::string simplemd::configurations::Adios2Configuration::getTag() const { return "adios2-configuration"; }

bool simplemd::configurations::Adios2Configuration::isValid() const { return _isValid; }