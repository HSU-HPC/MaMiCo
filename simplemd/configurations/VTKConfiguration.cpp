// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#include "simplemd/configurations/VTKConfiguration.h"
#include "tarch/configuration/ParseConfiguration.h"

const std::string simplemd::configurations::VTKConfiguration::FILENAME("filename");
const std::string simplemd::configurations::VTKConfiguration::WRITE_EVERY_TIMESTEP("write-every-timestep");

simplemd::configurations::VTKConfiguration::VTKConfiguration() : _filename(""), _writeEveryTimestep(0), _isValid(true) {}

void simplemd::configurations::VTKConfiguration::parseSubtag(tinyxml2::XMLElement* node) {
  int buffer = -1;
  // parse write-every-timestep
  tarch::configuration::ParseConfiguration::readIntMandatory(buffer, node, WRITE_EVERY_TIMESTEP);
  if (buffer < 0) {
    std::cout << WRITE_EVERY_TIMESTEP << " is smaller than zero: " << buffer << std::endl;
    _isValid = false;
    exit(EXIT_FAILURE);
  }
  _writeEveryTimestep = (unsigned int)(buffer);

  // parse file stem
  tarch::configuration::ParseConfiguration::readStringMandatory(_filename, node, FILENAME);

#if (MD_DEBUG == MD_YES)
  std::cout << "VTK-filename:     " << _filename << std::endl;
  std::cout << "Write every step: " << _writeEveryTimestep << std::endl;
#endif
}

std::string simplemd::configurations::VTKConfiguration::getTag() const { return "vtk-configuration"; }

bool simplemd::configurations::VTKConfiguration::isValid() const { return _isValid; }
