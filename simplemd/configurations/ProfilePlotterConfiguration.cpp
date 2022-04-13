// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#include "simplemd/configurations/ProfilePlotterConfiguration.h"
#include "tarch/configuration/ParseConfiguration.h"

const std::string simplemd::configurations::ProfilePlotterConfiguration::START_CELL("start-cell");
const std::string simplemd::configurations::ProfilePlotterConfiguration::RANGE("range");
const std::string simplemd::configurations::ProfilePlotterConfiguration::WRITE_EVERY_TIMESTEP("write-every-timestep");
const std::string simplemd::configurations::ProfilePlotterConfiguration::SAMPLE_EVERY_TIMESTEP("sample-every-timestep");
const std::string simplemd::configurations::ProfilePlotterConfiguration::START_AT_TIMESTEP("start-at-timestep");

simplemd::configurations::ProfilePlotterConfiguration::ProfilePlotterConfiguration()
    : _writeEveryTimestep(0), _sampleEveryTimestep(0), _startAtTimestep(0), _startCell(0), _range(0), _isValid(false) {}

void simplemd::configurations::ProfilePlotterConfiguration::parseSubtag(tinyxml2::XMLElement* node) {
  _isValid = true;
  tarch::la::Vector<MD_DIM, int> buf(-1);
  int intBuf = -1;

  tarch::configuration::ParseConfiguration::readVector<MD_DIM, int>(buf, node, START_CELL);
  for (unsigned int d = 0; d < MD_DIM; d++) {
    if (buf[d] < 0) {
      std::cout << "ERROR simplemd::configurations::ProfilePlotterConfiguration::parseSubtag: " << START_CELL.c_str() << " (entry " << d
                << ") is smaller than zero!" << std::endl;
      exit(EXIT_FAILURE);
      _isValid = false;
    }
    _startCell[d] = (unsigned int)(buf[d]);
  }

  tarch::configuration::ParseConfiguration::readVector<MD_DIM, int>(buf, node, RANGE);
  for (unsigned int d = 0; d < MD_DIM; d++) {
    if (buf[d] < 0) {
      std::cout << "ERROR simplemd::configurations::ProfilePlotterConfiguration::parseSubtag: " << RANGE.c_str() << " (entry " << d << ") is smaller than zero!"
                << std::endl;
      exit(EXIT_FAILURE);
      _isValid = false;
    }
    _range[d] = (unsigned int)(buf[d]);
  }

  tarch::configuration::ParseConfiguration::readIntMandatory(intBuf, node, WRITE_EVERY_TIMESTEP);
  if (intBuf < 0) {
    std::cout << "ERROR simplemd::configurations::ProfilePlotterConfiguration::parseSubtag: " << WRITE_EVERY_TIMESTEP.c_str() << " is smaller than zero!"
              << std::endl;
    exit(EXIT_FAILURE);
    _isValid = false;
  }
  _writeEveryTimestep = (unsigned int)(intBuf);

  tarch::configuration::ParseConfiguration::readIntMandatory(intBuf, node, SAMPLE_EVERY_TIMESTEP);
  if (intBuf < 0) {
    std::cout << "ERROR simplemd::configurations::ProfilePlotterConfiguration::parseSubtag: " << SAMPLE_EVERY_TIMESTEP.c_str() << " is smaller than zero!"
              << std::endl;
    exit(EXIT_FAILURE);
    _isValid = false;
  }
  _sampleEveryTimestep = (unsigned int)(intBuf);

  tarch::configuration::ParseConfiguration::readIntMandatory(intBuf, node, START_AT_TIMESTEP);
  if (intBuf < 0) {
    std::cout << "ERROR simplemd::configurations::ProfilePlotterConfiguration::parseSubtag: " << START_AT_TIMESTEP.c_str() << " is smaller than zero!"
              << std::endl;
    exit(EXIT_FAILURE);
    _isValid = false;
  }
  _startAtTimestep = (unsigned int)(intBuf);
}
