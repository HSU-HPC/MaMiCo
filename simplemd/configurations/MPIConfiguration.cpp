// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#include "simplemd/configurations/MPIConfiguration.h"
#include "tarch/configuration/ParseConfiguration.h"

const std::string simplemd::configurations::MPIConfiguration::NUMBER_OF_PROCESSES("number-of-processes");

simplemd::configurations::MPIConfiguration::MPIConfiguration() : _numberOfProcesses(1), _isValid(true) {}

void simplemd::configurations::MPIConfiguration::parseSubtag(tinyxml2::XMLElement* node) {
  // parse number of processes
  tarch::la::Vector<MD_DIM, int> buffer(-1);
  tarch::configuration::ParseConfiguration::readVector<MD_DIM, int>(buffer, node, NUMBER_OF_PROCESSES);
  for (unsigned int d = 0; d < MD_DIM; d++) {
    if (buffer[d] <= 0) {
      std::cout << NUMBER_OF_PROCESSES << ": Entry " << d << " is smaller than or equal zero: " << buffer << std::endl;
      _isValid = false;
      exit(EXIT_FAILURE);
    }
    _numberOfProcesses[d] = (unsigned int)(buffer[d]);
  }

#if (MD_DEBUG == MD_YES)
  std::cout << "Number of processes: " << _numberOfProcesses << std::endl;
#endif

#if (MD_PARALLEL == MD_NO)
  unsigned int numProc = _numberOfProcesses[0];
  for (unsigned int d = 1; d < MD_DIM; d++) {
    numProc = numProc * _numberOfProcesses[d];
  }
  _isValid = _isValid && (numProc == 1);
  if (!_isValid)
    std::cout << "Invalid MPIConfiguration, MD_PARALLEL == MD_NO but _numberOfProcesses != 1" << std::endl;
#endif
}

std::string simplemd::configurations::MPIConfiguration::getTag() const { return "mpi-configuration"; }

bool simplemd::configurations::MPIConfiguration::isValid() const { return _isValid; }
