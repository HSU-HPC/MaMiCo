// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#include "simplemd/configurations/DomainConfiguration.h"

#include "simplemd/MolecularDynamicsDefinitions.h"
#include "tarch/configuration/ParseConfiguration.h"

const std::string simplemd::configurations::DomainConfiguration::MOLECULES_PER_DIRECTION("molecules-per-direction");
const std::string simplemd::configurations::DomainConfiguration::DOMAIN_SIZE("domain-size");
const std::string simplemd::configurations::DomainConfiguration::DOMAIN_OFFSET("domain-offset");
const std::string simplemd::configurations::DomainConfiguration::CUTOFF_RADIUS("cutoff-radius");
const std::string simplemd::configurations::DomainConfiguration::LINKED_CELL_SIZE("linked-cell-size");
const std::string simplemd::configurations::DomainConfiguration::K_B("k_B");
const std::string simplemd::configurations::DomainConfiguration::BLOCK_SIZE("block-size");
const std::string simplemd::configurations::DomainConfiguration::BOUNDARY[MD_LINKED_CELL_NEIGHBOURS]
#if (MD_DIM == 1)
    = {"west", "east"}
#elif (MD_DIM == 2)
    = {"south-west", "south", "south-east", "west", "east", "north-west", "north", "north-east"}
#elif (MD_DIM == 3)
    = {"bottom-south-west",
       "bottom-south",
       "bottom-south-east",
       "bottom-west",
       "bottom",
       "bottom-east",
       "bottom-north-west",
       "bottom-north",
       "bottom-north-east",
       "south-west",
       "south",
       "south-east",
       "west",
       "east",
       "north-west",
       "north",
       "north-east",
       "top-south-west",
       "top-south",
       "top-south-east",
       "top-west",
       "top",
       "top-east",
       "top-north-west",
       "top-north",
       "top-north-east"}
#endif
;
const std::string simplemd::configurations::DomainConfiguration::PERIODIC_BOUNDARY("periodic");
const std::string simplemd::configurations::DomainConfiguration::GEOMETRY_BOUNDARY("geometry");
const std::string simplemd::configurations::DomainConfiguration::OPEN_BOUNDARY("open");
const std::string simplemd::configurations::DomainConfiguration::REFLECTING_BOUNDARY("reflecting");
const std::string simplemd::configurations::DomainConfiguration::RDF_FILENAME("radial-distribution-function");
const std::string simplemd::configurations::DomainConfiguration::CELLS_PER_LINKED_CELL("cells-per-linked-cell");
const std::string simplemd::configurations::DomainConfiguration::INIT_FROM_CHECKPOINT("init-from-checkpoint");
const std::string simplemd::configurations::DomainConfiguration::INIT_FROM_SEQUENTIAL_CHECKPOINT("init-from-sequential-checkpoint");
const std::string simplemd::configurations::DomainConfiguration::LINKED_CELLS_PER_NUMBER_DENSITY_EVALUATION("linked-cells-per-number-density-evaluation");

simplemd::configurations::DomainConfiguration::DomainConfiguration()
    : _moleculesPerDirection(0), _domainSize(0.0), _domainOffset(0.0), _cutoffRadius(0.0), _meshWidth(0.0), _kB(0.0), _blockSize(0),
      _boundary(simplemd::PERIODIC_BOUNDARY), _checkpointFilestem(""), _initFromCheckpoint(false), _initFromSequentialCheckpoint(false), _isValid(true) {}

void simplemd::configurations::DomainConfiguration::parseSubtag(tinyxml2::XMLElement* node) {
  std::string strBuf = "";
  int intBuf = -1;

  // get checkpoint info
  _initFromCheckpoint = false;
  _checkpointFilestem = "";
  tarch::configuration::ParseConfiguration::readStringOptional(_checkpointFilestem, node, INIT_FROM_CHECKPOINT);
  _initFromCheckpoint = (_checkpointFilestem != "");
  if (!_initFromCheckpoint) {
    // std::cout <<  INIT_FROM_CHECKPOINT << " not defined..." << std::endl;
  } else {
    // std::cout << "Initialise from checkpoint using filestem " <<
    // _checkpointFilestem << std::endl;
  }

  // get checkpoint info
  _initFromSequentialCheckpoint = false;
  strBuf = "";
  tarch::configuration::ParseConfiguration::readStringOptional(strBuf, node, INIT_FROM_SEQUENTIAL_CHECKPOINT);
  _initFromSequentialCheckpoint = (strBuf != "");
  if (!_initFromSequentialCheckpoint) {
    // std::cout <<  INIT_FROM_SEQUENTIAL_CHECKPOINT << " not defined..." <<
    // std::endl;
  } else {
    _checkpointFilestem = strBuf;
    // std::cout << "Initialise from sequential checkpoint using filestem " <<
    // _checkpointFilestem << std::endl;
  }

  // read molecules per direction (if no checkpoint is defined as input)
  if ((!_initFromCheckpoint) && (!_initFromSequentialCheckpoint)) {
    tarch::la::Vector<MD_DIM, int> buffer(-1);
    tarch::configuration::ParseConfiguration::readVector<MD_DIM, int>(buffer, node, MOLECULES_PER_DIRECTION);
    for (unsigned int d = 0; d < MD_DIM; d++) {
      if (buffer[d] <= 0) {
        std::cout << MOLECULES_PER_DIRECTION << ": Entry " << d << " is smaller than or equal zero: " << buffer << std::endl;
        _isValid = false;
        exit(EXIT_FAILURE);
      }
      _moleculesPerDirection[d] = (unsigned int)(buffer[d]);
    }
  }

  // read domain size
  tarch::configuration::ParseConfiguration::readVector<MD_DIM, double>(_domainSize, node, DOMAIN_SIZE);
  for (unsigned int d = 0; d < MD_DIM; d++) {
    if (_domainSize[d] <= 0.0) {
      std::cout << DOMAIN_SIZE << ": Entry " << d << " is smaller than or equal zero: " << _domainSize << std::endl;
      _isValid = false;
      exit(EXIT_FAILURE);
    }
  }

  // read domain offset
  tarch::configuration::ParseConfiguration::readVector<MD_DIM, double>(_domainOffset, node, DOMAIN_OFFSET);

  // read cutoff radius and initialise also meshwidth for linked cell scheme
  tarch::configuration::ParseConfiguration::readDoubleMandatory(_cutoffRadius, node, CUTOFF_RADIUS);
  if (_cutoffRadius <= 0.0) {
    std::cout << CUTOFF_RADIUS << " is smaller than or equal zero!" << std::endl;
    _isValid = false;
    exit(EXIT_FAILURE);
  }

  // read linked-cell size if existent; otherwise, use cut-off radius
  tarch::configuration::ParseConfiguration::readVector<MD_DIM, double>(_meshWidth, node, LINKED_CELL_SIZE);
  for (unsigned int d = 0; d < MD_DIM; d++) {
    if (_meshWidth[d] < _cutoffRadius) {
      std::cout << "ERROR: Meshwidth(" << d << ") is smaller than cut off radius!" << std::endl;
      _isValid = false;
      exit(EXIT_FAILURE);
    }
  }

  // read Boltzmann's constant
  tarch::configuration::ParseConfiguration::readDoubleMandatory(_kB, node, K_B);
  if (_kB <= 0.0) {
    std::cout << K_B << " is smaller than or equal zero!" << std::endl;
    _isValid = false;
    exit(EXIT_FAILURE);
  }

  // read block size
  tarch::configuration::ParseConfiguration::readIntMandatory(intBuf, node, BLOCK_SIZE);
  if (intBuf <= 0) {
    std::cout << BLOCK_SIZE << " is smaller than or equal zero!" << std::endl;
    _isValid = false;
    exit(EXIT_FAILURE);
  }
  _blockSize = (unsigned int)intBuf;

  // read all boundary information
  for (unsigned int i = 0; i < MD_LINKED_CELL_NEIGHBOURS; i++) {
    strBuf = "";
    tarch::configuration::ParseConfiguration::readStringMandatory(strBuf, node, BOUNDARY[i]);
    if (strBuf == PERIODIC_BOUNDARY) {
      _boundary[i] = simplemd::PERIODIC_BOUNDARY;
    } else if (strBuf == OPEN_BOUNDARY) {
      _boundary[i] = simplemd::OPEN_BOUNDARY;
    } else if (strBuf == GEOMETRY_BOUNDARY) {
      _boundary[i] = simplemd::GEOMETRY_BOUNDARY;
    } else if (strBuf == REFLECTING_BOUNDARY) {
      _boundary[i] = simplemd::REFLECTING_BOUNDARY;
    } else {
      std::cout << BOUNDARY[i] << ": Unknown boundary type " << strBuf << "!" << std::endl;
      _isValid = false;
      exit(EXIT_FAILURE);
    }
  }

  // test periodic boundary relations
  for (unsigned int i = 0; i < MD_LINKED_CELL_NEIGHBOURS; i++) {
    if (_boundary[i] == simplemd::PERIODIC_BOUNDARY) {
      if (_boundary[MD_LINKED_CELL_NEIGHBOURS - 1 - i] != simplemd::PERIODIC_BOUNDARY) {
        std::cout << "Conflict in periodic boundaries " << BOUNDARY[i] << "<->" << BOUNDARY[MD_LINKED_CELL_NEIGHBOURS - 1 - i] << std::endl;
        _isValid = false;
        exit(EXIT_FAILURE);
      }
    }
  }

#if (MD_DEBUG == MD_YES)
  std::cout << "Read from checkpoint            " << _initFromCheckpoint << ", filestem: " << _checkpointFilestem << std::endl;
  std::cout << "Read from sequential checkpoint " << _initFromSequentialCheckpoint << ", filestem: " << _checkpointFilestem << std::endl;
  std::cout << "Domain-size:                    " << _domainSize << std::endl;
  std::cout << "Domain-offset:                  " << _domainOffset << std::endl;
  std::cout << "Mol./direction:                 " << _moleculesPerDirection << std::endl;
  std::cout << "Cutoff radius:                  " << _cutoffRadius << std::endl;
  std::cout << "Linked-cell size:               " << _meshWidth << std::endl;
  std::cout << "K_B:                            " << _kB << std::endl;
  std::cout << "Boundary:                       " << _boundary << std::endl;
  std::cout << "Blocksize:                      " << _blockSize << std::endl;
#endif
}

std::string simplemd::configurations::DomainConfiguration::getTag() const { return "domain-configuration"; }

bool simplemd::configurations::DomainConfiguration::isValid() const { return _isValid; }
