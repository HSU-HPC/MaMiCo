// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#include "simplemd/services/MoleculeService.h"
#include "simplemd/cell-mappings/CopyMoleculesMapping.h"
#include "simplemd/molecule-mappings/UpdateLinkedCellListsMapping.h"
#include "simplemd/molecule-mappings/WriteCheckPointMapping.h"
#include "simplemd/services/LinkedCellService.h"

simplemd::services::MoleculeService::~MoleculeService() {
  for (unsigned int i = 0; i < _molecules.size(); i++) {
    if (_molecules[i] != NULL) {
      free(_molecules[i]);
      _molecules[i] = NULL;
    }
  }
  _molecules.clear();
  _freeMoleculePositions.clear();
  _numberMolecules = 0;
  _blockSize = 0;
}

simplemd::services::MoleculeService::MoleculeService(const tarch::la::Vector<MD_DIM, double>& domainSize, const tarch::la::Vector<MD_DIM, double>& domainOffset,
                                                     const tarch::la::Vector<MD_DIM, unsigned int>& moleculesPerDirection,
                                                     const tarch::la::Vector<MD_DIM, double>& meanVelocity, const double& kB, const double& temperature,
                                                     const unsigned int& blockSize,
                                                     const simplemd::services::MolecularPropertiesService& molecularPropertiesService) {
  tarch::la::Vector<MD_DIM, double> position(0.0);
  tarch::la::Vector<MD_DIM, double> velocity(0.0);
  unsigned int indexNumberMolecules = 0;
  unsigned int numberBlocks = 0;

  // delete old stuff and set new variables
  for (unsigned int i = 0; i < _molecules.size(); i++) {
    if (_molecules[i] != NULL) {
      free(_molecules[i]);
      _molecules[i] = NULL;
    }
  }
  _molecules.clear();
  _freeMoleculePositions.clear();
  _blockSize = blockSize;
  _meanVelocity = meanVelocity;
  _numberMolecules = moleculesPerDirection[0];
  for (unsigned int d = 1; d < MD_DIM; d++) {
    _numberMolecules = _numberMolecules * moleculesPerDirection[d];
  }
  numberBlocks = _numberMolecules / _blockSize + (_numberMolecules % _blockSize != 0);
  if (_numberMolecules < 1) {
    std::cout << "ERROR simplemd::services::MoleculeService::MoleculeService(): _numberMolecules < 1!" << std::endl;
    exit(EXIT_FAILURE);
  }

  // allocate memory and initialise it
  for (unsigned int i = 0; i < numberBlocks; i++) {
    _molecules.push_back((Molecule*)NULL);
    _molecules[i] = (Molecule*)malloc(sizeof(Molecule) * _blockSize);
    if (_molecules[i] == NULL) {
      std::cout << "ERROR simplemd::services::MoleculeService::MoleculeService(): _molecules == NULL!" << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  // loop over domain and determine position vector (place molecules initially on a grid)
#if (MD_DIM > 2)
  for (unsigned int z = 0; z < moleculesPerDirection[2]; z++) {
    position[2] = (0.5 + z) * domainSize[2] / moleculesPerDirection[2] + domainOffset[2];
#endif
#if (MD_DIM > 1)
    for (unsigned int y = 0; y < moleculesPerDirection[1]; y++) {
      position[1] = (0.5 + y) * domainSize[1] / moleculesPerDirection[1] + domainOffset[1];
#endif
      for (unsigned int x = 0; x < moleculesPerDirection[0]; x++) {
        position[0] = (0.5 + x) * domainSize[0] / moleculesPerDirection[0] + domainOffset[0];

        // get initial velocity
        getInitialVelocity(meanVelocity, kB, temperature, molecularPropertiesService, velocity);
        // initialise molecule in memory block and set ID
        Molecule* myNewMolecule = NULL;
        unsigned int blockId = indexNumberMolecules / _blockSize;
        unsigned int blockIndex = indexNumberMolecules % _blockSize;
        myNewMolecule = new (&_molecules[blockId][blockIndex]) Molecule(position, velocity);
#if (MD_DEBUG == MD_YES)
        if (myNewMolecule == NULL) {
          std::cout << "ERROR simplemd::services::MoleculeService::MoleculeService(): myNewMolecule==NULL!" << std::endl;
          exit(EXIT_FAILURE);
        }
#endif
        myNewMolecule->setID(indexNumberMolecules);
        // increment index counter
        indexNumberMolecules++;
      }
#if (MD_DIM > 1)
    }
#endif
#if (MD_DIM > 2)
  }
#endif

  // determine the number of free molecule positions in last block and set
  // number of molecules to the real value
  _numberMolecules = indexNumberMolecules;
  for (unsigned int i = indexNumberMolecules; i < _blockSize * _molecules.size(); i++) {
    _freeMoleculePositions.push_back(i);
  }

  // reset the mean velocity to exactly the velocity specified
  resetMeanVelocity();
}

simplemd::services::MoleculeService::MoleculeService(const tarch::la::Vector<MD_DIM, double>& domainSize, const tarch::la::Vector<MD_DIM, double>& domainOffset,
                                                     const std::string& checkPointFileStem, const unsigned int& blockSize,
                                                     const simplemd::services::ParallelTopologyService& parallelTopologyService) {
  tarch::la::Vector<MD_DIM, double> position(0.0);
  tarch::la::Vector<MD_DIM, double> velocity(0.0);
  tarch::la::Vector<MD_DIM, double> force(0.0);
  unsigned int idCounter = 0;
  unsigned int numberBlocks = 0;
  std::stringstream ss;
  ss << checkPointFileStem;
#if (MD_PARALLEL == MD_YES)
  ss << "_" << parallelTopologyService.getRank();
#else
  ss << "_0";
#endif
  ss << ".checkpoint";
  std::ifstream file(ss.str().c_str());
  if (!file.is_open()) {
    std::cout << "ERROR MoleculeService::MoleculeService: Could not open file " << ss.str() << "!" << std::endl;
    exit(EXIT_FAILURE);
  }

  // delete old stuff and set new variables
  for (unsigned int i = 0; i < _molecules.size(); i++) {
    if (_molecules[i] != NULL) {
      free(_molecules[i]);
      _molecules[i] = NULL;
    }
  }
  _molecules.clear();
  _freeMoleculePositions.clear();
  _blockSize = blockSize;
  file >> _numberMolecules;
  file >> numberBlocks;
  if (numberBlocks != MD_DIM) {
    std::cout << "ERROR MoleculeService::MoleculeService: File " << ss.str() << " is meant for problems of dimension " << numberBlocks << "!" << std::endl;
    exit(EXIT_FAILURE);
  }
  numberBlocks = _numberMolecules / _blockSize + (_numberMolecules % _blockSize != 0);
  if (_numberMolecules < 1) {
    std::cout << "ERROR MoleculeService::MoleculeService: _numberMolecules < 1!" << std::endl;
    exit(EXIT_FAILURE);
  }

  // allocate memory and initialise it
  for (unsigned int i = 0; i < numberBlocks; i++) {
    _molecules.push_back((Molecule*)NULL);
    _molecules[i] = (Molecule*)malloc(sizeof(Molecule) * _blockSize);
    if (_molecules[i] == NULL) {
      std::cout << "ERROR simplemd::services::MoleculeService::MoleculeService(): _molecules == NULL!" << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  idCounter = 0;
  for (unsigned int i = 0; i < _numberMolecules; i++) {
    for (unsigned int d = 0; d < MD_DIM; d++) {
      file >> position[d];
    }
    for (unsigned int d = 0; d < MD_DIM; d++) {
      file >> velocity[d];
    }
    for (unsigned int d = 0; d < MD_DIM; d++) {
      file >> force[d];
    }

    Molecule* myNewMolecule = NULL;
    unsigned int blockId = idCounter / _blockSize;
    unsigned int blockIndex = idCounter % _blockSize;
    myNewMolecule = new (&_molecules[blockId][blockIndex]) Molecule(position, velocity);
#if (MD_DEBUG == MD_YES)
    if (myNewMolecule == NULL) {
      std::cout << "ERROR simplemd::services::MoleculeService::MoleculeService(): myNewMolecule==NULL!" << std::endl;
      exit(EXIT_FAILURE);
    }
#endif
    myNewMolecule->setForceOld(force);
    myNewMolecule->setID(idCounter);
    idCounter++;
  }
  file.close();
  // add free molecule positions
  for (unsigned int i = idCounter; i < _blockSize * _molecules.size(); i++) {
    _freeMoleculePositions.push_back(i);
  }

  // _meanVelocity was not initialised!
  // resetMeanVelocity();
}

simplemd::services::MoleculeService::MoleculeService(const tarch::la::Vector<MD_DIM, double>& localDomainSize,
                                                     const tarch::la::Vector<MD_DIM, double>& localDomainOffset, const std::string& checkPointFileStem,
                                                     const unsigned int& blockSize) {
  // buffer for single molecule
  tarch::la::Vector<MD_DIM, double> position(0.0);
  tarch::la::Vector<MD_DIM, double> velocity(0.0);
  tarch::la::Vector<MD_DIM, double> force(0.0);
  unsigned int idCounter = 0;
  unsigned int numberBlocks = 0;
  // buffer for all relevant molecules
  std::vector<tarch::la::Vector<MD_DIM, double>> positions;
  std::vector<tarch::la::Vector<MD_DIM, double>> velocities;
  std::vector<tarch::la::Vector<MD_DIM, double>> forces;
  std::stringstream ss;
  // use information of rank 0
  ss << checkPointFileStem << "_0.checkpoint";
  std::ifstream file(ss.str().c_str());
  if (!file.is_open()) {
    std::cout << "ERROR MoleculeService::MoleculeService: Could not open file " << ss.str() << "!" << std::endl;
    exit(EXIT_FAILURE);
  }

  // delete old stuff and set new variables
  for (unsigned int i = 0; i < _molecules.size(); i++) {
    if (_molecules[i] != NULL) {
      free(_molecules[i]);
      _molecules[i] = NULL;
    }
  }
  _molecules.clear();
  _freeMoleculePositions.clear();
  _blockSize = blockSize;

  file >> _numberMolecules;
  file >> numberBlocks;
  if (numberBlocks != MD_DIM) {
    std::cout << "ERROR MoleculeService::MoleculeService: File " << ss.str() << " is meant for problems of dimension " << numberBlocks << "!" << std::endl;
    exit(EXIT_FAILURE);
  }
  for (unsigned int i = 0; i < _numberMolecules; i++) {
    for (unsigned int d = 0; d < MD_DIM; d++) {
      file >> position[d];
    }
    for (unsigned int d = 0; d < MD_DIM; d++) {
      file >> velocity[d];
    }
    for (unsigned int d = 0; d < MD_DIM; d++) {
      file >> force[d];
    }
    // check if molecule is inside this local domain and if so: add it to
    // dynamic vectors
    bool isInside = true;
    for (unsigned d = 0; d < MD_DIM; d++) {
      isInside = isInside && (position[d] >= localDomainOffset[d]) && (position[d] < localDomainOffset[d] + localDomainSize[d]);
    }
    if (isInside) {
      positions.push_back(position);
      velocities.push_back(velocity);
      forces.push_back(force);
    }
  }
  file.close();
  // update number of molecules
  _numberMolecules = (unsigned int)positions.size();
  numberBlocks = _numberMolecules / _blockSize + (_numberMolecules % _blockSize != 0);
  if (_numberMolecules < 1) {
    std::cout << "ERROR MoleculeService::MoleculeService: _numberMolecules < 1!" << std::endl;
    exit(EXIT_FAILURE);
  }

  // allocate memory and initialise it
  for (unsigned int i = 0; i < numberBlocks; i++) {
    _molecules.push_back((Molecule*)NULL);
    _molecules[i] = (Molecule*)malloc(sizeof(Molecule) * _blockSize);
    if (_molecules[i] == NULL) {
      std::cout << "ERROR simplemd::services::MoleculeService::MoleculeService(): _molecules == NULL!" << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  idCounter = 0;
  // write molecule data from dynamic vectors to allocated memory
  for (unsigned int i = 0; i < _numberMolecules; i++) {
    position = positions[i];
    velocity = velocities[i];
    force = forces[i];

    Molecule* myNewMolecule = NULL;
    unsigned int blockId = idCounter / _blockSize;
    unsigned int blockIndex = idCounter % _blockSize;
    myNewMolecule = new (&_molecules[blockId][blockIndex]) Molecule(position, velocity);
#if (MD_DEBUG == MD_YES)
    if (myNewMolecule == NULL) {
      std::cout << "ERROR simplemd::services::MoleculeService::MoleculeService(): myNewMolecule==NULL!" << std::endl;
      exit(EXIT_FAILURE);
    }
#endif
    myNewMolecule->setForceOld(force);
    myNewMolecule->setID(idCounter);
    idCounter++;
  }
  // add free molecule positions
  for (unsigned int i = idCounter; i < _blockSize * _molecules.size(); i++) {
    _freeMoleculePositions.push_back(i);
  }

  // empty vectors
  positions.clear();
  velocities.clear();
  forces.clear();

  // _meanVelocity was not initialised!
  // resetMeanVelocity();
}

void simplemd::services::MoleculeService::getInitialVelocity(const tarch::la::Vector<MD_DIM, double>& meanVelocity, const double& kB, const double& temperature,
                                                             const simplemd::services::MolecularPropertiesService& molecularPropertiesService,
                                                             tarch::la::Vector<MD_DIM, double>& initialVelocity) const {
  tarch::la::Vector<MD_DIM, double> randomNumbers(0.0);
  // standard deviation for fluctuation of velocity around meanVelocity
  double stdDeviation = std::sqrt(MD_DIM * kB * temperature / molecularPropertiesService.getMolecularProperties().getMass());

  // random number (unit mean
  randomNumbers[0] = tarch::utils::RandomNumberService::getInstance().getGaussianRandomNumber();
  // then, get D-1 uniform random numbers over 2PI
  for (unsigned int d = 1; d < MD_DIM; d++) {
    randomNumbers[d] = 2.0 * MD_PI * tarch::utils::RandomNumberService::getInstance().getUniformRandomNumber();
  }

  // put initial velocity together. For 2D/ 3D, we use the polar/ spherical coordinates to generate the fluctuation around the mean velocity
#if (MD_DIM == 1)
  initialVelocity = meanVelocity + stdDeviation * randomNumbers;
#elif (MD_DIM == 2)
  initialVelocity[0] = meanVelocity[0] + stdDeviation * (randomNumbers[0] * std::cos(randomNumbers[1]));
  initialVelocity[1] = meanVelocity[1] + stdDeviation * (randomNumbers[0] * std::sin(randomNumbers[1]));
#elif (MD_DIM == 3)
  initialVelocity[0] = meanVelocity[0] + stdDeviation * (randomNumbers[0] * std::sin(randomNumbers[1]) * std::cos(randomNumbers[2]));
  initialVelocity[1] = meanVelocity[1] + stdDeviation * (randomNumbers[0] * std::sin(randomNumbers[1]) * std::sin(randomNumbers[2]));
  initialVelocity[2] = meanVelocity[2] + stdDeviation * (randomNumbers[0] * std::cos(randomNumbers[1]));
#endif
}

void simplemd::services::MoleculeService::shutdown() {
  for (unsigned int i = 0; i < _molecules.size(); i++) {
    if (_molecules[i] != NULL) {
      free(_molecules[i]);
      _molecules[i] = NULL;
    }
  }
  _molecules.clear();
  _freeMoleculePositions.clear();
  _numberMolecules = 0;
  _blockSize = 0;
}

simplemd::Molecule* simplemd::services::MoleculeService::addMolecule(const Molecule& molecule) {
  _numberMolecules++;

  // if there is a free position within the available memory, just put molecule to the last entry
  // and erase that entry afterwards. Important: set ID of the new molecule to its current position.
  if (!_freeMoleculePositions.empty()) {
    // determine block and index within block
    const unsigned int front = _freeMoleculePositions.front();
    const unsigned int blockId = front / _blockSize;
    const unsigned int blockIndex = front % _blockSize;
    _molecules[blockId][blockIndex] = molecule;
    _molecules[blockId][blockIndex].setID(front);
    _freeMoleculePositions.pop_front();
    return &_molecules[blockId][blockIndex];

    // otherwise: reallocate memory and adapt _numberMolecules and _freeMoleculePositions respectively. Besides,
    // set ID of the added molecule.
  } else {
    _molecules.push_back((Molecule*)NULL);
    _molecules[_molecules.size() - 1] = (Molecule*)malloc(sizeof(Molecule) * _blockSize);
    if (_molecules[_molecules.size() - 1] == NULL) {
      std::cout << "ERROR simplemd::services::MoleculeService::addMolecule(): malloc!" << std::endl;
      exit(EXIT_FAILURE);
    }
    const unsigned int blockSize = _blockSize;
    const unsigned int size = (unsigned int)(_molecules.size() - 1);
    const unsigned int id = size * blockSize;
    _molecules[size][0] = molecule;
    _molecules[size][0].setID(id);
    for (unsigned int i = 1; i < blockSize; i++) {
      _freeMoleculePositions.push_back(i + id);
    }
    return &_molecules[size][0];
  }
}

const unsigned int& simplemd::services::MoleculeService::getNumberMolecules() const { return _numberMolecules; }

void simplemd::services::MoleculeService::deleteMolecule(Molecule& molecule) {
  // add free position at the front
  _freeMoleculePositions.push_front(molecule.getID());
  _numberMolecules--;
}

void simplemd::services::MoleculeService::reorganiseMemory(const simplemd::services::ParallelTopologyService& parallelTopologyService,
                                                           simplemd::services::LinkedCellService& linkedCellService) {
  unsigned int numberBlocks = 0;
  // loop counter
  unsigned int i = 0;
  simplemd::cellmappings::CopyMoleculesMapping copyMoleculesMapping(*this);
  simplemd::moleculemappings::UpdateLinkedCellListsMapping updateLinkedCellListsMapping(parallelTopologyService, linkedCellService);

  // create copy of molecules and store the copies in a std::list
  std::list<simplemd::Molecule>& copyBuffer = copyMoleculesMapping.getCopyOfMolecules();
  linkedCellService.iterateCells(copyMoleculesMapping, false);
  // std::cout <<  "Number molecules: " << copyBuffer.size() << std::endl;

  // reset all values and delete molecules
  _numberMolecules = (unsigned int)(copyBuffer.size());
  _freeMoleculePositions.clear();
  for (i = 0; i < _molecules.size(); i++) {
    if (_molecules[i] != NULL) {
      free(_molecules[i]);
      _molecules[i] = NULL;
    }
  }
  _molecules.clear();

  // create new memory
  numberBlocks = _numberMolecules / _blockSize + (_numberMolecules % _blockSize != 0);
  for (i = 0; i < numberBlocks; i++) {
    _molecules.push_back((Molecule*)NULL);
    _molecules[i] = (Molecule*)malloc(sizeof(Molecule) * _blockSize);
    if (_molecules[i] == NULL) {
      std::cout << "ERROR simplemd::services::MoleculeService::reorganiseMemory(): _molecules == NULL!" << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  // copy molecules into correct positions and delete temporary copies
  i = 0;
  for (std::list<simplemd::Molecule>::iterator it = copyBuffer.begin(); it != copyBuffer.end(); it++) {
    const unsigned int block = i / _blockSize;
    const unsigned int blockEntry = i % _blockSize;
    // initialise object structure and memory
    simplemd::Molecule* myMolecule = new (&_molecules[block][blockEntry]) Molecule();
    *myMolecule = (*it);
    // set current molecule ID
    myMolecule->setID(i);

    // increment counter
    i++;
  }
  copyMoleculesMapping.removeCopy();

  // add empty entries to list
  for (i = _numberMolecules; i < _blockSize * _molecules.size(); i++) {
    _freeMoleculePositions.push_back(i);
  }

  // populate linked cell lists again
  iterateMolecules(updateLinkedCellListsMapping, false);
}

void simplemd::services::MoleculeService::writeCheckPoint(const simplemd::services::ParallelTopologyService& parallelTopologyService,
                                                          const std::string& filestem, const unsigned int& t) {
  simplemd::moleculemappings::WriteCheckPointMapping writeCheckPointMapping(parallelTopologyService, filestem, t);
  iterateMolecules(writeCheckPointMapping, false);
}

void simplemd::services::MoleculeService::resetMeanVelocity() {
  simplemd::moleculemappings::ComputeMeanVelocityMapping compute;
  iterateMolecules(compute, false);
  tarch::la::Vector<MD_DIM, double> currentVel = compute.getMeanVelocity();
  simplemd::moleculemappings::SetMeanVelocityMapping set(currentVel, _meanVelocity);
  iterateMolecules(set, false);

  // check again
  iterateMolecules(compute, false);
}
