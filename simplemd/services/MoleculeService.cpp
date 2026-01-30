// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#include "simplemd/services/MoleculeService.h"
#include "simplemd/molecule-mappings/WriteCheckPointMapping.h"

simplemd::services::MoleculeService::~MoleculeService() { shutdown(); }

void simplemd::services::MoleculeService::initContainer(ParallelTopologyService parallelTopologyService, size_t moleculeCount, double capacityFactor) {
#if (MD_DEBUG == MD_YES)
  if (_moleculeContainer != nullptr) {
    std::cout << "simplemd::services::MoleculeService::initContainer: _moleculeContainer has already been initialized." << std::endl;
    exit(1);
  }
#endif
  auto averageCellMoleculeCount = (double)moleculeCount / parallelTopologyService.getLocalNumberOfCellsLinear();
  unsigned int cellCapacity = averageCellMoleculeCount * capacityFactor;
  cellCapacity = cellCapacity < 5 ? 5 : cellCapacity; // let the capacity be at least 5
  _moleculeContainer = new MoleculeContainer(parallelTopologyService, cellCapacity);
#if (MD_DEBUG == MD_YES)
  if (_moleculeContainer == nullptr) {
    std::cout << "simplemd::services::MoleculeService::initContainer: _moleculeContainer could not be initialized." << std::endl;
    exit(1);
  }
#endif
}

simplemd::services::MoleculeService::MoleculeService(const tarch::la::Vector<MD_DIM, double>& localDomainSize,
                                                     const tarch::la::Vector<MD_DIM, double>& domainOffset,
                                                     const tarch::la::Vector<MD_DIM, unsigned int>& moleculesPerDirection,
                                                     const tarch::la::Vector<MD_DIM, double>& meanVelocity, const double& kB, const double& temperature,
                                                     const double capacityFactor,
                                                     const simplemd::services::MolecularPropertiesService& molecularPropertiesService,
                                                     const simplemd::services::ParallelTopologyService& parallelTopologyService)
    : _localDomainSize(localDomainSize) {
  _meanVelocity = meanVelocity;
  size_t moleculeCount = 1;
  for (int d = 0; d < MD_DIM; d++) {
    moleculeCount *= moleculesPerDirection[d];
  }
  initContainer(parallelTopologyService, moleculeCount, capacityFactor);
  tarch::la::Vector<MD_DIM, double> position(0.0);
  tarch::la::Vector<MD_DIM, double> velocity(0.0);
  unsigned int indexNumberMolecules = 0;

  Molecule tmpMolecule;
  // loop over domain and determine position vector (place molecules initially on a grid)
#if (MD_DIM > 2)
  for (unsigned int z = 0; z < moleculesPerDirection[2]; z++) {
    position[2] = (0.5 + z) * _localDomainSize[2] / moleculesPerDirection[2] + domainOffset[2];
#endif
#if (MD_DIM > 1)
    for (unsigned int y = 0; y < moleculesPerDirection[1]; y++) {
      position[1] = (0.5 + y) * _localDomainSize[1] / moleculesPerDirection[1] + domainOffset[1];
#endif
      for (unsigned int x = 0; x < moleculesPerDirection[0]; x++) {
        position[0] = (0.5 + x) * _localDomainSize[0] / moleculesPerDirection[0] + domainOffset[0];

        // get initial velocity
        getInitialVelocity(meanVelocity, kB, temperature, molecularPropertiesService, velocity);
        // initialise molecule in memory block and set ID
        tmpMolecule.setPosition(position);
        tmpMolecule.setVelocity(velocity);
        tmpMolecule.setID(indexNumberMolecules);
        _moleculeContainer->insert(tmpMolecule);
        // increment index counter
        indexNumberMolecules++;
      }
#if (MD_DIM > 1)
    }
#endif
#if (MD_DIM > 2)
  }
#endif
  // reset the mean velocity to exactly the velocity specified
  resetMeanVelocity();
}

simplemd::services::MoleculeService::MoleculeService(const std::string& checkPointFileStem, const double capacityFactor,
                                                     const simplemd::services::ParallelTopologyService& parallelTopologyService) {
  tarch::la::Vector<MD_DIM, double> position(0.0);
  tarch::la::Vector<MD_DIM, double> velocity(0.0);
  tarch::la::Vector<MD_DIM, double> force(0.0);
  unsigned int idCounter = 0;
  unsigned int dimensions = 0;
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

  unsigned int moleculeCount;
  file >> moleculeCount;
  initContainer(parallelTopologyService, moleculeCount, capacityFactor);
  file >> dimensions;
  if (dimensions != MD_DIM) {
    std::cout << "ERROR MoleculeService::MoleculeService: File " << ss.str() << " is meant for problems of dimension " << dimensions << "!" << std::endl;
    exit(EXIT_FAILURE);
  }

  idCounter = 0;
  Molecule tmpMolecule;
  for (unsigned int i = 0; i < moleculeCount; i++) {
    for (unsigned int d = 0; d < MD_DIM; d++) {
      file >> position[d];
    }
    for (unsigned int d = 0; d < MD_DIM; d++) {
      file >> velocity[d];
    }
    for (unsigned int d = 0; d < MD_DIM; d++) {
      file >> force[d];
    }

    tmpMolecule.setPosition(position);
    tmpMolecule.setVelocity(velocity);
    tmpMolecule.setForceOld(force);
    tmpMolecule.setID(idCounter);
    _moleculeContainer->insert(tmpMolecule);

    idCounter++;
  }

  // _meanVelocity was not initialised!
  // resetMeanVelocity();
}

simplemd::services::MoleculeService::MoleculeService(const tarch::la::Vector<MD_DIM, double>& localDomainSize,
                                                     const tarch::la::Vector<MD_DIM, double>& localDomainOffset, const std::string& checkPointFileStem,
                                                     const double capacityFactor, const simplemd::services::ParallelTopologyService& parallelTopologyService)
    : _localDomainSize(localDomainSize) {
  // buffer for single molecule
  tarch::la::Vector<MD_DIM, double> position(0.0);
  tarch::la::Vector<MD_DIM, double> velocity(0.0);
  tarch::la::Vector<MD_DIM, double> force(0.0);
  unsigned int idCounter = 0;
  unsigned int dimensions = 0;
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

  unsigned int moleculeCount;
  file >> moleculeCount;
  initContainer(parallelTopologyService, moleculeCount, capacityFactor);
  file >> dimensions;
  if (dimensions != MD_DIM) {
    std::cout << "ERROR MoleculeService::MoleculeService: File " << ss.str() << " is meant for problems of dimension " << dimensions << "!" << std::endl;
    exit(EXIT_FAILURE);
  }

  for (unsigned int i = 0; i < moleculeCount; i++) {
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
      isInside = isInside && (position[d] >= localDomainOffset[d]) && (position[d] < localDomainOffset[d] + _localDomainSize[d]);
    }
    if (isInside) {
      positions.push_back(position);
      velocities.push_back(velocity);
      forces.push_back(force);
    }
  }
  file.close();
  // update number of molecules
  moleculeCount = (unsigned int)positions.size();
  if (moleculeCount < 1) {
    std::cout << "ERROR MoleculeService::MoleculeService: moleculeCount < 1!" << std::endl;
    exit(EXIT_FAILURE);
  }

  idCounter = 0;
  Molecule tmpMolecule;
  // write molecule data from dynamic vectors to allocated memory
  for (unsigned int i = 0; i < moleculeCount; i++) {
    position = positions[i];
    velocity = velocities[i];
    force = forces[i];

    tmpMolecule.setPosition(position);
    tmpMolecule.setVelocity(velocity);
    tmpMolecule.setForceOld(force);
    tmpMolecule.setID(idCounter);
    _moleculeContainer->insert(tmpMolecule);
    idCounter++;
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
  initialVelocity[0] = meanVelocity[0] + stdDeviation * (randomNumbers[0] * TARCH_COS(randomNumbers[1]));
  initialVelocity[1] = meanVelocity[1] + stdDeviation * (randomNumbers[0] * TARCH_SIN(randomNumbers[1]));
#elif (MD_DIM == 3)
  initialVelocity[0] = meanVelocity[0] + stdDeviation * (randomNumbers[0] * TARCH_SIN(randomNumbers[1]) * TARCH_COS(randomNumbers[2]));
  initialVelocity[1] = meanVelocity[1] + stdDeviation * (randomNumbers[0] * TARCH_SIN(randomNumbers[1]) * TARCH_SIN(randomNumbers[2]));
  initialVelocity[2] = meanVelocity[2] + stdDeviation * (randomNumbers[0] * TARCH_COS(randomNumbers[1]));
#endif
}

void simplemd::services::MoleculeService::shutdown() {
  if (_moleculeContainer != nullptr) {
    delete _moleculeContainer;
    _moleculeContainer = nullptr;
  }
}

unsigned int simplemd::services::MoleculeService::getLocalNumberOfMoleculesWithGhost() const {
  return _moleculeContainer->getLocalNumberOfMoleculesWithGhost();
}

void simplemd::services::MoleculeService::writeCheckPoint(const simplemd::services::ParallelTopologyService& parallelTopologyService,
                                                          const std::string& filestem, const unsigned int& t) {
  simplemd::moleculemappings::WriteCheckPointMapping writeCheckPointMapping(parallelTopologyService, filestem, t);
  _moleculeContainer->iterateMolecules(writeCheckPointMapping);
}

void simplemd::services::MoleculeService::resetMeanVelocity() {
  simplemd::moleculemappings::ComputeMeanVelocityMapping compute;
  _moleculeContainer->iterateMolecules(compute);
  tarch::la::Vector<MD_DIM, double> currentVel = compute.getMeanVelocity();
  simplemd::moleculemappings::SetMeanVelocityMapping set(currentVel, _meanVelocity);
  _moleculeContainer->iterateMolecules(set);

  // check again
  _moleculeContainer->iterateMolecules(compute);
}

bool simplemd::services::MoleculeService::tarchDebugIsOn() {
#if (TARCH_DEBUG == TARCH_YES)
  return true;
#else
  return false;
#endif
}
