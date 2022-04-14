// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#include "tarch/utils/RandomNumberService.h"
#include <ctime>
#if (TARCH_PARALLEL == TARCH_YES)
#include <mpi.h>
#endif

tarch::utils::RandomNumberService& tarch::utils::RandomNumberService::getInstance() {
  static tarch::utils::RandomNumberService singleton;
  return singleton;
}

void tarch::utils::RandomNumberService::init(bool fixSeed) {
  if (_isInitialized)
    return;
  // for testing purpose: fix seed to a certain number
  if (fixSeed) {
    srand(10000);
  } else {
    int rank = 0;
#if (TARCH_PARALLEL == TARCH_YES)
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
    srand(static_cast<unsigned int>(time(NULL)) * (rank + 1));
  }
  _randomNumbers.assign(0.0);
  _isFirstRandomNumber = true;
  _isInitialized = true;
}

void tarch::utils::RandomNumberService::shutdown() {}

double tarch::utils::RandomNumberService::getUniformRandomNumber() const {
  // We do not want this method to return 1.0, as this can lead to invalid behaviour
  // in some situations. E.g., the Usher particle insertion relies on this method
  // in order to generate valid positions within the local subdomain, which only
  // the lower edge is part of, but not the upper edge.
  return (0.5 * (double)rand()) / (0.5 * (double)RAND_MAX + 1.0);
}

double tarch::utils::RandomNumberService::getGaussianRandomNumber() {

  // if the first random number is going to be called, create Gaussian random numbers
  if (_isFirstRandomNumber) {

    // in this case: generate new numbers
    double s = 2.0;
    tarch::la::Vector<2, double> v(0.0);

    while (s >= 1.0) {
      _randomNumbers[0] = ((double)rand()) / ((double)RAND_MAX);
      _randomNumbers[1] = ((double)rand()) / ((double)RAND_MAX);
      v[0] = (2.0 * _randomNumbers[0] - 1.0);
      v[1] = (2.0 * _randomNumbers[1] - 1.0);
      s = v[0] * v[0] + v[1] * v[1];
    }
    _randomNumbers[0] = v[0] * sqrt((-2.0 * log(s)) / s);
    _randomNumbers[1] = v[1] * sqrt((-2.0 * log(s)) / s);

    // change to the other variable in the next call
    _isFirstRandomNumber = false;
    return _randomNumbers[0];

    // otherwise: change to the other random number
  } else {
    _isFirstRandomNumber = true;
    return _randomNumbers[1];
  }
}
