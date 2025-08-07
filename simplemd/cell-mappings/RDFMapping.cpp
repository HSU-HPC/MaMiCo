// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#include "simplemd/cell-mappings/RDFMapping.h"

#include "simplemd/MolecularDynamicsDefinitions.h"
#include <fstream>

simplemd::cellmappings::RDFMapping::RDFMapping(const simplemd::services::ParallelTopologyService& parallelTopologyService,
                                               simplemd::services::LinkedCellService& linkedCellService, const double& cutoffRadius,
                                               const unsigned int& numberIntervals)
    : _parallelTopologyService(parallelTopologyService), _linkedCellService(linkedCellService), _cutoffRadius(cutoffRadius), _numberIntervals(numberIntervals),
      _meshsize(cutoffRadius / ((double)numberIntervals)), _particleCounter(0.0), _evaluationCounter(0.0) {
  _particlesPerInterval.clear();
  for (unsigned int i = 0; i < numberIntervals; i++) {
    _particlesPerInterval.push_back(0.0);
  }
}

simplemd::cellmappings::RDFMapping::~RDFMapping() {}

void simplemd::cellmappings::RDFMapping::beginCellIteration() {
  // increment evaluation counter
  _evaluationCounter += 1.0;
}

void simplemd::cellmappings::RDFMapping::endCellIteration() {}

void simplemd::cellmappings::RDFMapping::evaluateRDF(const unsigned int& localMDSimulation) {
  std::stringstream ss;
  ss << "RDF_" << localMDSimulation << "_";
#if (MD_PARALLEL == MD_YES)
  ss << "_" << _parallelTopologyService.getRank();
#endif
  ss << "_" << (int)_evaluationCounter;
  ss << ".dat";
  std::ofstream file(ss.str().c_str());
  if (!file.is_open()) {
    std::cout << "ERROR RDFMapping: Could not open file!" << std::endl;
    exit(EXIT_FAILURE);
  }
  unsigned int startWritingPosition = 0;

  // compute avg. number density so far
  double numberDensity = 1.0;
  for (unsigned int d = 0; d < MD_DIM; d++) {
    numberDensity = numberDensity * _linkedCellService.getLocalDomainSize()[d];
  }
  numberDensity = _particleCounter / (_evaluationCounter * numberDensity);

  // find starting position for writing
  for (unsigned int i = 0; i < _numberIntervals - 1; i++) {
    if ((_particlesPerInterval[i] == 0) && (_particlesPerInterval[i + 1] == 0)) {
      startWritingPosition = i + 1;
    } else {
      break;
    }
  }

  // write number of intervals and meshsize
  file << _numberIntervals - startWritingPosition << " " << _meshsize << std::endl;
  for (unsigned int i = startWritingPosition; i < _numberIntervals; i++) {
    const double volume
#if (MD_DIM == 1)
        = _meshsize
#elif (MD_DIM == 2)
        = MD_PI * ((i + 1) * (i + 1) * _meshsize * _meshsize - i * i * _meshsize * _meshsize)
#elif (MD_DIM == 3)
        = MD_PI * 4.0 / 3.0 * ((i + 1) * (i + 1) * (i + 1) * _meshsize * _meshsize * _meshsize - i * i * i * _meshsize * _meshsize * _meshsize)
#endif
        ;
    double buf = _particlesPerInterval[i] / (volume * numberDensity * _particleCounter);
    file << i * _meshsize + 0.5 * _meshsize << " " << buf << std::endl;
  }
  file.close();
}

void simplemd::cellmappings::RDFMapping::handleCell(LinkedCell& cell, const unsigned int& cellIndex) {
  // only consider non-ghost cells
  if (cell.isGhostCell()) {
    return;
  }

  double dist;
  unsigned int interval;

  for (auto m1 = cell.begin(); m1 != cell.end()--; m1++) {
    auto m2 = m1;

    // iterate over all other molecules not touched so far
    m2++;
    while (m2 != cell.end()) {

      dist = std::sqrt(tarch::la::dot((m2->getConstPosition() - m1->getConstPosition()), (m2->getConstPosition() - m1->getConstPosition())));

      if (dist < _cutoffRadius) {
        interval = (unsigned int)(dist / _meshsize);
        _particlesPerInterval[interval] += 2.0;
      }

      m2++;
    }
    _particleCounter += 1.0;
  }
}
void simplemd::cellmappings::RDFMapping::handleCellPair(LinkedCell& cell1, LinkedCell& cell2, const unsigned int& cellIndex1, const unsigned int& cellIndex2) {
  double dist;
  unsigned int interval;
  for (auto m1 = cell1.begin(); m1 != cell1.end(); m1++) {

    for (auto m2 = cell2.begin(); m2 != cell2.end(); m2++) {
      dist = std::sqrt(tarch::la::dot((m2->getConstPosition() - m1->getConstPosition()), (m2->getConstPosition() - m1->getConstPosition())));
      interval = (unsigned int)(dist / _meshsize);
      if (dist < _cutoffRadius) {
        interval = (unsigned int)(dist / _meshsize);
        _particlesPerInterval[interval] += ((double)(!cell1.isGhostCell()) + (double)(!cell2.isGhostCell()));
      }
    }
  }
}
