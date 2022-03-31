// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_CELLMAPPINGS_RDFMAPPING_H_
#define _MOLECULARDYNAMICS_CELLMAPPINGS_RDFMAPPING_H_

#include "simplemd/LinkedCell.h"
#include "simplemd/services/MolecularPropertiesService.h"
#include "simplemd/services/LinkedCellService.h"
#include "simplemd/services/ParallelTopologyService.h"
#include <iostream>
#include <sstream>

namespace simplemd {
namespace cellmappings { class RDFMapping; }
}

/** computes the radial distribution function for the molecules.
 *  Should be called using the iterateCellPairs()-method from the
 * LinkedCellService.
 *
 *  @author Philipp Neumann
 */
class simplemd::cellmappings::RDFMapping {
public:
  RDFMapping(const simplemd::services::ParallelTopologyService &
                 parallelTopologyService,
             simplemd::services::LinkedCellService &linkedCellService,
             const double &cutoffRadius, const unsigned int &numberIntervals);

  ~RDFMapping();

  void beginCellIteration();

  void endCellIteration();

  void evaluateRDF(const unsigned int &localMDSimulation);

  void handleCell(LinkedCell &cell, const unsigned int &cellIndex);
  void handleCellPair(LinkedCell &cell1, LinkedCell &cell2,
                      const unsigned int &cellIndex1,
                      const unsigned int &cellIndex2);

private:
  const simplemd::services::ParallelTopologyService &_parallelTopologyService;
  simplemd::services::LinkedCellService &_linkedCellService;
  // cut-off radius
  const double _cutoffRadius;
  // number of shells that are considered
  const unsigned int _numberIntervals;
  // width of each shell
  const double _meshsize;
  // contains the total number of particles in each shell
  std::vector<double> _particlesPerInterval;
  // contains the total number of particles
  double _particleCounter;

  // counts the evaluations of this mapping
  double _evaluationCounter;
};

#endif // _MOLECULARDYNAMICS_CELLMAPPINGS_EMPTYMAPPING_H_
