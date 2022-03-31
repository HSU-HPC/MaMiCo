// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef COUPLING_INTERFACE_IMPL_MARDYN_LINKEDCELLSFORCOUPLING_H_
#define COUPLING_INTERFACE_IMPL_MARDYN_LINKEDCELLSFORCOUPLING_H_

#include "tarch/la/Vector.h"

#include "particleContainer/LinkedCells.h"
#include "particleContainer/ParticleCell.h"
#include "particleContainer/ParticleContainer.h"
#include "utils/Logger.h"
using Log::global_log;

/*
 * 	Helper class extending MarDyn's LinkedCells
 * 	Gives a Constructor that provides the creation of cells based on a given
 * cell size
 * 	Other LinkedCell parameters are set according to this cell size
 * 	@author Hanno Flohr
 */
class LinkedCellsForCoupling : public LinkedCells {

public:
  LinkedCellsForCoupling(tarch::la::Vector<3, double> mamicoCellSize, tarch::la::Vector<3, unsigned int> linkedCellsPerMacroscopicCell, double bBoxMin[3],
                         double bBoxMax[3], double cutoffRadius) {
    _cutoffRadius = cutoffRadius;
    _LJCutoffRadius = cutoffRadius;

    for (int d = 0; d < 3; d++) {
      _boundingBoxMin[d] = bBoxMin[d];
      _boundingBoxMax[d] = bBoxMax[d];
    }

    global_log->debug() << "cutoff: " << cutoffRadius << std::endl;
    global_log->debug() << "LJ cutoff" << _LJCutoffRadius << std::endl;

    int cellCount = 1;
    for (int d = 0; d < 3; d++) {
      // calculate cell length
      _cellLength[d] = mamicoCellSize[d] / linkedCellsPerMacroscopicCell[d];

      // set number of cells
      _boxWidthInNumCells[d] = (_boundingBoxMax[d] - _boundingBoxMin[d]) / _cellLength[d];
      if (_boxWidthInNumCells[d] == 0)
        _boxWidthInNumCells[d] = 1;

      // set halo parameters
      _haloWidthInNumCells[d] = ceil(cutoffRadius / _cellLength[d]);
      _haloLength[d] = _haloWidthInNumCells[d] * _cellLength[d];
      _haloBoundingBoxMin[d] = _boundingBoxMin[d] - _haloLength[d];
      _haloBoundingBoxMax[d] = _boundingBoxMax[d] + _haloLength[d];

      _cellsPerDimension[d] = _boxWidthInNumCells[d] + 2 * _haloWidthInNumCells[d];

      cellCount *= _cellsPerDimension[d];
      assert(cellCount > 0);
    }
    global_log->debug() << "Cell size (" << _cellLength[0] << ", " << _cellLength[1] << ", " << _cellLength[2] << ")" << std::endl;
    global_log->debug() << "Cells per dimension (incl. halo): " << _cellsPerDimension[0] << ", " << _cellsPerDimension[1] << ", " << _cellsPerDimension[2]
                        << std::endl;

    _cellsInCutoff = floor(_cellLength[0] / cutoffRadius);
    if (_cellsInCutoff < 1) {
      global_log->error() << "LinkedCellsForCoupling (constructor): cellsInCutOff < 1! cell "
                             "length should be >= cutoff radius"
                          << std::endl;
      exit(1);
    }

    _cells.resize(cellCount);

    // If the width of the inner region is less than the width of the halo
    // region a parallelization is not possible (with the used algorithms).
    // If a particle leaves the box, it would need to be communicated to the two
    // next neighbors.
    if (_boxWidthInNumCells[0] < 2 * _haloWidthInNumCells[0] || _boxWidthInNumCells[1] < 2 * _haloWidthInNumCells[1] ||
        _boxWidthInNumCells[2] < 2 * _haloWidthInNumCells[2]) {
      global_log->error() << "LinkedCellsForCoupling (constructor): width of "
                             "inner region is less than the width of halo."
                          << std::endl;
      global_log->error() << "_cellsPerDimension: " << _cellsPerDimension[0] << " / " << _cellsPerDimension[1] << " / " << _cellsPerDimension[2] << std::endl;
      global_log->error() << "_haloWidthInNumCells: " << _haloWidthInNumCells[0] << " / " << _haloWidthInNumCells[1] << " / " << _haloWidthInNumCells[2]
                          << std::endl;
      exit(5);
    }
    this->_localInsertionsMinusDeletions = 0;

    initializeCells();
    calculateNeighbourIndices();
    _cellsValid = false;
  }

  ~LinkedCellsForCoupling() {}

  // returns a pointer to the ParticleCell at the provided index
  ParticleCell *getCellPointer(int index) { return &_cells[index]; }
};

#endif /* COUPLING_INTERFACE_IMPL_MARDYN_LINKEDCELLSFORCOUPLING_H_ */
