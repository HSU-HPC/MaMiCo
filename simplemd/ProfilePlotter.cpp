// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#include "simplemd/ProfilePlotter.h"

simplemd::ProfilePlotter::ProfilePlotter(const std::vector<simplemd::configurations::ProfilePlotterConfiguration>& configurations,
                                         const simplemd::services::ParallelTopologyService& parallelTopologyService,
                                         simplemd::services::LinkedCellService& linkedCellService, const double& linkedCellVolume,
                                         const unsigned int& localMDSimulation)
    : _linkedCellService(linkedCellService) {
  for (unsigned int i = 0; i < _plotters.size(); i++) {
    if (_plotters[i] != NULL) {
      delete _plotters[i];
      _plotters[i] = NULL;
    }
  }
  _plotters.clear();

  for (unsigned int i = 0; i < configurations.size(); i++) {
    _plotters.push_back(new simplemd::cellmappings::ProfilePlotterMapping(parallelTopologyService, linkedCellService, configurations[i].getWriteEveryTimestep(),
                                                                          configurations[i].getSampleEveryTimestep(), configurations[i].getStartAtTimestep(),
                                                                          linkedCellVolume, localMDSimulation));
    _ranges.push_back(configurations[i].getRange());
    _startCells.push_back(configurations[i].getStartCell());
  }

#if (MD_PARALLEL == MD_YES)
  adjustProfilesInParallel(parallelTopologyService);
#endif
}

simplemd::ProfilePlotter::~ProfilePlotter() {
  for (unsigned int i = 0; i < _plotters.size(); i++) {
    if (_plotters[i] != NULL) {
      delete _plotters[i];
      _plotters[i] = NULL;
    }
  }
  _plotters.clear();
}

#if (MD_PARALLEL == MD_YES)
void simplemd::ProfilePlotter::adjustProfilesInParallel(const simplemd::services::ParallelTopologyService& parallelTopologyService) {
  const unsigned int size = (unsigned int)_plotters.size();
  tarch::la::Vector<MD_DIM, unsigned int> localStartCell;
  tarch::la::Vector<MD_DIM, unsigned int> localRange;
  for (unsigned int i = 0; i < size; i++) {
    _isActive.push_back(parallelTopologyService.globalToLocalRegionOfInterest(_startCells[i], _ranges[i], localStartCell, localRange));
    _startCells[i] = localStartCell; // will be 0 if not active
    _ranges[i] = localRange;         // will be 0 if not active
  }
}
#endif

void simplemd::ProfilePlotter::accumulateAndPlotInformation(const unsigned int& t) {
  const unsigned int size = (unsigned int)_plotters.size();
  for (unsigned int i = 0; i < size; i++) {
#if (MD_PARALLEL == MD_YES)
    if (!_isActive[i])
      continue; // skip profiles which have no intersection with local domain
#endif
    _plotters[i]->setCurrentTimestep(t);
    _linkedCellService.iterateCells<simplemd::cellmappings::ProfilePlotterMapping>(*_plotters[i], _startCells[i], _ranges[i], false);
  }
}
