// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#include "simplemd/BoundaryTreatment.h"

void simplemd::BoundaryTreatment::putBoundaryParticlesToInnerCells(const tarch::la::Vector<MD_LINKED_CELL_NEIGHBOURS, simplemd::BoundaryType>& boundary,
                                                                   simplemd::services::ParallelTopologyService& parallelTopologyService) {
  _periodicBoundaryMapping.setDomainSize(parallelTopologyService.getGlobalDomainSize());
  _periodicBoundaryMapping.setProcessCoordinates(parallelTopologyService.getProcessCoordinates());
  _periodicBoundaryMapping.setNumberOfProcesses(parallelTopologyService.getNumberOfProcesses());
  applyMappingToBoundaryCells(boundary, simplemd::PERIODIC_BOUNDARY, false, _periodicBoundaryMapping);
  // collect molecules in open-boundary region and remove them from the simulation.
  _collectMoleculesMapping.reset();
  applyMappingToBoundaryCells(boundary, simplemd::OPEN_BOUNDARY, false, _collectMoleculesMapping);

#if (MD_PARALLEL == MD_YES)
  applyMappingToBoundaryCells(boundary, simplemd::PARALLEL_BOUNDARY, false, _parallelBoundaryMapping);

  // after data from periodic/ parallel boundaries have been sent, receive them and put them into the local
  // data management systems
  parallelTopologyService.communicationSteps_1_2();
  parallelTopologyService.communicationSteps_3_4(_moleculeService, _linkedCellService);
#endif
}

void simplemd::BoundaryTreatment::fillBoundaryCells(const tarch::la::Vector<MD_LINKED_CELL_NEIGHBOURS, simplemd::BoundaryType>& boundary,
                                                    simplemd::services::ParallelTopologyService& parallelTopologyService) {
  tarch::la::Vector<MD_DIM, unsigned int> startOuter;
  tarch::la::Vector<MD_DIM, unsigned int> numberCellsOuter;
  _fillCellsMapping.setDomainSize(parallelTopologyService.getGlobalDomainSize());

  applyMappingToOutermostNonBoundaryCells(false, _fillCellsMapping);

#if (MD_PARALLEL == MD_YES)
  parallelTopologyService.communicationSteps_1_2();
  parallelTopologyService.communicationSteps_3_4(_moleculeService, _linkedCellService);
#endif
}

void simplemd::BoundaryTreatment::emptyGhostBoundaryCells() {
  tarch::la::Vector<MD_DIM, unsigned int> startOuter;
  tarch::la::Vector<MD_DIM, unsigned int> numberCellsOuter;
#if (MD_DIM == 1)
  startOuter[0] = 0;
  numberCellsOuter[0] = 1;
  _linkedCellService.iterateCells(_deleteMoleculesMapping, startOuter, numberCellsOuter, false);
  startOuter[0] = _linkedCellService.getLocalNumberOfCells()[0] + 2 * _linkedCellService.getLocalIndexOfFirstCell()[0] - 1;
  _linkedCellService.iterateCells(_deleteMoleculesMapping, startOuter, numberCellsOuter, false);
#endif

#if (MD_DIM == 2)
  // lower edge
  startOuter[0] = 0;
  startOuter[1] = 0;
  numberCellsOuter[0] = _linkedCellService.getLocalNumberOfCells()[0] + 2 * _linkedCellService.getLocalIndexOfFirstCell()[0];
  numberCellsOuter[1] = 1;
  _linkedCellService.iterateCells(_deleteMoleculesMapping, startOuter, numberCellsOuter, false);
  // upper edge
  startOuter[0] = 0;
  startOuter[1] = _linkedCellService.getLocalNumberOfCells()[1] + 2 * _linkedCellService.getLocalIndexOfFirstCell()[1] - 1;
  _linkedCellService.iterateCells(_deleteMoleculesMapping, startOuter, numberCellsOuter, false);
  // left edge
  startOuter[0] = 0;
  startOuter[1] = 1;
  numberCellsOuter[0] = 1;
  numberCellsOuter[1] = _linkedCellService.getLocalNumberOfCells()[1];
  _linkedCellService.iterateCells(_deleteMoleculesMapping, startOuter, numberCellsOuter, false);
  // right edge
  startOuter[0] = _linkedCellService.getLocalNumberOfCells()[0] + 2 * _linkedCellService.getLocalIndexOfFirstCell()[0] - 1;
  startOuter[1] = 1;
  _linkedCellService.iterateCells(_deleteMoleculesMapping, startOuter, numberCellsOuter, false);
#endif

#if (MD_DIM == 3)
  // lower face
  numberCellsOuter[0] = _linkedCellService.getLocalNumberOfCells()[0] + 2 * _linkedCellService.getLocalIndexOfFirstCell()[0];
  numberCellsOuter[1] = _linkedCellService.getLocalNumberOfCells()[1] + 2 * _linkedCellService.getLocalIndexOfFirstCell()[1];
  numberCellsOuter[2] = 1;
  startOuter = tarch::la::Vector<MD_DIM, unsigned int>(0);
  _linkedCellService.iterateCells(_deleteMoleculesMapping, startOuter, numberCellsOuter, false);
  // top face
  startOuter[0] = 0;
  startOuter[1] = 0;
  startOuter[2] = _linkedCellService.getLocalNumberOfCells()[2] + 2 * _linkedCellService.getLocalIndexOfFirstCell()[2] - 1;
  _linkedCellService.iterateCells(_deleteMoleculesMapping, startOuter, numberCellsOuter, false);
  // left face
  numberCellsOuter[0] = 1;
  numberCellsOuter[1] = _linkedCellService.getLocalNumberOfCells()[1] + 2 * _linkedCellService.getLocalIndexOfFirstCell()[1];
  numberCellsOuter[2] = _linkedCellService.getLocalNumberOfCells()[2];
  startOuter[0] = 0;
  startOuter[1] = 0;
  startOuter[2] = 1;
  _linkedCellService.iterateCells(_deleteMoleculesMapping, startOuter, numberCellsOuter, false);
  // right face
  startOuter[0] = _linkedCellService.getLocalNumberOfCells()[0] + 2 * _linkedCellService.getLocalIndexOfFirstCell()[0] - 1;
  startOuter[1] = 0;
  startOuter[2] = 1;
  _linkedCellService.iterateCells(_deleteMoleculesMapping, startOuter, numberCellsOuter, false);
  // front face
  numberCellsOuter[0] = _linkedCellService.getLocalNumberOfCells()[0];
  numberCellsOuter[1] = 1;
  numberCellsOuter[2] = _linkedCellService.getLocalNumberOfCells()[2];
  startOuter[0] = 1;
  startOuter[1] = 0;
  startOuter[2] = 1;
  _linkedCellService.iterateCells(_deleteMoleculesMapping, startOuter, numberCellsOuter, false);
  // back face
  startOuter[0] = 1;
  startOuter[1] = _linkedCellService.getLocalNumberOfCells()[1] + 2 * _linkedCellService.getLocalIndexOfFirstCell()[1] - 1;
  startOuter[2] = 1;
  _linkedCellService.iterateCells(_deleteMoleculesMapping, startOuter, numberCellsOuter, false);
#endif
}

std::list<simplemd::Molecule> simplemd::BoundaryTreatment::getEscapedMolecules() const { return _collectMoleculesMapping.getCollectedMolecules(); }

void simplemd::BoundaryTreatment::putBoundaryParticlesToInnerCellsAndFillBoundaryCells(
    const tarch::la::Vector<MD_LINKED_CELL_NEIGHBOURS, simplemd::BoundaryType>& boundary,
    simplemd::services::ParallelTopologyService& parallelTopologyService) {

  // putBoundaryParticlesToInnerCells
  _fillCellsMapping.setDomainSize(parallelTopologyService.getGlobalDomainSize());

  // iterating with _fillCellsMapping, not _periodicBoundaryMapping or _parallelBoundaryMapping!
  applyMappingToBoundaryCells(boundary, simplemd::PERIODIC_BOUNDARY, false, _fillCellsMapping);
  // collect molecules in open-boundary region and remove them from the simulation.
  _collectMoleculesMapping.reset();
  applyMappingToBoundaryCells(boundary, simplemd::OPEN_BOUNDARY, false, _collectMoleculesMapping);

#if (MD_PARALLEL == MD_YES)
  applyMappingToBoundaryCells(boundary, simplemd::PARALLEL_BOUNDARY, false, _fillCellsMapping);
#endif

  // fillBoundaryCells
  applyMappingToOutermostNonBoundaryCells(false, _fillCellsMapping);

  // handle local periodic ghost cells
  parallelTopologyService.unpackLocalBuffer(_moleculeService, _linkedCellService);

#if (MD_PARALLEL == MD_YES)
  parallelTopologyService.communicationSteps_1_2();
  parallelTopologyService.communicationSteps_3_4(_moleculeService, _linkedCellService);
#endif
}

void simplemd::BoundaryTreatment::putBoundaryParticlesToInnerCellsFillBoundaryCellsAndOverlapWithForceComputations(
    const tarch::la::Vector<MD_LINKED_CELL_NEIGHBOURS, simplemd::BoundaryType>& boundary, simplemd::services::ParallelTopologyService& parallelTopologyService,
#if (MD_BODY == 2)
    simplemd::cellmappings::LennardJonesForceMapping& lennardJonesForce,
#else
    simplemd::cellmappings::AxilrodTellerForceMapping& axilrodTellerForce,
#endif
    const bool& useOpenMP) {
  // prevent execution in case domain is too small
  const tarch::la::Vector<MD_DIM, unsigned int> localNumberOfCells(_linkedCellService.getLocalNumberOfCells());
  for (unsigned int d = 0; d < MD_DIM; d++) {
    if (localNumberOfCells[d] < 3) {
      std::cout << "Fatal Error: BoundaryTreatment::putBoundaryParticlesToInnerCellsFillBoundaryCellsAndOverlapWithForceComputations: " << std::endl;
      std::cout << "in order for overlapping process communications with computations algorithm to work, there have to be at least 3 local linked cells along "
                   "every dimension on every processor"
                << std::endl;
      std::cout << "Current local number of linked cells: " << localNumberOfCells << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  // putBoundaryParticlesToInnerCells
  _fillCellsMapping.setDomainSize(parallelTopologyService.getGlobalDomainSize());

  applyMappingToBoundaryCells(boundary, simplemd::PERIODIC_BOUNDARY, false, _fillCellsMapping);
  // collect molecules in open-boundary region and remove them from the simulation.
  _collectMoleculesMapping.reset();
  applyMappingToBoundaryCells(boundary, simplemd::OPEN_BOUNDARY, false, _collectMoleculesMapping);

#if (MD_PARALLEL == MD_YES)
  applyMappingToBoundaryCells(boundary, simplemd::PARALLEL_BOUNDARY, false, _fillCellsMapping);
#endif

  // fillBoundaryCells
  applyMappingToOutermostNonBoundaryCells(false, _fillCellsMapping);

// send out buffers, but don't wait for them
#if (MD_PARALLEL == MD_YES)
#if (MD_DEBUG == MD_YES)
  std::cout << "BoundaryTreatment::putBoundaryParticlesToInnerCellsFillBoundaryCellsAndOverlapWithForceComputations: issuing send and receive calls on "
               "communication buffers."
            << std::endl;
#endif
  parallelTopologyService.communicationSteps_1_2();
#endif

// handle local periodic ghost cells, if any
#if (MD_DEBUG == MD_YES)
  std::cout << "BoundaryTreatment::putBoundaryParticlesToInnerCellsFillBoundaryCellsAndOverlapWithForceComputations: unpacking local buffer" << std::endl;
#endif
  parallelTopologyService.unpackLocalBuffer(_moleculeService, _linkedCellService);

  // compute forces in inner cells, where boundary and process-leaving particles are not needed.
#if (MD_BODY == 2)
  applyMappingToCommunicationIndependentCells(useOpenMP, lennardJonesForce);
#else
  applyMappingToCommunicationIndependentCells(useOpenMP, axilrodTellerForce);
#endif

// wait for buffers
#if (MD_PARALLEL == MD_YES)
#if (MD_DEBUG == MD_YES)
  std::cout << "BoundaryTreatment::putBoundaryParticlesToInnerCellsFillBoundaryCellsAndOverlapWithForceComputations: begin waiting for buffer requests to be "
               "fulfilled."
            << std::endl;
#endif
  parallelTopologyService.communicationSteps_3_4(_moleculeService, _linkedCellService);
#endif

  // compute rest of forces
#if (MD_BODY == 2)
  applyMappingToCommunicationDependentCells(useOpenMP, lennardJonesForce);
#else
  applyMappingToCommunicationDependentCells(useOpenMP, axilrodTellerForce);
#endif
}
