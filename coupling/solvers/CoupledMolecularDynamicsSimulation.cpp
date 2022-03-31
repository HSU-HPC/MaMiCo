// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#include "coupling/solvers/CoupledMolecularDynamicsSimulation.h"

coupling::solvers::CoupledMolecularDynamicsSimulation::
    CoupledMolecularDynamicsSimulation(
        const simplemd::configurations::MolecularDynamicsConfiguration &
            configuration)
    : simplemd::MolecularDynamicsSimulation(configuration),
      _macroscopicCellService(NULL), _couplingSwitchedOn(true) {}

void coupling::solvers::CoupledMolecularDynamicsSimulation::
    simulateOneCouplingTimestep(const unsigned int &t) {
  // if coupling is switched off, perform "normal" MD timestep
  if (!_couplingSwitchedOn) {
    simulateOneTimestep(t);
    return;
  }
  if (_parallelTopologyService->isIdle()) {
    return;
  }

  _boundaryTreatment->putBoundaryParticlesToInnerCellsAndFillBoundaryCells(
      _localBoundary, *_parallelTopologyService);

  // call to synchronise data in cells; needs to be at this point of the
  // coupling algorithm as the particles need to be placed inside
  // the correct sampling volumes (hence: after communication with neighbours
  // and molecule updates);
  // do it BEFORE quantities are manipulated as we can then also do some
  // pre-processing here.
  _macroscopicCellService->processInnerMacroscopicCellAfterMDTimestep();

  // ------------ coupling step: distribute mass ---------------------
  _macroscopicCellService->distributeMass(t);

  // for isothermal simulations: apply thermostat
  _macroscopicCellService->applyTemperatureToMolecules(t);

  // ---------- from here: go on with usual MD algorithm
  // ------------------------------

  // compute forces. After this step, each molecule has received all force
  // contributions from its neighbors.
  _linkedCellService->iterateCellPairs(*_lennardJonesForce, false);

  // distribute momentum -> some methods require modification of force terms,
  // therefore we call it AFTER the force computation and before everything else
  _macroscopicCellService->distributeMomentum(t);

  // apply boundary forces
  _macroscopicCellService->applyBoundaryForce(t);

  // evaluate statistics
  evaluateStatistics(t);

  _boundaryTreatment->emptyGhostBoundaryCells();

  // plot VTK output
  if ((_configuration.getVTKConfiguration().getWriteEveryTimestep() != 0) &&
      (t % _configuration.getVTKConfiguration().getWriteEveryTimestep() == 0)) {
    _vtkMoleculeWriter->setTimestep(t);
    _moleculeService->iterateMolecules(*_vtkMoleculeWriter, false);
  }
  // write checkpoint
  if ((_configuration.getCheckpointConfiguration().getWriteEveryTimestep() !=
       0) && (t % _configuration.getCheckpointConfiguration()
                      .getWriteEveryTimestep() == 0)) {
    _moleculeService->writeCheckPoint(
        *_parallelTopologyService,
        _configuration.getCheckpointConfiguration().getFilename(), t);
  }

  // reorganise memory if needed
  if ((_configuration.getSimulationConfiguration()
           .getReorganiseMemoryEveryTimestep() != 0) &&
      (t % _configuration.getSimulationConfiguration()
               .getReorganiseMemoryEveryTimestep() == 0)) {
    _moleculeService->reorganiseMemory(*_parallelTopologyService,
                                       *_linkedCellService);
  }

  // plot also macroscopic cell information
  _macroscopicCellService->plotEveryMicroscopicTimestep(t);

  _linkedCellService->iterateCells(*_emptyLinkedListsMapping, false);

  // time integration. After this step, the velocities and the positions of the
  // molecules have been updated.
  _moleculeService->iterateMolecules(*_timeIntegrator, false);

  // sort molecules into linked cells
  _moleculeService->iterateMolecules(*_updateLinkedCellListsMapping, false);

  if (_parallelTopologyService->getProcessCoordinates() ==
      tarch::la::Vector<MD_DIM, unsigned int>(0)) {
    //if(t%50==0) std::cout <<"Finish MD timestep " << t << "..." << std::endl;
  }
}
