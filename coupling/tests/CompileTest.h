// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_TESTS_COMPILETEST_H_
#define _MOLECULARDYNAMICS_COUPLING_TESTS_COMPILETEST_H_

#include "Test.h"
#include "coupling/configurations/MacroscopicCellConfiguration.h"
#include "coupling/services/MacroscopicCellService.h"
#include "simplemd/LinkedCell.h"

class CompileTest : public Test {
public:
  CompileTest() : Test("CompileTest") {}
  virtual ~CompileTest() {}

  virtual void run() {
    const bool youAreStupid = false;

    // never execute the functions, but include them in build
    if (youAreStupid) {
      compileTest<2>();
      compileTest<3>();
    }
  }

  template <unsigned int dim> void compileTest() {
    coupling::interface::MDSolverInterface<simplemd::LinkedCell, dim>* mdSolverInterface = NULL;
    coupling::interface::MacroscopicSolverInterface<dim>* macroscopicSolverInterface = NULL;
    tarch::la::Vector<dim, unsigned int> numberProcesses(1);
    unsigned int rank = 0;
    const coupling::configurations::ParticleInsertionConfiguration particleInsertionConfiguration;
    const coupling::configurations::MomentumInsertionConfiguration momentumInsertionConfiguration;
    const coupling::configurations::BoundaryForceConfiguration<dim> boundaryForceConfiguration;
    const coupling::configurations::TransferStrategyConfiguration<dim> transferStrategyConfiguration;
    const coupling::configurations::ParallelTopologyConfiguration parallelTopologyConfiguration;
    unsigned int numberMDTimestepsPerCouplingCycle = 100;
    const coupling::configurations::MacroscopicCellConfiguration<dim> macroscopicCellConfiguration;
    std::vector<coupling::datastructures::MacroscopicCell<dim>*> cells;
    unsigned int* globalCellIndices = NULL;

    coupling::services::MacroscopicCellService<dim>* macroscopicCellService = new coupling::services::MacroscopicCellServiceImpl<simplemd::LinkedCell, dim>(
        0, mdSolverInterface, macroscopicSolverInterface, numberProcesses, rank, particleInsertionConfiguration, momentumInsertionConfiguration,
        boundaryForceConfiguration, transferStrategyConfiguration, parallelTopologyConfiguration, numberMDTimestepsPerCouplingCycle,
        macroscopicCellConfiguration);

    macroscopicCellService->sendFromMacro2MD(cells, globalCellIndices);
    macroscopicCellService->sendFromMD2Macro(cells, globalCellIndices);
    macroscopicCellService->processInnerMacroscopicCellAfterMDTimestep();
    macroscopicCellService->computeAndStoreTemperature(0.0);
    macroscopicCellService->applyTemperatureToMolecules(0);
    macroscopicCellService->distributeMass(0);
    macroscopicCellService->distributeMomentum(0);
    macroscopicCellService->applyBoundaryForce(0);
    delete macroscopicCellService;
  }
};
#endif // _MOLECULARDYNAMICS_COUPLING_TESTS_COMPILETEST_H_
