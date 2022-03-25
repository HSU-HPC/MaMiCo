// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef TESTMARDYN_H_
#define TESTMARDYN_H_

#include "Simulation.h"

#include "particleContainer/LinkedCells.h"
#include "utils/FileUtils.h"

#include "coupling/CouplingMDDefinitions.h"
#include "coupling/configurations/MaMiCoConfiguration.h"
#include "coupling/interface/MamicoInterfaceProvider.h"
#include "coupling/interface/impl/macroscopictestsolvers/VoidMacroscopicSolverInterface.h"
#include "coupling/services/MacroscopicCellService.h"
#include "coupling/tests/Test.h"

#include "coupling/interface/impl/MarDyn/MarDynCell.h"
#include "coupling/interface/impl/MarDyn/MarDynMDSolverInterface.h"

using Log::global_log;

/**
 * 	test class for the MarDyn MD coupling interface
 * 	initializes all parts of the coupling and tests the get methods of the
 *MDSolverInterface
 *	@author Hanno Flohr
 */
class TestMarDyn : public Test {

public:
  TestMarDyn(int argc, char **argv, std::string name)
      : Test(name),
        //_mamicoInterface(NULL),
        _marDyn(NULL), _marDynMDsolver(NULL), _macroSolverInterface(NULL), _macroscopicCellService(NULL), _argc(argc), _argv(argv) {}

  virtual ~TestMarDyn() {
    if (_macroscopicCellService != NULL) {
      delete _macroscopicCellService;
      _macroscopicCellService = NULL;
    }
    if (_macroSolverInterface != NULL) {
      delete _macroSolverInterface;
      _macroSolverInterface = NULL;
    }

    if (_marDynMDsolver != NULL) {
      delete _marDynMDsolver;
      _marDynMDsolver = NULL;
    }
    if (_marDyn != NULL) {
      delete _marDyn;
      _marDyn = NULL;
    }
  };

  virtual void run() {
    std::cout << "Load MarDyn test configuration.. " << std::endl;
    loadMarDynTestConfiguration("mardyn_dummy_coupling.cfg", 10);

    std::cout << "Testing MD solver get methods: " << std::endl;
    MarDynMDSolverInterface *mdsi =
        (MarDynMDSolverInterface *)coupling::interface::MamicoInterfaceProvider<MarDynCell, 3>::getInstance().getMDSolverInterface();
    if (mdsi == NULL) {
      std::cout << "ERROR: MD solver is null! Aborting." << std::endl;
      exit(EXIT_FAILURE);
    }

    // test variables
    bool validResult;
    double testValue;
    tarch::la::Vector<3, double> testVec(0.0);

    // test getMoleculeMass()
    testValue = mdsi->getMoleculeMass();
    std::cout << "Molecule mass: " << testValue << std::endl;
    if (testValue != 1.0)
      std::cout << "Wrong mass, should be 1.0!" << std::endl;

    // test getKB()
    testValue = mdsi->getKB();
    std::cout << "kB: " << testValue << std::endl;
    if (testValue != 1.0)
      std::cout << "Wrong Boltzmann's constant, should be 1.0!" << std::endl;

    // test getMoleculeSigma()
    testValue = mdsi->getMoleculeSigma();
    std::cout << "Molecule sigma: " << testValue << std::endl;
    if (testValue != 1.0)
      std::cout << "Wrong molecule sigma value, should be 1.0!" << std::endl;

    // test getMoleculeEpsilon()
    testValue = mdsi->getMoleculeEpsilon();
    std::cout << "Molecule epsilon: " << testValue << std::endl;
    if (testValue != 1.0)
      std::cout << "Wrong molecule epsilon value, should be 1.0!" << std::endl;

    // test getDt
    testValue = mdsi->getDt();
    std::cout << "Timestep: " << testValue << std::endl;
    if (testValue != 0.002)
      std::cout << "Wrong timestep size, should be 0.002!" << std::endl;

    // test getGlobalMDDomainSize()
    testVec = mdsi->getGlobalMDDomainSize();
    std::cout << "GlobalMDDomainSize: " << testVec << std::endl;
    validResult = (testVec[0] == 30.0);
    for (int d = 0; d < 3; d++) {
      validResult = validResult && (testVec[d] == 30.0);
    }
    if (!validResult)
      std::cout << "Wrong global MD domain size! Should be 30.0 in each component!" << std::endl;

    // test getGlobalMDDomainOffset()
    testVec = mdsi->getGlobalMDDomainOffset();
    std::cout << "GlobalMDDomainOffset: " << testVec << std::endl;
    validResult = (testVec[0] == 0.0);
    for (int d = 0; d < 3; d++) {
      validResult = validResult && (testVec[d] == 0.0);
    }
    if (!validResult)
      std::cout << "Wrong global MD domain offset! Should be 0.0 in each component!" << std::endl;

    std::cout << "Load Macroscopic solver test configuration.." << std::endl;
    loadMacroscopicSolverConfiguration();

    std::cout << "Load MaMiCo test configuration and init macroscopic cell service.." << std::endl;
    loadMamicoTestConfiguration();

    std::cout << "Initialization done." << std::endl;
  }

  void loadTestConfiguration() {
    std::cout << "Initializing test configuration: " << std::endl;
    std::cout << "Load MarDyn test configuration.." << std::endl;
    loadMarDynTestConfiguration("test100.cfg", 1);

    std::cout << "Load Macroscopic solver test configuration.." << std::endl;
    loadMacroscopicSolverConfiguration();

    std::cout << "Load MaMiCo test configuration and init macroscopic cell service.." << std::endl;
    loadMamicoTestConfiguration();

    std::cout << "Initialization done." << std::endl;
  }

protected:
  void loadMarDynTestConfiguration(std::string testConfiguration, unsigned long int timesteps) {
    global_log = new Log::Logger(Log::Info);
    if (_marDyn != NULL) {
      delete _marDyn;
      _marDyn = NULL;
    }
    tarch::la::Vector<3, double> mamicoCellSize(2.5);
    tarch::la::Vector<3, unsigned int> linkedCellsPerMacroscopicCell(1);
    _marDyn = new MarDynCoupledSimulation(mamicoCellSize, linkedCellsPerMacroscopicCell);

    // initialize mardyn
    if (fileExists(testConfiguration.c_str())) {
      _marDyn->initConfigOldstyle(testConfiguration);
    } else {
      std::cout << "ERROR: Cannot open input file: '" << testConfiguration << "'" << std::endl;
      exit(EXIT_FAILURE);
    }
    _marDyn->setNumTimesteps(timesteps);
    _marDyn->prepare_start();

    // initialize MardynMDSolver
    if (_marDynMDsolver != NULL) {
      delete _marDynMDsolver;
      _marDynMDsolver = NULL;
    }
    _marDynMDsolver = new MarDynMDSolverInterface(_marDyn);

    // set MarDynMDSolver as MDSolverInterface
    coupling::interface::MamicoInterfaceProvider<MarDynCell, 3>::getInstance().setMDSolverInterface(_marDynMDsolver);
    std::cout << "MarDyn MD solver initialized and set in interface." << std::endl;
  }

  /** loads a dummy solver interface which basically does nothing */
  void loadMacroscopicSolverConfiguration() {
    if (_macroSolverInterface != NULL) {
      delete _macroSolverInterface;
      _macroSolverInterface = NULL;
    }
    _macroSolverInterface = new coupling::interface::VoidMacroscopicSolverInterface<3>();
    if (_macroSolverInterface == NULL) {
      std::cout << "ERROR in TestMarDyn: macroSolverInterface==NULL!" << std::endl;
      exit(EXIT_FAILURE);
    }
    coupling::interface::MamicoInterfaceProvider<MarDynCell, 3>::getInstance().setMacroscopicSolverInterface(_macroSolverInterface);
  }

  /** loads the MaMiCo configuration and initializes the macroscopic cells
   *service; should be called after loadMarDynTestConfiguration() and
   *loadMacroscopicSolverConfiguration().
   */
  void loadMamicoTestConfiguration() {
    const std::string mamicoTestConfig = "mamico_mardyn_test_configuration.xml";
    tarch::la::Vector<3, unsigned int> numberProcesses(1);
    int rank = 0;

    const unsigned int numberMDTimestepsPerCouplingCycle = 10;
    coupling::configurations::MaMiCoConfiguration<3> config;
    std::cout << "Parse config: " << mamicoTestConfig << std::endl;
    tarch::configuration::ParseConfiguration::parseConfiguration<coupling::configurations::MaMiCoConfiguration<3>>(mamicoTestConfig, "mamico", config);

    std::cout << "Init macroscopic cell service.." << std::endl;
    if (_macroscopicCellService != NULL) {
      delete _macroscopicCellService;
      _macroscopicCellService = NULL;
    }
    _macroscopicCellService = new coupling::services::MacroscopicCellServiceImpl<MarDynCell, 3>(
        0, coupling::interface::MamicoInterfaceProvider<MarDynCell, 3>::getInstance().getMDSolverInterface(),
        coupling::interface::MamicoInterfaceProvider<MarDynCell, 3>::getInstance().getMacroscopicSolverInterface(), numberProcesses, rank,
        config.getParticleInsertionConfiguration(), config.getMomentumInsertionConfiguration(), config.getBoundaryForceConfiguration(),
        config.getTransferStrategyConfiguration(), config.getParallelTopologyConfiguration(), numberMDTimestepsPerCouplingCycle,
        config.getMacroscopicCellConfiguration());
    if (_macroscopicCellService == NULL) {
      std::cout << "ERROR TestMarDyn: _macroscopicCellService==NULL!" << std::endl;
      exit(EXIT_FAILURE);
    }
    std::cout << "Set macroscopic cell service in MaMiCoInterfaceProvider.." << std::endl;
    coupling::interface::MamicoInterfaceProvider<MarDynCell, 3>::getInstance().setMacroscopicCellService(_macroscopicCellService);

    coupling::services::MacroscopicCellService<3> *macroCellService =
        coupling::interface::MamicoInterfaceProvider<MarDynCell, 3>::getInstance().getMacroscopicCellService();
    std::cout << "init test macro cell service: " << (macroCellService == NULL) << std::endl;
  }

  // MarDyn simulation
  MarDynCoupledSimulation *_marDyn;
  // MarDyn MD solver interface
  MarDynMDSolverInterface *_marDynMDsolver;
  // test solver interface
  coupling::interface::VoidMacroscopicSolverInterface<3> *_macroSolverInterface;
  coupling::services::MacroscopicCellService<3> *_macroscopicCellService;
  // command line arguments
  int _argc;
  char **_argv;
};

#endif /* TESTMARDYN_H_ */
