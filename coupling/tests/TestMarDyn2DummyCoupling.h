// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef TESTMARDYN2DUMMYCOUPLING_H_
#define TESTMARDYN2DUMMYCOUPLING_H_

#include <string>
#include <iostream>

#include "tarch/la/Vector.h"
//MarDyn
#include "utils/FileUtils.h"
#include "utils/Timer.h"
#include "particleContainer/LinkedCells.h"
//MaMiCo
#include "coupling/tests/Test.h"
#include "tarch/configuration/ParseConfiguration.h"
#include "coupling/solvers/DummySolver.h"
#include "coupling/solvers/DummySolverInterfaceService.h"
#include "coupling/configurations/MaMiCoConfiguration.h"
#include "coupling/services/MacroscopicCellService.h"
#include "coupling/CouplingMDDefinitions.h"
#include "coupling/datastructures/MacroscopicCell.h"
#include "coupling/datastructures/MacroscopicCells.h"
#include "coupling/interface/impl/MarDyn/MarDynMDSolverInterface.h"
#include "coupling/interface/MamicoInterfaceProvider.h"

//#include "MDSimulationFactory.h"

#ifdef ENABLE_MPI
#include <mpi.h>
#include "parallel/DomainDecomposition.h"
#endif

using Log::global_log;

/*	Test class. Provides an implementation of the coupling of Mardyn with input
 * values created by a Dummy solver.
 * 	@author Hanno Flohr
 */
class TestMarDyn2DummyCoupling : public Test {
public:
  TestMarDyn2DummyCoupling(std::string name, int argc, char **argv)
      : Test(name), _marDynMDSolverInterface(NULL), _dummySolver(NULL),
        _macroscopicCellService(NULL), _argc(argc), _argv(argv) {}

  virtual ~TestMarDyn2DummyCoupling() {}

  virtual void run() {
    std::cout << "Testing MarDyn coupled with Dummy Macroscopic solver: "
              << std::endl;

    //initialize solvers and macroscopic cell service
    init();

    Timer mamicoTimer;

    //test coupled simulation
    mamicoTimer.start();
    simulate();
    mamicoTimer.stop();

    global_log->info() << "The coupling step took: " << mamicoTimer.get_etime()
                       << " sec." << std::endl;

    shutdown();
  }

private:

  //executes one coupling step
  void simulate() {
    int rank = 0;
#ifdef ENABLE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    //add information to buffers for send operation from DummySolver to MD
    if (rank == 0)
      storeDummySolverDataInSendBuffer();

    //send DummySolver information to MD side
    _macroscopicCellService->sendFromMacro2MD(
        DummySolverInterfaceService::getInstance().getSendBuffer(),
        DummySolverInterfaceService::getInstance()
            .getGlobalCellIndices4SendBuffer());

    _marDynMDSolverInterface->getSimulation()->setInitValuesForCoupling();
    //execute md timesteps
    unsigned int timesteps = 250000;
    for (unsigned int t = 0; t < timesteps; t++)
      _marDynMDSolverInterface->getSimulation()->simulateOneMDTimestep(t);

    //send MD information to DummySolver
    _macroscopicCellService->sendFromMD2Macro(
        DummySolverInterfaceService::getInstance().getReceiveBuffer(),
        DummySolverInterfaceService::getInstance()
            .getGlobalCellIndices4ReceiveBuffer());

#ifdef ENABLE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    //add information to buffers for send operation from MD to DummySolver
    if (rank == 0)
      writeReceiveBufferDataToDummySolver();

#ifdef ENABLE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    if (rank == 0)
      writeDummySolverInfo(timesteps);

    //does NOT work with activated vectorization!
    //_macroscopicCellService->plotEveryMacroscopicTimestep(timesteps);
  }

  void init() {
    //parse mamico config file
    const std::string mamicoConfigFile = "mamico_mardyn_test_configuration.xml";
    coupling::configurations::MaMiCoConfiguration<3> mamicoConfig;
    tarch::configuration::ParseConfiguration::parseConfiguration<
        coupling::configurations::MaMiCoConfiguration<3> >(
        mamicoConfigFile, "mamico", mamicoConfig);

    //initialize mardyn simulation
    const std::string marDynConfigFile = "mardyn_dummy_coupling.cfg";
    global_log = new Log::Logger(Log::Info);

#ifdef ENABLE_MPI
    global_log->set_mpi_output_root(0);
#endif

    std::cout.precision(8);
    MarDynCoupledSimulation *marDynSimulation = new MarDynCoupledSimulation(
        mamicoConfig.getMacroscopicCellConfiguration().getMacroscopicCellSize(),
        mamicoConfig.getMacroscopicCellConfiguration()
            .getNumberLinkedCellsPerMacroscopicCell());

    if (fileExists(marDynConfigFile.c_str())) {
      marDynSimulation->initConfigOldstyle(marDynConfigFile);
    } else {
      std::cout << "ERROR: Cannot open MarDyn input file: '" << marDynConfigFile
                << "'" << std::endl;
      exit(EXIT_FAILURE);
    }
    marDynSimulation->setNumTimesteps(250000);
    marDynSimulation->prepare_start();

    //initialize mardyn md solver interface
    if (_marDynMDSolverInterface != NULL) {
      delete _marDynMDSolverInterface;
      _marDynMDSolverInterface = NULL;
    }
    _marDynMDSolverInterface = new MarDynMDSolverInterface(marDynSimulation);

    coupling::interface::MamicoInterfaceProvider<MarDynCell, 3>::getInstance()
        .setMDSolverInterface(_marDynMDSolverInterface);

#ifdef ENABLE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    //initialize dummy macroscopic solver
    if (_dummySolver != NULL) {
      delete _dummySolver;
      _dummySolver = NULL;
    }
    _dummySolver = new DummySolver(18, 18, 18, 0.48);

    //initialize macroscopic solver interface
    tarch::la::Vector<3, unsigned int> numberProcesses(1);
    int rank = 0;

    //#ifdef ENABLE_MPI
    //		int size;
    //		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //		MPI_Comm_size(MPI_COMM_WORLD, &size);
    //
    //		//print the decomposition to file "decomp_test"
    //		marDynSimulation->domainDecomposition().printDecomp("decomp_test",marDynSimulation->getDomain());
    //
    //		//set the number of processes using the same method as done in MarDyn
    //(MPI_Dims_create)
    //		int gridSize [3];
    //		for(int d=0;d<3;d++) gridSize[d] = 0;
    //		MPI_CHECK( 	MPI_Dims_create(size,3,(int*) &gridSize ) );
    //		for(int d=0;d<3;d++) numberProcesses[d] = gridSize[d];
    //		std::cout << "nP: " << numberProcesses << std::endl;
    //#endif

    tarch::la::Vector<3, double> globalMDDomainSize =
        _marDynMDSolverInterface->getGlobalMDDomainSize();
    tarch::la::Vector<3, double> globalMDDomainOffset =
        _marDynMDSolverInterface->getGlobalMDDomainOffset();
    tarch::la::Vector<3, double> macroscopicCellSize =
        mamicoConfig.getMacroscopicCellConfiguration().getMacroscopicCellSize();
    //tarch::la::Vector<3,unsigned int> linkedCellsPerMacroscopicCell =
    //mamicoConfig.getMacroscopicCellConfiguration().getNumberLinkedCellsPerMacroscopicCell();

    //tarch::la::Vector<3,unsigned int> nP(1);
    DummySolverInterfaceService::getInstance().init(
        numberProcesses, rank, globalMDDomainSize, globalMDDomainOffset,
        macroscopicCellSize);
    std::cout << "DummySolverInterfaceService set!" << std::endl;

    //initialize macroscopic cell service
    const unsigned int numberMDTimestepsPerCouplingCycle = 250000;

    if (_macroscopicCellService != NULL) {
      delete _macroscopicCellService;
      _macroscopicCellService = NULL;
    }
    _macroscopicCellService =
        new coupling::services::MacroscopicCellServiceImpl<MarDynCell, 3>(
            0, _marDynMDSolverInterface,
            DummySolverInterfaceService::getInstance().getInterface(),
            numberProcesses, rank,
            mamicoConfig.getParticleInsertionConfiguration(),
            mamicoConfig.getMomentumInsertionConfiguration(),
            mamicoConfig.getBoundaryForceConfiguration(),
            mamicoConfig.getTransferStrategyConfiguration(),
            mamicoConfig.getParallelTopologyConfiguration(),
            numberMDTimestepsPerCouplingCycle,
            mamicoConfig.getMacroscopicCellConfiguration());
    //set macroscopic cell service in interface provider
    coupling::interface::MamicoInterfaceProvider<MarDynCell, 3>::getInstance()
        .setMacroscopicCellService(_macroscopicCellService);
    //set temperature for the coupled simulation
    _macroscopicCellService->computeAndStoreTemperature(1.8);

#ifdef ENABLE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    std::cout << "MacroscopicCellService set!" << std::endl;
  }

  //shut down simulation and interfaces
  void shutdown() {
    if (_marDynMDSolverInterface != NULL) {
      delete _marDynMDSolverInterface;
      _marDynMDSolverInterface = NULL;
    }
    if (_dummySolver != NULL) {
      delete _dummySolver;
      _dummySolver = NULL;
    }
    DummySolverInterfaceService::getInstance().shutdown();
    delete global_log;
  }

  //store dummy solver information in send buffers
  void storeDummySolverDataInSendBuffer() {
    tarch::la::Vector<3, unsigned int> loop(2);
    unsigned int sendCounter = 0;
    for (loop[2] = 2; loop[2] < 16; loop[2]++) {
      for (loop[1] = 2; loop[1] < 16; loop[1]++) {
        for (loop[0] = 2; loop[0] < 16; loop[0]++) {
          const tarch::la::Vector<3, unsigned int> index =
              coupling::initDimVector<3>(loop);
          const double density =
              _dummySolver->getDensity(index[0], index[1], index[2]);
          const tarch::la::Vector<3, double> velocity =
              _dummySolver->getVelocity(index[0], index[1], index[2]);

          bool send =
              DummySolverInterfaceService::getInstance().addToSendBuffer(
                  density, velocity, index);
          if (send == true) {
            sendCounter++;
          }
        }
      }
    }
    std::cout << "send counter: " << sendCounter << std::endl;
  }

  //write data from the receive buffer to the dummy solver
  void writeReceiveBufferDataToDummySolver() {
    tarch::la::Vector<3, unsigned int> loop(2);
    unsigned int receiveCounter = 0;
    for (loop[2] = 2; loop[2] < 16; loop[2]++) {
      for (loop[1] = 2; loop[1] < 16; loop[1]++) {
        for (loop[0] = 2; loop[0] < 16; loop[0]++) {
          const tarch::la::Vector<3, unsigned int> index =
              coupling::initDimVector<3>(loop);
          double density;
          tarch::la::Vector<3, double> velocity;
          bool receive =
              DummySolverInterfaceService::getInstance().getFromReceiveBuffer(
                  density, velocity, index);
          if (receive == true) {
            receiveCounter++;
            _dummySolver->setDensity(density, index[0], index[1], index[2]);
            _dummySolver->setVelocity(velocity, index[0], index[1], index[2]);
          }
        }
      }
    }
    std::cout << "receive counter: " << receiveCounter << std::endl;
  }

  //write dummy solver information to a vtk file
  void writeDummySolverInfo(unsigned int timestep) {
    std::stringstream momentum;

    momentum << "# vtk DataFile Version 2.0" << std::endl;
    momentum << "generated by MaMiCo/MarDyn-Coupling" << std::endl;
    momentum << "ASCII" << std::endl << std::endl;

    momentum << "DATASET STRUCTURED_POINTS" << std::endl;
    momentum << "DIMENSIONS 18 18 18" << std::endl;
    momentum << "ORIGIN -1.25 -1.25 -1.25" << std::endl;
    momentum << "SPACING 2.5 2.5 2.5" << std::endl;
    momentum << "POINT_DATA 5832" << std::endl << std::endl;

    momentum << "VECTORS dummyMomentum float" << std::endl;

    tarch::la::Vector<3, unsigned int> loop(0);
    for (loop[2] = 0; loop[2] < 18; loop[2]++) {
      for (loop[1] = 0; loop[1] < 18; loop[1]++) {
        for (loop[0] = 0; loop[0] < 18; loop[0]++) {
          tarch::la::Vector<3, double> vel =
              _dummySolver->getVelocity(loop[0], loop[1], loop[2]);
          momentum << vel[0] << " " << vel[1] << " " << vel[2] << std::endl;
        }
      }
    }

    std::stringstream filename;
    filename << "dummySolverMomentum_" << timestep << ".vtk";

    std::ofstream output;
    output.open(filename.str().c_str());
    if (!output.is_open()) {
      std::cout << "ERROR writeDummySolverInfo: could not open file '"
                << filename.str() << "'!" << std::endl;
      exit(EXIT_FAILURE);
    }

    output << momentum.str() << std::endl;
    momentum.clear();
    momentum.str("");
  }

  //computes the mean input velocity for the 12 overlapping cells (3-14)
  void meanInputVelocity() {
    std::cout.precision(4);

    double meanInputVelocity = 0.0;

    tarch::la::Vector<3, unsigned int> loop(0);
    loop[0] = 4;
    loop[2] = 9;
    //sum up input velocities of the 12 cells
    for (loop[1] = 3; loop[1] < 15; loop[1]++) {
      tarch::la::Vector<3, double> velocity =
          _dummySolver->getVelocity(loop[0], loop[1], loop[2]);
      std::cout << loop[1] << ": " << velocity[0] << std::endl;
      meanInputVelocity += velocity[0];
    }

    //average of 12 cells
    meanInputVelocity = meanInputVelocity / 12;

    std::cout << "computed mean input velocity is: " << meanInputVelocity
              << std::endl;
  }

  //output for velocity and density fields of the dummy solver
  void outputVelocityAndDensityFields() {
    std::cout.precision(4);

    std::cout << "Velocity field (x-direction) of dummy solver: " << std::endl;
    tarch::la::Vector<3, unsigned int> loop(0);
    loop[2] = 9;
    for (loop[1] = 2; loop[1] < 16; loop[1]++) { //4-14
      std::cout << loop[1] << "| ";
      for (loop[0] = 2; loop[0] < 16; loop[0]++) {
        tarch::la::Vector<3, double> velocity =
            _dummySolver->getVelocity(loop[0], loop[1], loop[2]);
        std::cout << velocity[0] << "| ";
      }
      std::cout << std::endl;
    }

    std::cout << "Density field of dummy solver: " << std::endl;
    tarch::la::Vector<3, unsigned int> loopD(0);
    loopD[2] = 9;
    for (loopD[1] = 4; loopD[1] < 14; loopD[1]++) {
      std::cout << loopD[1] << "\t | ";
      for (loopD[0] = 4; loopD[0] < 14; loopD[0]++) {
        std::cout << _dummySolver->getDensity(loopD[0], loopD[1], loopD[2])
                  << "\t | ";
      }
      std::cout << std::endl;
    }

  }

  MarDynMDSolverInterface *_marDynMDSolverInterface;
  DummySolver *_dummySolver;
  coupling::services::MacroscopicCellService<3> *_macroscopicCellService;

  int _argc;
  char **_argv;
};

#endif /* TESTMARDYN2DUMMYCOUPLING_H_ */
