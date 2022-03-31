// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_TESTS_NIETEST_H_
#define _MOLECULARDYNAMICS_COUPLING_TESTS_NIETEST_H_

#include "coupling/CouplingMDDefinitions.h"
#include "coupling/configurations/MaMiCoConfiguration.h"
#include "coupling/interface/MDSimulationFactory.h"
#include "coupling/interface/impl/SimpleMD/SimpleMDSolverInterface.h"
#include "coupling/services/MultiMDCellService.h"
#include "coupling/solvers/CouetteSolver.h"
#include "coupling/solvers/CouetteSolverInterface.h"
#include "coupling/solvers/LBCouetteSolver.h"
#include "coupling/solvers/LBCouetteSolverInterface.h"
#include "coupling/tests/Test.h"
#include "simplemd/configurations/MolecularDynamicsConfiguration.h"
#include "tarch/configuration/ParseConfiguration.h"
#include "tarch/utils/MultiMDService.h"
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
#include <mpi.h>
#endif
#include <sys/time.h>

/** tests couette flow for Nie coupling scheme. We therefore establish one-way
 * coupling to an analytical Couette flow solver.
 *  @author Philipp Neumann
 */
class NieTest : public Test {
private:
  /** choose between analytical and LB Couette solver */
  enum SolverType { COUETTE_ANALYTICAL = 0, COUETTE_LB = 1 };

public:
  NieTest() : Test("NieTest") {}
  virtual ~NieTest() {}

  virtual void run() {
    int rank = 0;
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size); // for debugging only
#endif

    const double density = 0.81;           // density of particle system; DPD value: 3.0
    const double kinVisc = 2.14 / density; // kinematic viscosity (cf. paper by Nie et al.); DPD
                                           // value: 0.54/density
    const double temp = 1.1;               // temperature of MD system; DPD value: 1.0
    const int equSteps = 0;                // number of equilibration steps ; DPD value: 1000, MD value: 10000
    unsigned int mdStepCounter = 0;        // time step counter for MD

    const double channelheight = 50.0;                              // channel is always expected to have origin at (0.0,0.0,0.0) and
                                                                    // to be cubic (MD 30: 50.0, MD 60: 100.0, MD 120: 200.0)
    const tarch::la::Vector<3, double> wallVelocity(0.5, 0.0, 0.0); // velocity of moving wall (lower boundary moves);
                                                                    // analytic solver only supports flow in x-direction
    const int plotEveryTimestep = 1;                                // only for LB couette solver: VTK plotting per time step
    tarch::la::Vector<3, unsigned int> lbNumberProcesses(1, 1, 1);  // only for LB couette solver: number of processes

    const unsigned int couplingCycles = 105;         // number of coupling cycles, that is continuum time steps; MD/DPD:
                                                     // 1000
    const unsigned int totalNumberMDSimulations = 1; // total number of MD simulations; MD/DPD: 64 (DPD), 128 (MD
                                                     // scalability), 144 (MD)
    const SolverType solverType = COUETTE_LB; // LB or analytical couette solver
    const bool twoWayCoupling = true;
    const int twoWayCouplingInitCycles = 0;

    // for time measurements
    timeval start;
    timeval end;

    // parse configurations and start time measurement for initialisation
    simplemd::configurations::MolecularDynamicsConfiguration simpleMDConfig;
    tarch::configuration::ParseConfiguration::parseConfiguration<simplemd::configurations::MolecularDynamicsConfiguration>(
        "test_nie_simplemd.xml", "molecular-dynamics", simpleMDConfig);
    if (!simpleMDConfig.isValid()) {
      std::cout << "ERROR NieTest: Invalid SimpleMD config!" << std::endl;
      exit(EXIT_FAILURE);
    }
    coupling::configurations::MaMiCoConfiguration<3> mamicoConfig;
    tarch::configuration::ParseConfiguration::parseConfiguration<coupling::configurations::MaMiCoConfiguration<3>>("test_nie.xml", "mamico", mamicoConfig);
    if (!mamicoConfig.isValid()) {
      std::cout << "ERROR NieTest: Invalid MaMiCo config!" << std::endl;
      exit(EXIT_FAILURE);
    }

    // allocate solvers
    coupling::solvers::AbstractCouetteSolver<3> *couetteSolver = NULL;
    couetteSolver = getCouetteSolver(channelheight, tarch::la::Vector<3, double>(0.0, 0.0, 0.0), kinVisc,
                                     mamicoConfig.getMacroscopicCellConfiguration().getMacroscopicCellSize()[0],
                                     simpleMDConfig.getSimulationConfiguration().getDt() * simpleMDConfig.getSimulationConfiguration().getNumberOfTimesteps(),
                                     plotEveryTimestep, "LBCouette", lbNumberProcesses, solverType, rank);

    coupling::solvers::AbstractCouetteSolver<3> *couetteSolver2 = NULL;
    couetteSolver2 = getCouetteSolver(channelheight, wallVelocity, kinVisc, mamicoConfig.getMacroscopicCellConfiguration().getMacroscopicCellSize()[0],
                                      simpleMDConfig.getSimulationConfiguration().getDt() * simpleMDConfig.getSimulationConfiguration().getNumberOfTimesteps(),
                                      0, "Couette2", lbNumberProcesses, COUETTE_ANALYTICAL, rank);

    tarch::utils::MultiMDService<3> multiMDService(simpleMDConfig.getMPIConfiguration().getNumberOfProcesses(), totalNumberMDSimulations);
    const unsigned int mdInstances = multiMDService.getLocalNumberOfMDSimulations(); // local number of MD simulations
    // std::cout << "Rank: " << rank << " # total MD instances: " <<
    // totalNumberMDSimulations << ", # local MD instances: " << mdInstances <<
    // std::endl;

    std::vector<coupling::interface::MDSimulation *> simpleMD;
    for (unsigned int i = 0; i < mdInstances; i++) {
      simpleMD.push_back(coupling::interface::SimulationAndInterfaceFactory::getInstance().getMDSimulation(simpleMDConfig, mamicoConfig
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
                                                                                                           ,
                                                                                                           multiMDService.getLocalCommunicator()
#endif
                                                                                                               ));
      if (simpleMD[i] == NULL) {
        std::cout << "ERROR NieTest: simpleMD[" << i << "]==NULL!" << std::endl;
        exit(EXIT_FAILURE);
      }
    }

    // equilibrate MD
    if (rank == 0) {
      gettimeofday(&start, NULL);
    }
    for (unsigned int i = 0; i < mdInstances; i++) {
      simpleMD[i]->init(multiMDService, multiMDService.getGlobalNumberOfLocalMDSimulation(i));

      simpleMD[i]->switchOffCoupling();
      simpleMD[i]->simulateTimesteps(equSteps, mdStepCounter);
      mdStepCounter += equSteps;

      // reset counter for all but the last MD instance
      if (i < mdInstances - 1) {
        mdStepCounter = 0;
      }
    }
    // finish time measurement for initialisation(incl equilibration) and start
    // time measurement for coupling
    if (rank == 0) {
      gettimeofday(&end, NULL);
      double runtime = (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec - start.tv_usec);
      std::cout << "NieTest-Equilibration: " << (int)(runtime / 1000) << "ms" << std::endl;
    }

    // allocate coupling interfaces
    std::vector<coupling::interface::MDSolverInterface<MY_LINKEDCELL, 3> *> mdSolverInterface;
    for (unsigned int i = 0; i < mdInstances; i++) {
      // switch on coupling; in case of LAMMPS: initialse MDSolverInterface and
      // buffer result in MamicoInterfaceProvider; it will be delivered to
      // mdSolverInterface from there
      simpleMD[i]->switchOnCoupling();
      mdSolverInterface.push_back(
          coupling::interface::SimulationAndInterfaceFactory::getInstance().getMDSolverInterface(simpleMDConfig, mamicoConfig, simpleMD[i]));
      if (mdSolverInterface[i] == NULL) {
        std::cout << "ERROR main_nie: mdSolverInterface[" << i << "]==NULL!" << std::endl;
        exit(EXIT_FAILURE);
      }
    }
    coupling::interface::MacroscopicSolverInterface<3> *couetteSolverInterface = getCouetteSolverInterface(
        couetteSolver, simpleMDConfig.getDomainConfiguration().getGlobalDomainOffset(), mamicoConfig.getMacroscopicCellConfiguration().getMacroscopicCellSize(),
        getGlobalNumberMacroscopicCells(simpleMDConfig, mamicoConfig), mamicoConfig.getMomentumInsertionConfiguration()._innerOverlap, solverType, rank);

    // initialise macroscopic cell service for multi-MD case and set single cell
    // services in each MD simulation
    coupling::services::MultiMDCellService<MY_LINKEDCELL, 3> multiMDCellService(
        mdSolverInterface, couetteSolverInterface, simpleMDConfig.getMPIConfiguration().getNumberOfProcesses(), (unsigned int)rank, totalNumberMDSimulations,
        mamicoConfig.getParticleInsertionConfiguration(), mamicoConfig.getMomentumInsertionConfiguration(), mamicoConfig.getBoundaryForceConfiguration(),
        mamicoConfig.getTransferStrategyConfiguration(), mamicoConfig.getNoiseReductionConfiguration(), mamicoConfig.getParallelTopologyConfiguration(),
        simpleMDConfig.getSimulationConfiguration().getNumberOfTimesteps(), mamicoConfig.getMacroscopicCellConfiguration(), multiMDService);
    // set couette solver interface in MamicoInterfaceProvider
    coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL, 3>::getInstance().setMacroscopicSolverInterface(couetteSolverInterface);
    if ((solverType == COUETTE_LB) && twoWayCoupling) {
      static_cast<coupling::solvers::LBCouetteSolver *>(couetteSolver)
          ->setMDBoundary(simpleMDConfig.getDomainConfiguration().getGlobalDomainOffset(), simpleMDConfig.getDomainConfiguration().getGlobalDomainSize(),
                          mamicoConfig.getMomentumInsertionConfiguration()._innerOverlap);
    }

    for (unsigned int i = 0; i < mdInstances; i++) {
      simpleMD[i]->setMacroscopicCellService(&(multiMDCellService.getMacroscopicCellService(i)));
      // compute and store temperature in macroscopic cells (temp=1.1
      // everywhere)
      multiMDCellService.getMacroscopicCellService(i).computeAndStoreTemperature(temp);
    }

    // allocate buffers for send/recv operations
    std::vector<coupling::datastructures::MacroscopicCell<3> *> sendBuffer;
    unsigned int *globalCellIndices4SendBuffer =
        allocateSendBuffer(sendBuffer, multiMDCellService.getMacroscopicCellService(0).getIndexConversion(), rank, *couetteSolverInterface);
    std::vector<coupling::datastructures::MacroscopicCell<3> *> recvBuffer;
    unsigned int *globalCellIndices4RecvBuffer =
        allocateRecvBuffer(recvBuffer, multiMDCellService.getMacroscopicCellService(0).getIndexConversion(), rank, *couetteSolverInterface);

    std::default_random_engine generator(1234);

    // create noise reduction for analytical solver
    coupling::noisereduction::NoiseReduction<3> *noiseReduction(mamicoConfig.getNoiseReductionConfiguration().interpreteConfiguration<3>(
        multiMDCellService.getMacroscopicCellService(0).getIndexConversion(), multiMDService));

    // coupling
    if (rank == 0) {
      gettimeofday(&start, NULL);
    }
    for (unsigned int cycles = 0; cycles < couplingCycles; cycles++) {
      // run one time step for couette solver
      if (couetteSolver != NULL) {
        couetteSolver->advance(simpleMDConfig.getSimulationConfiguration().getDt() * simpleMDConfig.getSimulationConfiguration().getNumberOfTimesteps());

        couetteSolver2->advance(simpleMDConfig.getSimulationConfiguration().getDt() * simpleMDConfig.getSimulationConfiguration().getNumberOfTimesteps());
      }
      // // extract data from couette solver and send them to MD (can take any
      // index-conversion object)
      // fillSendBuffer(density,*couetteSolver,multiMDCellService.getMacroscopicCellService(0).getIndexConversion(),sendBuffer,globalCellIndices4SendBuffer);
      // multiMDCellService.sendFromMacro2MD(sendBuffer,globalCellIndices4SendBuffer);
      // // run MD instances
      // for (unsigned int i = 0; i < mdInstances; i++){
      //   // set macroscopic cell service and interfaces in
      // MamicoInterfaceProvider
      //   coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL,3>::getInstance().setMacroscopicCellService(&(multiMDCellService.getMacroscopicCellService(i)));
      //   coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL,3>::getInstance().setMDSolverInterface(mdSolverInterface[i]);

      //   simpleMD[i]->simulateTimesteps(simpleMDConfig.getSimulationConfiguration().getNumberOfTimesteps(),mdStepCounter);
      //   mdStepCounter+=simpleMDConfig.getSimulationConfiguration().getNumberOfTimesteps();
      //   // reset mdStepCounter unless this is the last MD instance
      //   if (i<mdInstances-1){mdStepCounter -=
      // simpleMDConfig.getSimulationConfiguration().getNumberOfTimesteps();}

      //   // plot macroscopic time step info in multi md service
      //   multiMDCellService.getMacroscopicCellService(i).plotEveryMacroscopicTimestep(cycles);
      // }
      // // send back data from MD instances and merge it
      // multiMDCellService.sendFromMD2Macro(recvBuffer,globalCellIndices4RecvBuffer);

      fillRecvBuffer(density, *couetteSolver2, multiMDCellService.getMacroscopicCellService(0).getIndexConversion(), recvBuffer, globalCellIndices4RecvBuffer,
                     generator);

      // call noise filter on recvBuffer
      noiseReduction->beginProcessInnerMacroscopicCells();
      for (unsigned int i = 0; i < recvBuffer.size(); i++) {
        noiseReduction->processInnerMacroscopicCell(*recvBuffer[i], i);
      }
      noiseReduction->endProcessInnerMacroscopicCells();
      if (noiseReduction->_doubleTraversal) {
        noiseReduction->beginProcessInnerMacroscopicCells();
        for (unsigned int i = 0; i < recvBuffer.size(); i++) {
          noiseReduction->processInnerMacroscopicCell(*recvBuffer[i], i);
        }
        noiseReduction->endProcessInnerMacroscopicCells();
      }

      if ((solverType == COUETTE_LB) && twoWayCoupling && cycles >= twoWayCouplingInitCycles) {
        static_cast<coupling::solvers::LBCouetteSolver *>(couetteSolver)
            ->setMDBoundaryValues(recvBuffer, globalCellIndices4RecvBuffer, multiMDCellService.getMacroscopicCellService(0).getIndexConversion());
      }
      // write data to csv-compatible file for evaluation
      // write2CSV(recvBuffer,globalCellIndices4RecvBuffer,multiMDCellService.getMacroscopicCellService(0).getIndexConversion(),rank,cycles);
      if (rank == 0) {
        std::cout << "Finish coupling cycle " << cycles << std::endl;
      }
    }

    // finish time measurement for coupled simulation; start time measurement
    // for shut down of simulation
    if (rank == 0) {
      gettimeofday(&end, NULL);
      double runtime = (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec - start.tv_usec);
      std::cout << "NieTest-Coupling: " << (int)(runtime / 1000) << "ms for " << couplingCycles << " coupling cycles and " << totalNumberMDSimulations
                << " MD simulations" << std::endl;
      gettimeofday(&start, NULL);
    }

    // shutdown simulation and free buffers/arrays
    deleteBuffer(sendBuffer);
    if (globalCellIndices4SendBuffer != NULL) {
      delete[] globalCellIndices4SendBuffer;
      globalCellIndices4SendBuffer = NULL;
    }
    deleteBuffer(recvBuffer);
    if (globalCellIndices4RecvBuffer != NULL) {
      delete[] globalCellIndices4RecvBuffer;
      globalCellIndices4RecvBuffer = NULL;
    }
    for (unsigned int i = 0; i < mdInstances; i++) {
      // the shutdown operation may also delete the md solver interface;
      // therefore, we update the MD solver interface in the vector
      // mdSolverInteface after the shutdown is completed
      coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL, 3>::getInstance().setMDSolverInterface(mdSolverInterface[i]);
      simpleMD[i]->shutdown();
      delete simpleMD[i];
      simpleMD[i] = NULL;
      mdSolverInterface[i] = coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL, 3>::getInstance().getMDSolverInterface();
    }
    simpleMD.clear();
    for (unsigned int i = 0; i < mdInstances; i++) {
      if (mdSolverInterface[i] != NULL) {
        delete mdSolverInterface[i];
        mdSolverInterface[i] = NULL;
      }
    }
    mdSolverInterface.clear();
    if (couetteSolverInterface != NULL) {
      delete couetteSolverInterface;
      couetteSolverInterface = NULL;
    }
    if (couetteSolver != NULL) {
      delete couetteSolver;
      couetteSolver = NULL;
    }

    // end time measurement for shutdown
    if (rank == 0) {
      gettimeofday(&end, NULL);
      double runtime = (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec - start.tv_usec);
      std::cout << "NieTest-Shutdown: " << (int)(runtime / 1000) << "ms" << std::endl;
    }
  }

private:
  /** computes global number of macroscopic cells from configs. Required by
   * couette solver interface before MacroscopicCellService is initialised! */
  tarch::la::Vector<3, unsigned int> getGlobalNumberMacroscopicCells(const simplemd::configurations::MolecularDynamicsConfiguration &simpleMDConfig,
                                                                     const coupling::configurations::MaMiCoConfiguration<3> &mamicoConfig) const {
    tarch::la::Vector<3, double> domainSize(simpleMDConfig.getDomainConfiguration().getGlobalDomainSize());
    tarch::la::Vector<3, double> dx(mamicoConfig.getMacroscopicCellConfiguration().getMacroscopicCellSize());
    tarch::la::Vector<3, unsigned int> globalNumberMacroscopicCells(0);
    for (unsigned int d = 0; d < 3; d++) {
      int buf = floor(domainSize[d] / dx[d] + 0.5);
      globalNumberMacroscopicCells[d] = (unsigned int)buf;
    }
    return globalNumberMacroscopicCells;
  }

  /** allocates the send buffer (with values for all macroscopic cells) and
   * returns indices. This is only done on rank 0. */
  unsigned int *allocateSendBuffer(std::vector<coupling::datastructures::MacroscopicCell<3> *> &sendBuffer, const coupling::IndexConversion<3> &indexConversion,
                                   int rank, coupling::interface::MacroscopicSolverInterface<3> &couetteSolverInterface) const {
    // determine global number of cells
    const tarch::la::Vector<3, unsigned int> cells(indexConversion.getGlobalNumberMacroscopicCells() + tarch::la::Vector<3, unsigned int>(2));
    const unsigned int num = cells[0] * cells[1] * cells[2];

    // delete all potential entries of sendBuffer
    deleteBuffer(sendBuffer);

    // count number of cells to be sent from this process; therefore, loop over
    // all global macroscopic cells...
    unsigned int numCellsSent = 0;
    for (unsigned int i = 0; i < num; i++) {
      // ... and find out, if the current cell should be send to MD from this
      // couette solver process
      if (couetteSolverInterface.sendMacroscopicQuantityToMDSolver(indexConversion.getGlobalVectorCellIndex(i))) {
        std::vector<unsigned int> ranks = couetteSolverInterface.getSourceRanks(indexConversion.getGlobalVectorCellIndex(i));
        bool containsThisRank = false;
        for (unsigned int k = 0; k < ranks.size(); k++) {
          containsThisRank = containsThisRank || (ranks[k] == (unsigned int)rank);
        }
        if (containsThisRank) {
          numCellsSent++;
        }
      }
    }

    // allocate array for cell indices
    unsigned int *indices = new unsigned int[numCellsSent];
    if (indices == NULL) {
      std::cout << "ERROR NieTest::allocateSendBuffer(): indices==NULL!" << std::endl;
      exit(EXIT_FAILURE);
    }

    // allocate sendBuffer and initialise all entries, incl. indices
    for (unsigned int i = 0; i < num; i++) {
      if (couetteSolverInterface.sendMacroscopicQuantityToMDSolver(indexConversion.getGlobalVectorCellIndex(i))) {
        std::vector<unsigned int> ranks = couetteSolverInterface.getSourceRanks(indexConversion.getGlobalVectorCellIndex(i));
        bool containsThisRank = false;
        for (unsigned int k = 0; k < ranks.size(); k++) {
          containsThisRank = containsThisRank || (ranks[k] == (unsigned int)rank);
        }
        if (containsThisRank) {
          sendBuffer.push_back(new coupling::datastructures::MacroscopicCell<3>());
          if (sendBuffer[sendBuffer.size() - 1] == NULL) {
            std::cout << "ERROR NieTest::allocateSendBuffer: sendBuffer[" << sendBuffer.size() - 1 << "]==NULL!" << std::endl;
            exit(EXIT_FAILURE);
          }
          indices[sendBuffer.size() - 1] = i;
        }
      }
    }

#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
    for (unsigned int i = 0; i < numCellsSent; i++) {
      std::vector<unsigned int> ranks = couetteSolverInterface.getSourceRanks(indexConversion.getGlobalVectorCellIndex(indices[i]));
      std::cout << "Current rank= " << rank << ", Send cell " << indexConversion.getGlobalVectorCellIndex(indices[i]) << ", ranks=";
      for (unsigned int j = 0; j < ranks.size(); j++) {
        std::cout << " " << ranks[j];
      }
      std::cout << std::endl;
    }
#endif
    return indices;
  }

  /** allocates the recv-buffer. This buffer contains all global inner
   * macroscopic cells, but only on rank 0. On all other ranks, no cells are
   * stored and a NULL ptr is returned */
  unsigned int *allocateRecvBuffer(std::vector<coupling::datastructures::MacroscopicCell<3> *> &recvBuffer, const coupling::IndexConversion<3> &indexConversion,
                                   int rank, coupling::interface::MacroscopicSolverInterface<3> &couetteSolverInterface) const {

    // determine global number of cells
    const tarch::la::Vector<3, unsigned int> cells(indexConversion.getGlobalNumberMacroscopicCells() + tarch::la::Vector<3, unsigned int>(2));
    const unsigned int num = cells[0] * cells[1] * cells[2];

    // delete all potential entries of sendBuffer
    deleteBuffer(recvBuffer);

    // determine number of cells that should be received
    unsigned int numCellsRecv = 0;
    for (unsigned int i = 0; i < num; i++) {
      if (couetteSolverInterface.receiveMacroscopicQuantityFromMDSolver(indexConversion.getGlobalVectorCellIndex(i))) {
        std::vector<unsigned int> ranks = couetteSolverInterface.getTargetRanks(indexConversion.getGlobalVectorCellIndex(i));
        bool containsThisRank = false;
        for (unsigned int k = 0; k < ranks.size(); k++) {
          containsThisRank = containsThisRank || (ranks[k] == (unsigned int)rank);
        }
        if (containsThisRank) {
          numCellsRecv++;
        }
      }
    }
    // allocate array for cell indices
    unsigned int *indices = new unsigned int[numCellsRecv];
    if (indices == NULL) {
      std::cout << "ERROR NieTest::allocateRecvBuffer(): indices==NULL!" << std::endl;
      exit(EXIT_FAILURE);
    }

    // allocate recvBuffer and initialise all entries, incl. indices
    for (unsigned int i = 0; i < num; i++) {
      if (couetteSolverInterface.receiveMacroscopicQuantityFromMDSolver(indexConversion.getGlobalVectorCellIndex(i))) {
        std::vector<unsigned int> ranks = couetteSolverInterface.getTargetRanks(indexConversion.getGlobalVectorCellIndex(i));
        bool containsThisRank = false;
        for (unsigned int k = 0; k < ranks.size(); k++) {
          containsThisRank = containsThisRank || (ranks[k] == (unsigned int)rank);
        }
        if (containsThisRank) {
          recvBuffer.push_back(new coupling::datastructures::MacroscopicCell<3>());
          if (recvBuffer[recvBuffer.size() - 1] == NULL) {
            std::cout << "ERROR NieTest::allocateRecvBuffer: recvBuffer[" << recvBuffer.size() - 1 << "]==NULL!" << std::endl;
            exit(EXIT_FAILURE);
          }
          // set linearized index
          indices[recvBuffer.size() - 1] = i;
        }
      }
    }

#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
    for (unsigned int i = 0; i < numCellsRecv; i++) {
      std::vector<unsigned int> ranks = couetteSolverInterface.getTargetRanks(indexConversion.getGlobalVectorCellIndex(indices[i]));
      std::cout << "Current rank= " << rank << ", Recv cell " << indexConversion.getGlobalVectorCellIndex(indices[i]) << ", ranks=";
      for (unsigned int j = 0; j < ranks.size(); j++) {
        std::cout << " " << ranks[j];
      }
      std::cout << std::endl;
    }
#endif
    return indices;
  }

  /** write cells that have been received from MD to csv file */
  void write2CSV(std::vector<coupling::datastructures::MacroscopicCell<3> *> &recvBuffer, const unsigned int *const recvIndices,
                 const coupling::IndexConversion<3> &indexConversion, int rank, int couplingCycle) const {
    // form file name and open file
    std::stringstream ss;
    ss << "CouetteAvgMultiMDCells_" << rank << "_" << couplingCycle << ".csv";
    std::ofstream file(ss.str().c_str());
    if (!file.is_open()) {
      std::cout << "ERROR NieTest::write2CSV(): Could not open file " << ss.str() << "!" << std::endl;
      exit(EXIT_FAILURE);
    }

    // loop over received cells; read macroscopic mass+momentum buffers and
    // write cell index, mass and velocity to one line in the csv-file
    const unsigned int numCellsRecv = recvBuffer.size();
    for (unsigned int i = 0; i < numCellsRecv; i++) {
      tarch::la::Vector<3, double> vel(recvBuffer[i]->getMacroscopicMomentum());
      if (recvBuffer[i]->getMacroscopicMass() != 0.0) {
        vel = (1.0 / recvBuffer[i]->getMacroscopicMass()) * vel;
      }
      const tarch::la::Vector<3, unsigned int> counter(indexConversion.getGlobalVectorCellIndex(recvIndices[i]));
      file << counter[0] << " ; " << counter[1] << " ; " << counter[2] << " ; " << vel[0] << " ; " << vel[1] << " ; " << vel[2] << " ; "
           << recvBuffer[i]->getMacroscopicMass() << ";";
      file << std::endl;
    }

    // close file
    file.close();
  }

  /** deletes the send buffer */
  void deleteBuffer(std::vector<coupling::datastructures::MacroscopicCell<3> *> &sendBuffer) const {
    // delete all potential entries of sendBuffer
    for (unsigned int i = 0; i < sendBuffer.size(); i++) {
      if (sendBuffer[i] != NULL) {
        delete sendBuffer[i];
        sendBuffer[i] = NULL;
      }
    }
    sendBuffer.clear();
  }

  /** fills send buffer with data from couette solver */
  void fillSendBuffer(const double density, const coupling::solvers::AbstractCouetteSolver<3> &couetteSolver,
                      const coupling::IndexConversion<3> &indexConversion, std::vector<coupling::datastructures::MacroscopicCell<3> *> &sendBuffer,
                      const unsigned int *const globalCellIndices4SendBuffer) const {
    const unsigned int size = sendBuffer.size();
    const tarch::la::Vector<3, double> domainOffset(indexConversion.getGlobalMDDomainOffset());
    const tarch::la::Vector<3, double> macroscopicCellSize(indexConversion.getMacroscopicCellSize());
    const double mass = density * macroscopicCellSize[0] * macroscopicCellSize[1] * macroscopicCellSize[2];

    for (unsigned int i = 0; i < size; i++) {
      // get global cell index vector
      const tarch::la::Vector<3, unsigned int> globalIndex(indexConversion.getGlobalVectorCellIndex(globalCellIndices4SendBuffer[i]));
      // determine cell midpoint
      tarch::la::Vector<3, double> cellMidPoint(domainOffset - 0.5 * macroscopicCellSize);
      for (unsigned int d = 0; d < 3; d++) {
        cellMidPoint[d] = cellMidPoint[d] + ((double)globalIndex[d]) * macroscopicCellSize[d];
      }
      // compute momentum
      const tarch::la::Vector<3, double> momentum(mass * couetteSolver.getVelocity(cellMidPoint));
      sendBuffer[i]->setMicroscopicMass(mass);
      sendBuffer[i]->setMicroscopicMomentum(momentum);
    }
  }

  void fillRecvBuffer(const double density, const coupling::solvers::AbstractCouetteSolver<3> &couetteSolver,
                      const coupling::IndexConversion<3> &indexConversion, std::vector<coupling::datastructures::MacroscopicCell<3> *> &sendBuffer,
                      const unsigned int *const globalCellIndices4SendBuffer, std::default_random_engine &generator) const {
    const unsigned int size = sendBuffer.size();
    const tarch::la::Vector<3, double> domainOffset(indexConversion.getGlobalMDDomainOffset());
    const tarch::la::Vector<3, double> macroscopicCellSize(indexConversion.getMacroscopicCellSize());
    const double mass = density * macroscopicCellSize[0] * macroscopicCellSize[1] * macroscopicCellSize[2];

    std::normal_distribution<double> distribution(0.0, 1.0 / 6.0);
    for (unsigned int i = 0; i < size; i++) {
      // get global cell index vector
      const tarch::la::Vector<3, unsigned int> globalIndex(indexConversion.getGlobalVectorCellIndex(globalCellIndices4SendBuffer[i]));
      // determine cell midpoint
      tarch::la::Vector<3, double> cellMidPoint(domainOffset - 0.5 * macroscopicCellSize);
      for (unsigned int d = 0; d < 3; d++) {
        cellMidPoint[d] = cellMidPoint[d] + ((double)globalIndex[d]) * macroscopicCellSize[d];
      }
      // compute momentum
      const tarch::la::Vector<3, double> noise(distribution(generator), distribution(generator), distribution(generator));
      const tarch::la::Vector<3, double> momentum(mass * (couetteSolver.getVelocity(cellMidPoint) /*+noise*/));
      sendBuffer[i]->setMacroscopicMass(mass);
      sendBuffer[i]->setMacroscopicMomentum(momentum);
    }
  }

  coupling::solvers::AbstractCouetteSolver<3> *getCouetteSolver(const double channelheight, tarch::la::Vector<3, double> wallVelocity, const double kinVisc,
                                                                const double dx, const double dt, const int plotEveryTimestep, const std::string filestem,
                                                                const tarch::la::Vector<3, unsigned int> processes, const SolverType solverType,
                                                                const int rank) {
    coupling::solvers::AbstractCouetteSolver<3> *solver = NULL;
    // analytical solver: is only active on rank 0
    if (solverType == COUETTE_ANALYTICAL) {
      if (rank == 0) {
        solver = new coupling::solvers::CouetteSolver<3>(channelheight, wallVelocity[0], kinVisc);
        if (solver == NULL) {
          std::cout << "ERROR NieTest::getCouetteSolver(): Analytic solver==NULL!" << std::endl;
          exit(EXIT_FAILURE);
        }
      }
      // LB solver: active on "processes"
    } else if (solverType == COUETTE_LB) {
      solver = new coupling::solvers::LBCouetteSolver(channelheight, wallVelocity, kinVisc, dx, dt, plotEveryTimestep, filestem, processes, 1);
      if (solver == NULL) {
        std::cout << "ERROR NieTest::getCouetteSolver(): LB solver==NULL!" << std::endl;
        exit(EXIT_FAILURE);
      }
    } else {
      std::cout << "ERROR NieTest::getCouetteSolver(): Unknown solver type!" << std::endl;
      exit(EXIT_FAILURE);
      return NULL;
    }
    return solver;
  }

  coupling::interface::MacroscopicSolverInterface<3> *getCouetteSolverInterface(coupling::solvers::AbstractCouetteSolver<3> *couetteSolver,
                                                                                tarch::la::Vector<3, double> mdOffset,
                                                                                tarch::la::Vector<3, double> mamicoMeshsize,
                                                                                tarch::la::Vector<3, unsigned int> globalNumberMacroscopicCells,
                                                                                unsigned int outerRegion, const SolverType solverType, const int rank) {
    coupling::interface::MacroscopicSolverInterface<3> *interface = NULL;
    if (solverType == COUETTE_ANALYTICAL) {
      interface = new coupling::solvers::CouetteSolverInterface<3>(globalNumberMacroscopicCells, outerRegion);
    } else if (solverType == COUETTE_LB) {
      coupling::solvers::LBCouetteSolver *lbSolver = static_cast<coupling::solvers::LBCouetteSolver *>(couetteSolver);
      if (lbSolver == NULL) {
        std::cout << "ERROR NieTest::getCouetteSolverInterface(...), rank=" << rank << ": Could not convert abstract to LB solver!" << std::endl;
        exit(EXIT_FAILURE);
      }
      // compute number of cells of MD offset; detect any mismatches!
      tarch::la::Vector<3, unsigned int> offsetMDDomain(0);
      for (unsigned int d = 0; d < 3; d++) {
        if (mdOffset[d] < 0.0) {
          std::cout << "ERROR NieTest::getCouetteSolverInterface(...): mdOffset[" << d << "]<0.0!" << std::endl;
          exit(EXIT_FAILURE);
        }
        offsetMDDomain[d] = floor(mdOffset[d] / mamicoMeshsize[d] + 0.5);
        if (fabs((offsetMDDomain[d] * mamicoMeshsize[d] - mdOffset[d]) / mamicoMeshsize[d]) > 1.0e-8) {
          std::cout << "ERROR NieTest::getCouetteSolverInterface: MD offset "
                       "and mesh size mismatch!"
                    << std::endl;
          exit(EXIT_FAILURE);
        }
      }
      interface = new coupling::solvers::LBCouetteSolverInterface(lbSolver->getAvgNumberLBCells(), lbSolver->getNumberProcesses(), offsetMDDomain,
                                                                  globalNumberMacroscopicCells, outerRegion);
    }
    if (interface == NULL) {
      std::cout << "ERROR NieTest::getCouetteSolverInterface(...), rank=" << rank << ": interface==NULL!" << std::endl;
      exit(EXIT_FAILURE);
    }
    return interface;
  }
};
#endif // _MOLECULARDYNAMICS_COUPLING_TESTS_NIETEST_H_
