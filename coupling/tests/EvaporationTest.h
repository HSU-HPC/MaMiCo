// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _COUPLING_TESTS_EVAPORATIONTEST_H_
#define _COUPLING_TESTS_EVAPORATIONTEST_H_

#include "coupling/CouplingMDDefinitions.h"
#include "coupling/tests/Test.h"
#include "simplemd/configurations/MolecularDynamicsConfiguration.h"
#include "tarch/configuration/ParseConfiguration.h"
#include "tarch/utils/MultiMDService.h"
#if (BUILD_WITH_OPENFOAM)
#include "coupling/solvers/RhoCentralInterface4Evaporation.h"
// #include "coupling/solvers/FoamSolverInterface.h"
#endif
#include "coupling/configurations/EvaporationConfiguration.h"
#include "coupling/configurations/MaMiCoConfiguration.h"
#include "coupling/indexing/IndexingService.cpp"
#include "coupling/interface/MDSimulationFactory.h"
#include "coupling/interface/impl/SimpleMD/SimpleMDSolverInterface.h"
#include "coupling/services/MultiMDCellService.h"
#include "coupling/solvers/CouetteSolverInterface.h"
#include "coupling/solvers/LBCouetteSolverInterface.h"
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
#include <mpi.h>
#endif
#include <random>
#include <sys/time.h>

/**
 * @brief Runs a evaporation scenario simulation
 * @author Helene Wittenberg
 */
class EvaporationTest : public Test {
public:
  /** @brief simple constructor */
  EvaporationTest() : Test("EvaporationTest") {}
  /** @brief a dummy destructor */
  virtual ~EvaporationTest() {}

  /** triggers void init(), runOneCouplingCycle() and shutdown()
   *  @brief runs the simulation */
  virtual void run() {
    init();
    for (int cycle = 0; cycle < _EvapConfig.couplingCycles; cycle++) {
      runOneCouplingCycle(cycle);
    }
    shutdown();
  }

private:
  /** triggers initMPI(), parseConfiguration(), and initSolvers()
   *  @brief initialises everthing necessary for the test */
  void init() {
    initMPI();
    parseConfigurations();
    initSolvers();
  }

  /** it advances the macro (advanceMacro()) and micro solver (advanceMicro),
   *  computes the signal to noise ratio (computeSNR()) and sends the data from the macro to
   *  the micro solver (twoWayCoupling())
   *  @brief combines the functioniality necessary for a cycle of the coupled simulation  */
  void runOneCouplingCycle(int cycle) {
    advanceMacro(cycle);
    advanceMicro(cycle);
  }

  /** @brief initialises all MPI variables  */
  void initMPI() {
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
    MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
#endif
  }

  /** @brief reads the configuration from the xml file and calls parseEvaporationTestConfiguration() */
  void parseConfigurations() {
    std::string filename("evaporation.xml");
    tarch::configuration::ParseConfiguration::parseConfiguration<simplemd::configurations::MolecularDynamicsConfiguration>(filename, "molecular-dynamics",
                                                                                                                           _MDSolverConfig);
    if (!_MDSolverConfig.isValid()) {
      std::cout << "ERROR EvaporationTest: Invalid MDSolver config!" << std::endl;
      exit(EXIT_FAILURE);
    }
    tarch::configuration::ParseConfiguration::parseConfiguration<coupling::configurations::MaMiCoConfiguration<3>>(filename, "mamico", _mamicoConfig);
    if (!_mamicoConfig.isValid()) {
      std::cout << "ERROR EvaporationTest: Invalid MaMiCo config!" << std::endl;
      exit(EXIT_FAILURE);
    }

    _EvapConfig = coupling::configurations::EvaporationConfig::parseConfiguration(filename);
  }

  /** @brief initialises the macro and micro solver according to the setup from the xml file and pre-proccses them */
  void initSolvers() {
    // allocate solvers
    _macroSolver = new coupling::solvers::RhoCentralInterface4Evaporation(_EvapConfig.foam.directory, _EvapConfig.foam.folder, _rank);
    if (_macroSolver != NULL)
      std::cout << "Macro solver not null on rank: " << _rank << std::endl; // TODO: remove debug

    _multiMDService = new tarch::utils::MultiMDService<3>(_MDSolverConfig.getMPIConfiguration().getNumberOfProcesses(), _EvapConfig.totalNumberMDSimulations);
    _localMDInstances = _multiMDService->getLocalNumberOfMDSimulations();
    for (unsigned int i = 0; i < _localMDInstances; i++) {
      _MDSolver.push_back(coupling::interface::SimulationAndInterfaceFactory::getInstance().getMDSimulation(_MDSolverConfig, _mamicoConfig
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
                                                                                                            ,
                                                                                                            _multiMDService->getLocalCommunicator()
#endif
                                                                                                                ));
      if (_MDSolver[i] == NULL) {
        std::cout << "ERROR EvaporationTest: _MDSolver[" << i << "]==NULL!" << std::endl;
        exit(EXIT_FAILURE);
      }
    }
    if (_rank == 0) {
      gettimeofday(&_tv.start, NULL);
    }
    for (unsigned int i = 0; i < _localMDInstances; i++) {
      _MDSolver[i]->init(*_multiMDService, _multiMDService->getGlobalNumberOfLocalMDSimulation(i));
    }
    // allocate coupling interfaces
    for (unsigned int i = 0; i < _localMDInstances; i++) {
      _mdSolverInterface.push_back(
          coupling::interface::SimulationAndInterfaceFactory::getInstance().getMDSolverInterface(_MDSolverConfig, _mamicoConfig, _MDSolver[i]));
      if (_mdSolverInterface[i] == NULL) {
        std::cout << "ERROR EvaporationTest: mdSolverInterface[" << i << "] == NULL!" << std::endl;
        exit(EXIT_FAILURE);
      }
    }
    coupling::interface::MacroscopicSolverInterface<3>* macroSolverInterface = getMacroSolverInterface();
    _multiMDCellService = new coupling::services::MultiMDCellService<MY_LINKEDCELL, 3>(_mdSolverInterface, macroSolverInterface, _MDSolverConfig,
                                                                                       _mamicoConfig, "evaporation.xml", *_multiMDService);

    // init indexing
    coupling::indexing::IndexingService<3>::getInstance().init(_MDSolverConfig, _mamicoConfig, macroSolverInterface, (unsigned int)_rank);
    // set couette solver interface in MamicoInterfaceProvider
    coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL, 3>::getInstance().setMacroscopicSolverInterface(macroSolverInterface);
    for (unsigned int i = 0; i < _localMDInstances; i++) {
      _MDSolver[i]->setMacroscopicCellService(&(_multiMDCellService->getMacroscopicCellService(i)));
      // compute and store temperature in macroscopic cells (temp=1.1 everywhere)
      _multiMDCellService->getMacroscopicCellService(i).computeAndStoreTemperature(_EvapConfig.temp);
    }
    // allocate buffers for CFDToMD/MDToCFD operations
    allocateCFDToMDBuffer(_multiMDCellService->getMacroscopicCellService(0).getIndexConversion(),*macroSolverInterface);
    allocateMDToCFDBuffer(_multiMDCellService->getMacroscopicCellService(0).getIndexConversion(),*macroSolverInterface);
    // finish time measurement for initialisation
    if (_rank == 0) {
      gettimeofday(&_tv.end, NULL);
      double runtime = (_tv.end.tv_sec - _tv.start.tv_sec) * 1000000 + (_tv.end.tv_usec - _tv.start.tv_usec);
      std::cout << "Initialization: " << (int)(runtime / 1000) << "ms" << std::endl;
    }
    if (_rank == 0) {
      gettimeofday(&_tv.start_total, NULL);
    }
    std::cout << "Finish EvaporationTest::initSolvers() " << std::endl;
  }

  /** @brief advances the continuum solver and collects data to send to md (fillCFDToMDBuffer())
   *  @param cycle the number of the current coupling time step */
  void advanceMacro(int cycle) {
    if (_macroSolver != NULL) {
      if (_rank == 0) {
        gettimeofday(&_tv.start, NULL);
      }
      _macroSolver->advance();
      if (_rank == 0) {
        gettimeofday(&_tv.end, NULL);
        _tv.macro += (_tv.end.tv_sec - _tv.start.tv_sec) * 1000000 + (_tv.end.tv_usec - _tv.start.tv_usec);
        // std::cout << "Finish _macroSolver->advance " << std::endl;
      }
      fillCFDToMDBuffer();
    }
    _multiMDCellService->sendFromMacro2MD(_buf.CFDToMDBuffer, _buf.globalCellIndices4CFDToMDBuffer);
    // std::cout << "Finish _multiMDCellService->sendFromMacro2MD " << std::endl;
  }

  /** @brief advances the md solver for one coupling time step and collect the data for the coupling
   *  @param cycle the number of the current coupling time step  */
  void advanceMicro(int cycle) {
    if (_rank == 0) {
      gettimeofday(&_tv.start, NULL);
    }
    // run MD instances
    for (unsigned int i = 0; i < _localMDInstances; i++) {
      // set macroscopic cell service and interfaces in MamicoInterfaceProvider
      coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL, 3>::getInstance().setMacroscopicCellService(
          &(_multiMDCellService->getMacroscopicCellService(i)));
      coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL, 3>::getInstance().setMDSolverInterface(_mdSolverInterface[i]);
      _MDSolver[i]->simulateTimestepsWithVacuum(_MDSolverConfig.getSimulationConfiguration().getNumberOfTimesteps(), _mdStepCounter);
      // plot macroscopic time step info in multi md service
      _multiMDCellService->getMacroscopicCellService(i).plotEveryMacroscopicTimestep(cycle);
    }
    _mdStepCounter += _MDSolverConfig.getSimulationConfiguration().getNumberOfTimesteps();

    if (_rank == 0) {
      gettimeofday(&_tv.end, NULL);
      _tv.micro += (_tv.end.tv_sec - _tv.start.tv_sec) * 1000000 + (_tv.end.tv_usec - _tv.start.tv_usec);
    }
    // send back data from MD instances and merge it
    //_multiMDCellService->sendFromMD2Macro(_buf.MDToCFDBuffer,_buf.globalCellIndices4MDToCFDBuffer);
  }

  /** @brief finalize the time measurement, and cleans up at the end of the simulation */
  void shutdown() {
    // finish time measurement for coupled simulation
    if (_rank == 0) {
      gettimeofday(&_tv.end, NULL);
      double time_total = (_tv.end.tv_sec - _tv.start_total.tv_sec) * 1000000 + (_tv.end.tv_usec - _tv.start_total.tv_usec);
      std::cout << "Finished all coupling cycles after " << time_total / 1000000 << " s" << std::endl;
      std::cout << "Time percentages Micro, Macro, Other: " << std::endl;
      std::cout << _tv.micro / time_total * 100 << ", " << _tv.macro / time_total * 100 << ",  " << (1 - (_tv.micro + _tv.macro) / time_total) * 100
                << std::endl;
    }
    // free buffers/arrays
    deleteBuffer(_buf.CFDToMDBuffer);
    if (_buf.globalCellIndices4CFDToMDBuffer != NULL) {
      delete[] _buf.globalCellIndices4CFDToMDBuffer;
      _buf.globalCellIndices4CFDToMDBuffer = NULL;
    }
    deleteBuffer(_buf.MDToCFDBuffer);
    if (_buf.globalCellIndices4MDToCFDBuffer != NULL) {
      delete[] _buf.globalCellIndices4MDToCFDBuffer;
      _buf.globalCellIndices4MDToCFDBuffer = NULL;
    }

    // shutdown MD simulation
    for (unsigned int i = 0; i < _localMDInstances; i++) {
      // the shutdown operation may also delete the md solver interface; therefore, we update the MD solver interface in the vector _mdSolverInteface after the
      // shutdown is completed
      coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL, 3>::getInstance().setMDSolverInterface(_mdSolverInterface[i]);
      _MDSolver[i]->shutdown();
      delete _MDSolver[i];
      _MDSolver[i] = NULL;
      _mdSolverInterface[i] = coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL, 3>::getInstance().getMDSolverInterface();
    }
    _MDSolver.clear();
    for (unsigned int i = 0; i < _localMDInstances; i++) {
      if (_mdSolverInterface[i] != NULL) {
        delete _mdSolverInterface[i];
        _mdSolverInterface[i] = NULL;
      }
    }
    _mdSolverInterface.clear();

    coupling::interface::MacroscopicSolverInterface<3>* macroSolverInterface =
        coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL, 3>::getInstance().getMacroscopicSolverInterface();

    if (_multiMDService != NULL) {
      delete _multiMDService;
      _multiMDService = NULL;
    }
    if (macroSolverInterface != NULL) {
      delete macroSolverInterface;
      macroSolverInterface = NULL;
    }
    if (_macroSolver != NULL) {
      delete _macroSolver;
      _macroSolver = NULL;
    }
    if (_multiMDCellService != NULL) {
      delete _multiMDCellService;
      _multiMDCellService = NULL;
    }

    std::cout << "Finish EvaporationTest::shutdown() " << std::endl;
  }

  /** computes global number of macroscopic cells from configs. Required by macro solver interface before MacroscopicCellService is initialised! */
  tarch::la::Vector<3, unsigned int> getGlobalNumberMacroscopicCells() const {
    return tarch::la::Vector<3, unsigned int>{1,1,5};
  }

  /** This is only done on rank 0.
   *  @brief allocates the CFDToMD buffer (with values for all macroscopic cells).
   *  @param indexConversion instance of the indexConversion
   *  @param macroSolverInterface interface for the continuum solver */
  void allocateCFDToMDBuffer(const coupling::IndexConversion<3>& indexConversion, coupling::interface::MacroscopicSolverInterface<3>& macroSolverInterface) {
    const unsigned int num = 3; // just 3 cells shall be sent in oneD
    // delete all potential entries of CFDToMDBuffer // todo: where should they come from, not really necessary
    deleteBuffer(_buf.CFDToMDBuffer);
    // TODO: Adapt for parallel setup, check where this is necessary
    // allocate array for cell indices
    unsigned int* indices = new unsigned int[num]{31, 40, 49}; // refers to third, fourth, fifth cell
    if (indices == NULL) {
      std::cout << "ERROR EvaporationTest::allocateCFDToMDBuffer(): indices==NULL!" << std::endl;
      exit(EXIT_FAILURE);
    }
    // allocate CFDToMDBuffer and initialise all entries, incl. indices
    for (unsigned int i = 0; i < num; i++) {
      _buf.CFDToMDBuffer.push_back(new coupling::datastructures::MacroscopicCell<3>());
      if (_buf.CFDToMDBuffer[_buf.CFDToMDBuffer.size() - 1] == NULL) {
        std::cout << "ERROR EvaporationTest::allocateCFDToMDBuffer: CFDToMDBuffer[" << _buf.CFDToMDBuffer.size() - 1 << "]==NULL!" << std::endl;
        exit(EXIT_FAILURE);
      }
    }
    _buf.globalCellIndices4CFDToMDBuffer = indices;
  }

  /** @brief allocates the MDToCFD-buffer.
   *  Exists only on rank 0. On all other ranks a NULL ptr is returned */
  void allocateMDToCFDBuffer(const coupling::IndexConversion<3>& indexConversion, coupling::interface::MacroscopicSolverInterface<3>& macroSolverInterface) {
    // determine global number of cells
    const unsigned int num = 2;
    // delete all potential entries of MDToCFDBuffer
    deleteBuffer(_buf.MDToCFDBuffer);
    // allocate array for cell indices
    unsigned int* indices = new unsigned int[num]{58,67};
    if (indices == NULL) {
      std::cout << "ERROR EvaporationTest::allocateMDToCFDBuffer(): indices==NULL!" << std::endl;
      exit(EXIT_FAILURE);
    }
    // allocate MDToCFDBuffer and initialise all entries, incl. indices
    for (unsigned int i = 0; i < num; i++) {
      _buf.MDToCFDBuffer.push_back(new coupling::datastructures::MacroscopicCell<3>());
      if (_buf.MDToCFDBuffer[_buf.MDToCFDBuffer.size() - 1] == NULL) {
        std::cout << "ERROR EvaporationTest::allocateMDToCFDBuffer: MDToCFDBuffer[" << _buf.MDToCFDBuffer.size() - 1 << "]==NULL!" << std::endl;
        exit(EXIT_FAILURE);
      }
    }
    _buf.globalCellIndices4MDToCFDBuffer = indices;
  }

  /** @brief deletes the data in the buffer for the macro to md transfer
   *  @param buffer the buffer to be cleaned */
  void deleteBuffer(std::vector<coupling::datastructures::MacroscopicCell<3>*>& buffer) const {
    // delete all potential entries of buffer
    for (unsigned int i = 0; i < buffer.size(); i++) {
      if (buffer[i] != NULL) {
        delete buffer[i];
        buffer[i] = NULL;
      }
    }
    buffer.clear();
  }

  /** @brief fills CFDToMD buffer with data from macro/continuum solver */
  void fillCFDToMDBuffer() const {
    const coupling::IndexConversion<3>& indexConversion{_multiMDCellService->getIndexConversion()};
    const unsigned int size = _buf.CFDToMDBuffer.size();
    const tarch::la::Vector<3, double> domainOffset(indexConversion.getGlobalMDDomainOffset());
    const tarch::la::Vector<3, double> macroscopicCellSize(indexConversion.getMacroscopicCellSize());

    for (unsigned int i = 0; i < size; i++) {
      // get global cell index vector
      const tarch::la::Vector<3, unsigned int> globalIndex(indexConversion.getGlobalVectorCellIndex(_buf.globalCellIndices4CFDToMDBuffer[i]));
      // determine cell midpoint
      tarch::la::Vector<3, double> cellMidPoint(domainOffset - 0.5 * macroscopicCellSize);
      for (unsigned int d = 0; d < 3; d++) {
        cellMidPoint[d] = cellMidPoint[d] + ((double)globalIndex[d]) * macroscopicCellSize[d];
      }

      double mass = macroscopicCellSize[0] * macroscopicCellSize[1] * macroscopicCellSize[2];
      mass *= _macroSolver->getDensity(cellMidPoint);

      // compute momentum
      tarch::la::Vector<3, double> momentum(mass * _macroSolver->getVelocity(cellMidPoint));
      _buf.CFDToMDBuffer[i]->setMicroscopicMass(mass);
      _buf.CFDToMDBuffer[i]->setMicroscopicMomentum(momentum);
    }
  }

  /** @returns the interface for the macro/continuum solve  */
  coupling::interface::MacroscopicSolverInterface<3>* getMacroSolverInterface() {
    const tarch::la::Vector<3, double> mdOffset{_MDSolverConfig.getDomainConfiguration().getGlobalDomainOffset()};
    const tarch::la::Vector<3, double> mamicoMeshsize{_mamicoConfig.getMacroscopicCellConfiguration().getMacroscopicCellSize()};
    const tarch::la::Vector<3, unsigned int> globalNumberMacroscopicCells{getGlobalNumberMacroscopicCells()};
    const unsigned int outerRegion{_mamicoConfig.getMomentumInsertionConfiguration().getInnerOverlap()};
    coupling::interface::MacroscopicSolverInterface<3>* interface = NULL;
    // compute number of cells of MD offset; detect any mismatches!
    tarch::la::Vector<3, unsigned int> offsetMDDomain(0);
    for (unsigned int d = 0; d < 3; d++) {
      if (mdOffset[d] < 0.0) {
        std::cout << "ERROR EvaporationTest::getMacroSolverInterface(...): mdOffset[" << d << "]<0.0!" << std::endl;
        exit(EXIT_FAILURE);
      }
      offsetMDDomain[d] = floor(mdOffset[d] / mamicoMeshsize[d] + 0.5);
      if (fabs((offsetMDDomain[d] * mamicoMeshsize[d] - mdOffset[d]) / mamicoMeshsize[d]) > 1.0e-8) {
        std::cout << "ERROR CouetteTest::getCouetteSolverInterface: MD offset and mesh size mismatch!" << std::endl;
        exit(EXIT_FAILURE);
      }
    }
    interface = new coupling::solvers::CouetteSolverInterface<3>(globalNumberMacroscopicCells, outerRegion);
    if (interface == NULL) {
      std::cout << "ERROR EvaporationTest::getMacroSolverInterface(...), rank=" << _rank << ": interface==NULL!" << std::endl;
      exit(EXIT_FAILURE);
    }
    return interface;
  }

  /** the buffers store macroscopic cells, so momentum and density will be transferred
   *  @brief holds the buffers for the data transfer */
  struct CouplingBuffer {
    /** @brief the buffer for data transfer from macro to md */
    std::vector<coupling::datastructures::MacroscopicCell<3>*> CFDToMDBuffer;
    /** @brief the global indices of the macroscopic cells in the CFDToMDBuffer */
    unsigned int* globalCellIndices4CFDToMDBuffer;
    /** @brief the buffer for data transfer from md to macro */
    std::vector<coupling::datastructures::MacroscopicCell<3>*> MDToCFDBuffer;
    /** @brief the global indices of the macroscopic cells in the MDToCFDBuffer*/
    unsigned int* globalCellIndices4MDToCFDBuffer;
  };

  /** @brief holds all the variables for the time measurement of a simulation
   *  @todo Piet*/
  struct TimingValues {
    /** @brief */
    timeval start_total;
    /** @brief */
    timeval start;
    /** @brief */
    timeval end;
    /** @brief */
    double micro{0};
    /** @brief */
    double macro{0};
  };


  /** @brief the rank of the current MPI process */
  int _rank{0};
  /** @brief the config data and information for MDSolver */
  simplemd::configurations::MolecularDynamicsConfiguration _MDSolverConfig;
  /** @brief the config data and information for MaMiCo*/
  coupling::configurations::MaMiCoConfiguration<3> _mamicoConfig;
  /** @brief the EvaporationConfig for the current setup */
  coupling::configurations::EvaporationConfig _EvapConfig;
  /** @brief the counter for the time steps, which are done within the md */
  unsigned int _mdStepCounter{0};
  /** @brief the macro/continuum solver */
  coupling::solvers::RhoCentralInterface4Evaporation* _macroSolver;
  /** @todo piet */
  tarch::utils::MultiMDService<3>* _multiMDService;
  /** @todo piet */
  coupling::services::MultiMDCellService<MY_LINKEDCELL, 3>* _multiMDCellService;
  /** @brief the current buffer for the coupling */
  CouplingBuffer _buf;
  /** @todo piet */
  unsigned int _localMDInstances;
  /** @brief the interface to the md solver */
  std::vector<coupling::interface::MDSolverInterface<MY_LINKEDCELL, 3>*> _mdSolverInterface;
  /** @brief the md solver */
  std::vector<coupling::interface::MDSimulation*> _MDSolver;
  TimingValues _tv;
};

#endif // _COUPLING_TESTS_EVAPORATIONTEST_H_
