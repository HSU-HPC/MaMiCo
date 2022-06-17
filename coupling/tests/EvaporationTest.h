// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _COUPLING_TESTS_EVAPORATIONTEST_H_
#define _COUPLING_TESTS_EVAPORATIONTEST_H_

#include "coupling/CouplingMDDefinitions.h"
#include "coupling/solvers/CouetteSolver.h"
#include "coupling/solvers/LBCouetteSolver.h"
#include "coupling/tests/Test.h"
#include "simplemd/configurations/MolecularDynamicsConfiguration.h"
#include "tarch/configuration/ParseConfiguration.h"
#include "tarch/utils/MultiMDService.h"
#if (BUILD_WITH_OPENFOAM)
#include "coupling/solvers/FoamClass.h"
#include "coupling/solvers/FoamSolverInterface.h"
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
    _macroSolver =
        getMacroSolver(_mamicoConfig.getMacroscopicCellConfiguration().getMacroscopicCellSize()[0],
                         _MDSolverConfig.getSimulationConfiguration().getDt() * _MDSolverConfig.getSimulationConfiguration().getNumberOfTimesteps());
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
    coupling::interface::MacroscopicSolverInterface<3>* macroSolverInterface = getMacroSolverInterface(
        _macroSolver, _MDSolverConfig.getDomainConfiguration().getGlobalDomainOffset(),
        _mamicoConfig.getMacroscopicCellConfiguration().getMacroscopicCellSize(), getGlobalNumberMacroscopicCells(_MDSolverConfig, _mamicoConfig),
        _mamicoConfig.getMomentumInsertionConfiguration().getInnerOverlap());

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
      _macroSolver->advance(_MDSolverConfig.getSimulationConfiguration().getDt() * _MDSolverConfig.getSimulationConfiguration().getNumberOfTimesteps());
      if (_rank == 0) {
        gettimeofday(&_tv.end, NULL);
        _tv.macro += (_tv.end.tv_sec - _tv.start.tv_sec) * 1000000 + (_tv.end.tv_usec - _tv.start.tv_usec);
        // std::cout << "Finish _macroSolver->advance " << std::endl;
      }
    }
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
  tarch::la::Vector<3, unsigned int> getGlobalNumberMacroscopicCells(const simplemd::configurations::MolecularDynamicsConfiguration& MDSolverConfig,
                                                                     const coupling::configurations::MaMiCoConfiguration<3>& mamicoConfig) const {
    tarch::la::Vector<3, double> domainSize(MDSolverConfig.getDomainConfiguration().getGlobalDomainSize());
    tarch::la::Vector<3, double> dx(mamicoConfig.getMacroscopicCellConfiguration().getMacroscopicCellSize());
    tarch::la::Vector<3, unsigned int> globalNumberMacroscopicCells(0);
    for (unsigned int d = 0; d < 3; d++) {
      int buf = floor(domainSize[d] / dx[d] + 0.5);
      globalNumberMacroscopicCells[d] = (unsigned int)buf;
    }
    return globalNumberMacroscopicCells;
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
      std::cout << "the global vector cell index is CFD2MD" << indexConversion.getGlobalVectorCellIndex(indices[i]) << std::endl;
      _buf.CFDToMDBuffer.push_back(new coupling::datastructures::MacroscopicCell<3>());
      if (_buf.CFDToMDBuffer[_buf.CFDToMDBuffer.size() - 1] == NULL) {
        std::cout << "ERROR EvaporationTest::allocateCFDToMDBuffer: CFDToMDBuffer[" << _buf.CFDToMDBuffer.size() - 1 << "]==NULL!" << std::endl;
        exit(EXIT_FAILURE);
      }
    }
    _buf.globalCellIndices4CFDToMDBuffer = indices;
  }

  /** allocates the MDToCFD-buffer. This buffer contains all global inner macroscopic cells, but only on rank 0. On all other ranks, no cells are stored and a NULL
   * ptr is returned */
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
      std::cout << "the global vector cell index is MD2CFD" << indexConversion.getGlobalVectorCellIndex(indices[i]) << std::endl;
      _buf.MDToCFDBuffer.push_back(new coupling::datastructures::MacroscopicCell<3>());
      if (_buf.MDToCFDBuffer[_buf.MDToCFDBuffer.size() - 1] == NULL) {
        std::cout << "ERROR EvaporationTest::allocateMDToCFDBuffer: MDToCFDBuffer[" << _buf.MDToCFDBuffer.size() - 1 << "]==NULL!" << std::endl;
        exit(EXIT_FAILURE);
      }
    }
    _buf.globalCellIndices4MDToCFDBuffer = indices;
  }

  /** @brief deletes the data in the buffer for the macro to md transfer
   *  @param CFDToMDBuffer the buffer to be cleaned */
  void deleteBuffer(std::vector<coupling::datastructures::MacroscopicCell<3>*>& CFDToMDBuffer) const {
    // delete all potential entries of CFDToMDBuffer
    for (unsigned int i = 0; i < CFDToMDBuffer.size(); i++) {
      if (CFDToMDBuffer[i] != NULL) {
        delete CFDToMDBuffer[i];
        CFDToMDBuffer[i] = NULL;
      }
    }
    CFDToMDBuffer.clear();
  }

  /** @brief fills CFDToMD buffer with data from macro/continuum solver
   *  @param density the general density of the fluid
   *  @param macroSolver the continuum solver
   *  @param indexConversion an instance of the indexConversion
   *  @param CFDToMDBuffer the bufffer to CFDToMD data from macro to micro
   *  @param globalCellIndices4CFDToMDBuffer the global linearized indices of the macroscopic cells in the buffer  */
  void fillCFDToMDBuffer(const double density, const coupling::solvers::AbstractCouetteSolver<3>& macroSolver,
                      const coupling::IndexConversion<3>& indexConversion, std::vector<coupling::datastructures::MacroscopicCell<3>*>& CFDToMDBuffer,
                      const unsigned int* const globalCellIndices4CFDToMDBuffer) const {
    const unsigned int size = CFDToMDBuffer.size();
    const tarch::la::Vector<3, double> domainOffset(indexConversion.getGlobalMDDomainOffset());
    const tarch::la::Vector<3, double> macroscopicCellSize(indexConversion.getMacroscopicCellSize());

    for (unsigned int i = 0; i < size; i++) {
      // get global cell index vector
      const tarch::la::Vector<3, unsigned int> globalIndex(indexConversion.getGlobalVectorCellIndex(globalCellIndices4CFDToMDBuffer[i]));
      // determine cell midpoint
      tarch::la::Vector<3, double> cellMidPoint(domainOffset - 0.5 * macroscopicCellSize);
      for (unsigned int d = 0; d < 3; d++) {
        cellMidPoint[d] = cellMidPoint[d] + ((double)globalIndex[d]) * macroscopicCellSize[d];
      }

      double mass = density * macroscopicCellSize[0] * macroscopicCellSize[1] * macroscopicCellSize[2];
      mass *= static_cast<const coupling::solvers::LBCouetteSolver*>(&macroSolver)->getDensity(cellMidPoint);

      // compute momentum
      tarch::la::Vector<3, double> momentum(mass * macroSolver.getVelocity(cellMidPoint));
      CFDToMDBuffer[i]->setMicroscopicMass(mass);
      CFDToMDBuffer[i]->setMicroscopicMomentum(momentum);
    }
  }

  /** @returns the correct marco/continuum solver for the given setup
   *  @param dx the grid size (equidistant mesh)
   *  @param dt the time step  */
  coupling::solvers::AbstractCouetteSolver<3>* getMacroSolver(const double dx, const double dt) {
    coupling::solvers::AbstractCouetteSolver<3>* solver = NULL;
    // LB solver: active on lbNumberProcesses
    solver = new coupling::solvers::LBCouetteSolver(_EvapConfig.channelheight, tarch::la::Vector<3, double>(0), 0, dx, dt, -1, "LBEvaporation",
                                                    tarch::la::Vector<3, unsigned int>(1));
    if (solver == NULL) {
      std::cout << "ERROR EvaporationTest::getEvaporationSolver(): LB solver==NULL!" << std::endl;
      exit(EXIT_FAILURE);
    }
    return solver;
  }

  /** @returns the interface for the macro/continuum solver
   *  @param macroSolver the macro/continuum solver
   *  @param mdOffset the offset of the md domain from (0.0.0)
   *  @param mamicoMeshsize
   *  @param globalNumberMacroscopicCells the total number macroscopic cells for the whole domain
   *  @param outerRegion
   *  @todo piet, what is the mamicoMeshsize & the outer layer  */
  coupling::interface::MacroscopicSolverInterface<3>* getMacroSolverInterface(coupling::solvers::AbstractCouetteSolver<3>* macroSolver,
                                                                                tarch::la::Vector<3, double> mdOffset,
                                                                                tarch::la::Vector<3, double> mamicoMeshsize,
                                                                                tarch::la::Vector<3, unsigned int> globalNumberMacroscopicCells,
                                                                                unsigned int outerRegion) {
    coupling::interface::MacroscopicSolverInterface<3>* interface = NULL;
    coupling::solvers::LBCouetteSolver* lbSolver = static_cast<coupling::solvers::LBCouetteSolver*>(macroSolver);
    if (lbSolver == NULL) {
      std::cout << "ERROR EvaporationTest::getMacroSolverInterface(...), rank=" << _rank << ": Could not convert abstract to LB solver!" << std::endl;
      exit(EXIT_FAILURE);
    }
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
  /** @brief the current macro/continuum solver */
  coupling::solvers::AbstractCouetteSolver<3>* _macroSolver{nullptr};
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
