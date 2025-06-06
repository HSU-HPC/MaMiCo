// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _COUPLING_SCENARIO_COUETTESCENARIO_H_
#define _COUPLING_SCENARIO_COUETTESCENARIO_H_

#include "coupling/ErrorEstimation.h"
#include "coupling/InstanceHandling.h"
#include "coupling/MultiMDMediator.h"
#include "coupling/services/ParallelTimeIntegrationService.h"
#include "coupling/solvers/CouetteSolver.h"
#include "coupling/solvers/LBCouetteSolver.h"
#include "simplemd/configurations/MolecularDynamicsConfiguration.h"
#include "tarch/configuration/ParseConfiguration.h"
#include "tarch/utils/MultiMDService.h"
#include "tarch/utils/RandomNumberService.h"
#include "tarch/utils/Utils.h"
#if (BUILD_WITH_OPENFOAM)
#include "coupling/solvers/IcoFoamInterface.h"
#endif
#include "coupling/configurations/CouetteConfiguration.h"
#include "coupling/configurations/MaMiCoConfiguration.h"
#include "coupling/indexing/IndexingService.h"
#include "coupling/interface/MDSimulationFactory.h"
#include "coupling/interface/impl/SimpleMD/SimpleMDSolverInterface.h"
#include "coupling/services/MultiMDCellService.h"
#include "coupling/solvers/CouetteSolverInterface.h"
#include "coupling/solvers/FDCouetteSolver.h"
#include "coupling/solvers/LBCouetteSolverInterface.h"

#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
#include <mpi.h>
#endif
#include <chrono>
#include <math.h>
#include <random>
#include <sys/time.h>

#if defined(LS1_MARDYN)
#include "coupling/interface/impl/ls1/LS1MDSolverInterface.h"
#include "coupling/interface/impl/ls1/LS1StaticCommData.h"
#include "utils/Logger.h"
using Log::global_log;
#endif

// This is ignored if you dont use synthetic MD. For further instructions cf.
// SYNTHETIC part of initSolvers().
#define SYNTHETICMD_SEQUENCE "SYNTHETIC-MD"

/** Features, depending on couette.xml config file:
 * -> one-way or two-way coupling
 * -> statup or oscillating flow scenario
 * -> using analytical or Lattice Bolzmann Couette flow solver
 * -> using (multi-instance) SimpleMD or synthetic MD data (with Gaussian noise)
 * -> CSV output option for filtered data
 * -> ability to delay scenario startup or two-way coupling for filter
 * initialisation time
 * -> signal-to-noise ratio (SNR) computation (between filter output and
 * macroscopic solver)
 * -> runtime measurements and logging separately for coupled simulation
 * components
 * -> 'tws-loop' feature for testing many POD time-window-size parameters in a
 * single run
 * -> dynamic multi-instance MD, adapt number of MD instances according to defined error
 * @brief Versatile configurable Couette flow test for noise-filtered
 * multi-instance Nie coupling.
 * @author Piet Jarmatz
 */
class CouetteScenario : public Scenario {
public:
  /** @brief simple constructor */
  CouetteScenario() : Scenario("CouetteScenario"), _generator(std::chrono::system_clock::now().time_since_epoch().count()) {}
  /** @brief a dummy destructor */
  virtual ~CouetteScenario() {}

  /** triggers void init(), runOneCouplingCycle() and shutdown()
   *  @brief runs the simulation */
  void run() override {
    init();
    if (_cfg.twsLoop) {
      twsLoop();
      return;
    }
    _timeIntegrationService->run(_cfg.couplingCycles);
    shutdown();
  }

  /** triggers initMPI(), parseConfiguration(), and initSolvers()
   *  @brief initialises everthing necessary for the test */
  void init() override {
#if defined(LS1_MARDYN)
    Log::global_log = std::make_unique<Log::Logger>(Log::Error); // Log::Info
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
    global_log->set_mpi_output_root(0);
#endif
#endif
    getRootRank();
    parseConfigurations();
    initSolvers();
  }

  /** it advances the macro (advanceMacro()) and micro solver (advanceMicro),
   *  computes the signal to noise ratio (computeSNR()) and sends the data from
   * the macro to the micro solver (twoWayCoupling())
   *  @brief combines the functioniality necessary for a cycle of the coupled
   * simulation  */
  void runOneCouplingCycle(int cycle) override {
    advanceMacro(cycle);
    varyMD(cycle);
    advanceMicro(cycle);
    computeSNR(cycle);
    twoWayCoupling(cycle);
    if (_cfg.totalNumberMDSimulations < 0) // dynamic MD
      _multiMDCellService->finishCycle(cycle, *_instanceHandling);
    if (_isRootRank) {
      // Output status info only every 10 seconds
      gettimeofday(&_tv.end, NULL);
      int runtime_ms = (int)(((_tv.end.tv_sec - _tv.output.tv_sec) * 1000000 + (_tv.end.tv_usec - _tv.output.tv_usec)) / 1000);
      if (runtime_ms > 10000) {
        std::cout << "Finish coupling cycle " << cycle << std::endl;
        gettimeofday(&_tv.output, NULL);
      }
    }
  }

  coupling::solvers::AbstractCouetteSolver<3>* getSolver() override { return _couetteSolver; }

protected:
  /** @brief initialises all MPI variables  */
  void getRootRank() {
    int rank = 0;
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
    _isRootRank = (rank == 0);
  }

  /** @brief Executes the entire test several times for a range of time-window-size parameters */
  void twsLoop() {
    for (_tws = _cfg.twsLoopMin; _tws <= _cfg.twsLoopMax; _tws += _cfg.twsLoopStep) {
      init();
      for (int cycle = 0; cycle < _cfg.couplingCycles; cycle++)
        runOneCouplingCycle(cycle);
      shutdown();
    }
  }

  /** @brief reads the configuration from the xml file and calls
   * parseCouetteScenarioConfiguration() */
  void parseConfigurations() {
    std::string filename("couette.xml");

    tarch::configuration::ParseConfiguration::parseConfiguration<simplemd::configurations::MolecularDynamicsConfiguration>(filename, "molecular-dynamics",
                                                                                                                           _simpleMDConfig);
    if (!_simpleMDConfig.isValid()) {
      std::cout << "ERROR CouetteScenario: Invalid SimpleMD config!" << std::endl;
      exit(EXIT_FAILURE);
    }
    tarch::configuration::ParseConfiguration::parseConfiguration<coupling::configurations::MaMiCoConfiguration<3>>(filename, "mamico", _mamicoConfig);
    if (!_mamicoConfig.isValid()) {
      std::cout << "ERROR CouetteScenario: Invalid MaMiCo config!" << std::endl;
      exit(EXIT_FAILURE);
    }

    _cfg = coupling::configurations::CouetteConfig::parseCouetteConfiguration(filename);

    if (_cfg.miSolverType != coupling::configurations::CouetteConfig::MicroSolverType::LS1 &&
        _simpleMDConfig.getDomainDecompConfiguration().getDecompType() ==
            simplemd::configurations::DomainDecompConfiguration::DecompositionType::STATIC_IRREG_RECT_GRID) {
      std::cout << "ERROR Currently, only LS1 supports irregular rectilinear domain decomposition! Please change md type to ls1, or decomposition to default!"
                << std::endl;
      exit(EXIT_FAILURE);
    }

#if defined(LS1_MARDYN)
    auto offset = _simpleMDConfig.getDomainConfiguration().getGlobalDomainOffset();
    coupling::interface::LS1StaticCommData::getInstance().setConfigFilename("ls1config.xml");
    coupling::interface::LS1StaticCommData::getInstance().setBoxOffsetAtDim(0, offset[0]); // temporary till ls1 offset is natively supported
    coupling::interface::LS1StaticCommData::getInstance().setBoxOffsetAtDim(1, offset[1]);
    coupling::interface::LS1StaticCommData::getInstance().setBoxOffsetAtDim(2, offset[2]);
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
    auto subdomainWeights = _simpleMDConfig.getDomainDecompConfiguration().getSubdomainWeights();
    coupling::interface::LS1StaticCommData::getInstance().setSubdomainWeights(subdomainWeights);
#endif
#endif
  }

  /** @brief initialises the macro and micro solver according to the setup from
   * the xml file and pre-proccses them */
  void initSolvers() {
    _timeIntegrationService = std::make_unique<coupling::services::ParallelTimeIntegrationService<3>>(_mamicoConfig, this);
    _rank = _timeIntegrationService->getRank(); // returns the rank inside local time domain

    coupling::interface::MacroscopicSolverInterface<3>* couetteSolverInterface = getCouetteSolverInterface(
        _mamicoConfig.getCouplingCellConfiguration().getCouplingCellSize()[0], _simpleMDConfig.getDomainConfiguration().getGlobalDomainOffset(),
        _mamicoConfig.getCouplingCellConfiguration().getCouplingCellSize(), getGlobalNumberCouplingCells(_simpleMDConfig, _mamicoConfig),
        _mamicoConfig.getMomentumInsertionConfiguration().getInnerOverlap());

    // init indexing
    if (_simpleMDConfig.getDomainDecompConfiguration().getDecompType() ==
        simplemd::configurations::DomainDecompConfiguration::DecompositionType::STATIC_IRREG_RECT_GRID) {
      coupling::indexing::IndexingService<3>::getInstance().initWithMDSize(
          _simpleMDConfig.getDomainDecompConfiguration().getSubdomainWeights(), _simpleMDConfig.getDomainConfiguration().getGlobalDomainSize(),
          _simpleMDConfig.getDomainConfiguration().getGlobalDomainOffset(), _simpleMDConfig.getMPIConfiguration().getNumberOfProcesses(),
          _mamicoConfig.getCouplingCellConfiguration().getCouplingCellSize(), _mamicoConfig.getParallelTopologyConfiguration().getParallelTopologyType(),
          _mamicoConfig.getMomentumInsertionConfiguration().getInnerOverlap(), (unsigned int)_rank
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
          ,
          _timeIntegrationService->getPintComm()
#endif
      );
    } else {
      coupling::indexing::IndexingService<3>::getInstance().initWithMDSize(
          _simpleMDConfig.getDomainConfiguration().getGlobalDomainSize(), _simpleMDConfig.getDomainConfiguration().getGlobalDomainOffset(),
          _simpleMDConfig.getMPIConfiguration().getNumberOfProcesses(), _mamicoConfig.getCouplingCellConfiguration().getCouplingCellSize(),
          _mamicoConfig.getParallelTopologyConfiguration().getParallelTopologyType(), _mamicoConfig.getMomentumInsertionConfiguration().getInnerOverlap(),
          (unsigned int)_rank
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
          ,
          _timeIntegrationService->getPintComm()
#endif
      );
    }

    // for timing measurements
    _tv.micro = 0;
    _tv.macro = 0;
    _tv.filter = 0;

    // even if _cfg.miSolverType == SYNTHETIC then
    // multiMDService, _mdSimulations, _mdSolverInterface etc need to be initialized

    unsigned int totNumMD;
    if (_cfg.totalNumberMDSimulations > 0)
      totNumMD = _cfg.totalNumberMDSimulations;
    else // dynamic case, start with _cfg.lowerBoundNumberMDSimulations MD
      totNumMD = _cfg.lowerBoundNumberMDSimulations;

    _multiMDService = new tarch::utils::MultiMDService<3>(_simpleMDConfig.getMPIConfiguration().getNumberOfProcesses(), totNumMD
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
                                                          ,
                                                          _timeIntegrationService->getPintComm()
#endif
    );

    _instanceHandling = new coupling::InstanceHandling<MY_LINKEDCELL, 3>(_simpleMDConfig, _mamicoConfig, *_multiMDService);
    if (_instanceHandling == nullptr) {
      std::cout << "ERROR CouetteScenario::initSolvers() : _instanceHandling == NULL!" << std::endl;
      std::exit(EXIT_FAILURE);
    }

    _mdStepCounter = 0;
    if (_isRootRank) {
      gettimeofday(&_tv.start, NULL);
      gettimeofday(&_tv.output, NULL);
    }

    if (_cfg.miSolverType == coupling::configurations::CouetteConfig::SIMPLEMD || _cfg.miSolverType == coupling::configurations::CouetteConfig::LS1) {
      // equilibrate MD
      _instanceHandling->switchOffCoupling();
      _instanceHandling->equilibrate(_cfg.equSteps, _mdStepCounter);
      _instanceHandling->switchOnCoupling();
      _mdStepCounter += _cfg.equSteps;
    }

    // allocate coupling interfaces
    _instanceHandling->setMDSolverInterface();
    _mdSolverInterface = _instanceHandling->getMDSolverInterface();

    if (_cfg.twsLoop) {
      // initialise coupling cell service for multi-MD case and set single
      // cell services in each MD simulation
      _multiMDCellService = new coupling::services::MultiMDCellService<MY_LINKEDCELL, 3>(_mdSolverInterface, couetteSolverInterface, _simpleMDConfig,
                                                                                         _mamicoConfig, "couette.xml", *_multiMDService, _tws);
    } else {
      // initialise coupling cell service for multi-MD case and set single
      // cell services in each MD simulation
      _multiMDCellService = new coupling::services::MultiMDCellService<MY_LINKEDCELL, 3>(_mdSolverInterface, couetteSolverInterface, _simpleMDConfig,
                                                                                         _mamicoConfig, "couette.xml", *_multiMDService);
    }

    // init filtering for all md instances
    _multiMDCellService->constructFilterPipelines();

    _multiMDMediator = new coupling::MultiMDMediator<MY_LINKEDCELL, 3>(*_multiMDCellService, *_instanceHandling, *_multiMDService, couetteSolverInterface);

    // allocate solvers
    _couetteSolver = NULL;
    _couetteSolver =
        getCouetteSolver(_mamicoConfig.getCouplingCellConfiguration().getCouplingCellSize()[0],
                         _simpleMDConfig.getSimulationConfiguration().getDt() * _simpleMDConfig.getSimulationConfiguration().getNumberOfTimesteps());

    if (_cfg.miSolverType == coupling::configurations::CouetteConfig::SIMPLEMD || _cfg.miSolverType == coupling::configurations::CouetteConfig::LS1) {
      // set couette solver interface in MamicoInterfaceProvider
      coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL, 3>::getInstance().setMacroscopicSolverInterface(couetteSolverInterface);

      _instanceHandling->setCouplingCellServices(*_multiMDCellService);
      // compute and store temperature in coupling cells (temp=1.1
      // everywhere)
      _multiMDCellService->computeAndStoreTemperature(_cfg.temp);
    }

    /*
     * A synthethic solver is modeled using a dynamically linked filter,
     * i.e. a lambda function producing artifical data in every filter step.
     *
     * This is how to properly instanciate and use a synthethic solver:
     * - Create a sequence named like the macro SYNTHETICMD_SEQUENCE in xml. Use
     * whatever input, but make sure the filter system's output is set to this
     * sequence.
     * - Set filtered-values = "macro-mass macro-momentum" for that sequence.
     * - Use that sequence as input for all sequences that want (unfiltered) MD
     * input.
     *
     * Note that this synthethic solver is designed to be used on sequential
     * mode only and with only one MD instance.
     *
     * TODO
     * - major bug when there is ONLY a FFF in a sequence
     * - reduce capture: most variables in lambda can be defined beforehand as
     * they are const (e.g. everything coming from cfg)
     * - totalNumberMDSimulations > 1 is theoretically possible with this
     * redesign. test it and remove the restriction for it to be 1
     */

    else if (_cfg.miSolverType == coupling::configurations::CouetteConfig::SYNTHETIC) {
      // Synthetic MD runs sequentially only, as described above.
      if (_cfg.macro2Md || _cfg.totalNumberMDSimulations != 1 || _cfg.lbNumberProcesses[0] != 1 || _cfg.lbNumberProcesses[1] != 1 ||
          _cfg.lbNumberProcesses[2] != 1) {
        throw std::runtime_error("ERROR: Syntethic MD is only available in sequential mode!");
      }

      // set couette solver interface in MamicoInterfaceProvider
      coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL, 3>::getInstance().setMacroscopicSolverInterface(couetteSolverInterface);

      coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL, 3>::getInstance().setCouplingCellService(&(_multiMDCellService->getCouplingCellService(0)));
      coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL, 3>::getInstance().setMDSolverInterface(_mdSolverInterface[0]);

      // Create new FilterFromFunction instance and insert it into Filtering
      // System.
      try {
        _multiMDCellService->getCouplingCellService(0)
            .getFilterPipeline()
            ->getSequence(SYNTHETICMD_SEQUENCE)
            ->addFilter(
                new std::function<std::vector<double>(std::vector<double>, std::vector<std::array<unsigned int, 3>>)>{
                    // applyScalar
                    [this](std::vector<double> inputScalars,                    // not actually used as input:
                                                                                // matching MCS's
                                                                                // addFilterToSequence(...)
                                                                                // signature
                           std::vector<std::array<unsigned int, 3>> cellIndices // not used either
                    ) {
                      if (_isRootRank) {
                        gettimeofday(&_tv.start, NULL);
                      }

                      // std::cout << "Entering synthetic MD scalar..." << std::endl;

                      const unsigned int size = inputScalars.size();
                      auto dx = coupling::indexing::IndexingService<3>::getInstance().getCouplingCellSize();
                      const tarch::la::Vector<3, double> couplingCellSize(dx);
                      const double mass = (_cfg.density) * couplingCellSize[0] * couplingCellSize[1] * couplingCellSize[2];

                      std::vector<double> syntheticMasses;
                      for (unsigned int i = 0; i < size; i++) {
                        syntheticMasses.push_back(mass);
                      }
                      // std::cout << "Generated masses!" << std::endl;

                      if (_isRootRank) {
                        gettimeofday(&_tv.end, NULL);
                        _tv.micro += (_tv.end.tv_sec - _tv.start.tv_sec) * 1000000 + (_tv.end.tv_usec - _tv.start.tv_usec);
                      }

                      return syntheticMasses;
                    }},
                new std::function<std::vector<std::array<double, 3>>(std::vector<std::array<double, 3>>, std::vector<std::array<unsigned int, 3>>)>{
                    // applyVector
                    [this](std::vector<std::array<double, 3>> inputVectors,     // same as above
                           std::vector<std::array<unsigned int, 3>> cellIndices // same as above
                    ) {
                      if (_isRootRank) {
                        gettimeofday(&_tv.start, NULL);
                      }

                      // std::cout << "Entering synthetic MD vector." << std::endl;

                      const unsigned int size = inputVectors.size();
                      auto dx = coupling::indexing::IndexingService<3>::getInstance().getCouplingCellSize();
                      const tarch::la::Vector<3, double> couplingCellSize(dx);
                      const double mass = (_cfg.density) * couplingCellSize[0] * couplingCellSize[1] * couplingCellSize[2];

                      using coupling::indexing::IndexTrait;
                      using coupling::indexing::CellIndex;

                      const tarch::la::Vector<3, double> md2MacroDomainOffset = {
                          CellIndex<3, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::lowerBoundary.get()[0] * couplingCellSize[0],
                          CellIndex<3, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::lowerBoundary.get()[1] * couplingCellSize[1],
                          CellIndex<3, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::lowerBoundary.get()[2] * couplingCellSize[2],
                      };

                      std::normal_distribution<double> distribution(0.0, _cfg.noiseSigma);
                      std::vector<std::array<double, 3>> syntheticMomenta;
                      for (unsigned int i = 0; i < size; i++) {
                        // determine cell midpoint
                        CellIndex<3, IndexTrait::vector> globalIndex(
                            CellIndex<3, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>{i}); // construct global CellIndex from buffer
                                                                                                            // and convert it to vector
                        tarch::la::Vector<3, double> cellMidPoint(md2MacroDomainOffset - 0.5 * couplingCellSize);
                        for (unsigned int d = 0; d < 3; d++) {
                          cellMidPoint[d] = cellMidPoint[d] + ((double)globalIndex.get()[d]) * couplingCellSize[d];
                        }

                        // compute momentum
                        const tarch::la::Vector<3, double> noise(distribution(_generator), distribution(_generator), distribution(_generator));
                        const tarch::la::Vector<3, double> momentum(mass * ((*_couetteSolver).getVelocity(cellMidPoint) + noise));

                        // conversion from tarch::la::Vector to std::array
                        syntheticMomenta.push_back({momentum[0], momentum[1], momentum[2]});
                      }
                      // std::cout << "Generated momenta!" << std::endl;

                      if (_isRootRank) {
                        gettimeofday(&_tv.end, NULL);
                        _tv.micro += (_tv.end.tv_sec - _tv.start.tv_sec) * 1000000 + (_tv.end.tv_usec - _tv.start.tv_usec);
                      }

                      return syntheticMomenta;
                    }},
                0 // filterIndex
            );
      } catch (std::runtime_error& e) {
        auto expectedError = std::string("ERROR: Could not find Filter Sequence named ").append(SYNTHETICMD_SEQUENCE);
        if (expectedError.compare(e.what()) == 0) {
          std::cout << "ERROR: Synthetic MD solver selected without providing "
                       "filter sequence '"
                    << SYNTHETICMD_SEQUENCE << "' in config." << std::endl;
          exit(EXIT_FAILURE);
        } else
          throw e;
      }
    }

    // allocate buffers for send/recv operations
    allocateMacro2mdBuffer(*couetteSolverInterface);
    allocateMd2macroBuffer(*couetteSolverInterface);

    if (_cfg.initAdvanceCycles > 0 && _couetteSolver != NULL)
      _couetteSolver->advance(_cfg.initAdvanceCycles * _simpleMDConfig.getSimulationConfiguration().getDt() *
                              _simpleMDConfig.getSimulationConfiguration().getNumberOfTimesteps());

    // finish time measurement for initialisation
    if (_isRootRank) {
      gettimeofday(&_tv.end, NULL);
      double runtime = (_tv.end.tv_sec - _tv.start.tv_sec) * 1000000 + (_tv.end.tv_usec - _tv.start.tv_usec);
      std::cout << "Initialization: " << (int)(runtime / 1000) << "ms" << std::endl;
    }

    if (_cfg.computeSNR) {
      std::cout << "Output for every coupling cycle, for the cell 87 in recvBuffer:" << std::endl;
      std::cout << "cycle number (after filter-init-cycles), vel_x macroscopic "
                   "solver, vel_x filter output"
                << std::endl;
      _sum_signal = 0;
      _sum_noise = 0;
    }

    if (_isRootRank) {
      gettimeofday(&_tv.start_total, NULL);
    }
    std::cout << "Finish CouetteScenario::initSolvers() " << std::endl;
  }

  /** @brief advances the continuum solver and collects data to send to md
   * (fillSendBuffer())
   *  @param cycle the number of the current coupling time step */
  void advanceMacro(int cycle) {
    if (_couetteSolver != NULL) {
      if (_isRootRank) {
        gettimeofday(&_tv.start, NULL);
      }

      // run one time step for macroscopic couette solver
      if (_cfg.wallInitCycles > 0 && cycle == _cfg.wallInitCycles) {
        _couetteSolver->setWallVelocity(_cfg.wallVelocity);
        // When using Synthetic MD,
      }
      if (_cfg.wallOscillations != 0) {
        tarch::la::Vector<3, double> vel = cycle < _cfg.wallInitCycles ? _cfg.wallInitVelocity : _cfg.wallVelocity;
        vel = vel * cos(2 * M_PI * _cfg.wallOscillations * cycle / _cfg.couplingCycles);
        _couetteSolver->setWallVelocity(vel);
      }
      _couetteSolver->advance(_simpleMDConfig.getSimulationConfiguration().getDt() * _simpleMDConfig.getSimulationConfiguration().getNumberOfTimesteps());
      if (_isRootRank) {
        gettimeofday(&_tv.end, NULL);
        _tv.macro += (_tv.end.tv_sec - _tv.start.tv_sec) * 1000000 + (_tv.end.tv_usec - _tv.start.tv_usec);
        // std::cout << "Finish _couetteSolver->advance " << std::endl;
      }

      // extract data from couette solver and send them to MD (can take any
      // index-conversion object)
      fillSendBuffer(_cfg.density, *_couetteSolver, _couplingBuffer.macro2MDBuffer);
    }
    if (_cfg.macro2Md) {
#ifdef USE_COLLECTIVE_MPI
      _multiMDCellService->bcastFromMacro2MD(_couplingBuffer.macro2MDBuffer);
#else
      _multiMDCellService->sendFromMacro2MD(_couplingBuffer.macro2MDBuffer);
#endif
      // std::cout << "Finish _multiMDCellService->sendFromMacro2MD " <<
      // std::endl;
    }
  }

  void varyMD(int cycle) {
    // Non-dynamic case: Nothing to do.
    if (_cfg.totalNumberMDSimulations > 0)
      return;

    // modify number of instances only once every 10 cycles
    if (cycle % 10 != 0)
      return;

    int addMDInstances = 0;

    if (_rank == 0) {
      tarch::la::Vector<3, double> vel(0, 0, 0);
      double mass = 0;
      for (auto pair : _couplingBuffer.md2macroBuffer) {
        I01 idx;
        coupling::datastructures::CouplingCell<3>* couplingCell;
        std::tie(couplingCell, idx) = pair;
        vel += couplingCell->getMacroscopicMomentum();
        mass += couplingCell->getMacroscopicMass();
      }
      vel = vel / (double)_couplingBuffer.md2macroBuffer.size();
      mass /= _couplingBuffer.md2macroBuffer.size();

      double soundSpeed =
          (1 / std::sqrt(3)) * (_mamicoConfig.getCouplingCellConfiguration().getCouplingCellSize()[0] /
                                (_simpleMDConfig.getSimulationConfiguration().getDt() * _simpleMDConfig.getSimulationConfiguration().getNumberOfTimesteps()));
      double cellVolume = _mamicoConfig.getCouplingCellConfiguration().getCouplingCellSize()[0] *
                          _mamicoConfig.getCouplingCellConfiguration().getCouplingCellSize()[1] *
                          _mamicoConfig.getCouplingCellConfiguration().getCouplingCellSize()[2];

      coupling::error::ErrorEstimation errorControl(vel[0], _cfg.temp, mass, _simpleMDConfig.getMoleculeConfiguration().getMass(), soundSpeed,
                                                    _multiMDService->getTotalNumberOfMDSimulations(), cellVolume);

      double simtime = (double)cycle / _cfg.couplingCycles;
      double targetError = _cfg.absVelErrEnd * simtime + _cfg.absVelErrStart * (1 - simtime);
      errorControl.setAbsVelocityError(targetError);
      double NoMD = errorControl.getCorrectorNumberOfSamples(coupling::error::ErrorEstimation::Velocity, coupling::error::ErrorEstimation::Absolute);

      addMDInstances = NoMD - _multiMDService->getTotalNumberOfMDSimulations();
    }

#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
    MPI_Bcast(&addMDInstances, 1, MPI_INT, 0, _timeIntegrationService->getPintComm());
#endif
    if (addMDInstances < 0) {
      if ((int)_multiMDMediator->getNumberOfActiveMDSimulations() + addMDInstances < _cfg.lowerBoundNumberMDSimulations) {
        addMDInstances = (_cfg.lowerBoundNumberMDSimulations - (int)_multiMDMediator->getNumberOfActiveMDSimulations());
      }
    }
    if (addMDInstances < 0) {
      if (_rank == 0)
        std::cout << "Removing " << -addMDInstances << " of " << _multiMDService->getTotalNumberOfMDSimulations() << " MD simulations" << std::endl;
      _multiMDMediator->rmNMDSimulations(-addMDInstances);
    } else if (addMDInstances > 0) {
      if (_rank == 0)
        std::cout << "Adding " << addMDInstances << " to " << _multiMDService->getTotalNumberOfMDSimulations() << " MD simulations" << std::endl;
      _multiMDMediator->addNMDSimulations(addMDInstances);
    }
  }

  /** @brief advances the md solver for one coupling time step and collect the
   * data for the coupling
   *  @param cycle the number of the current coupling time step  */
  void advanceMicro(int cycle) {
    if (_isRootRank) {
      gettimeofday(&_tv.start, NULL);
    }
    if (_cfg.miSolverType == coupling::configurations::CouetteConfig::SIMPLEMD || _cfg.miSolverType == coupling::configurations::CouetteConfig::LS1) {
      // run MD instances
      _instanceHandling->simulateTimesteps(_simpleMDConfig.getSimulationConfiguration().getNumberOfTimesteps(), _mdStepCounter, *_multiMDCellService);
      // plot macroscopic time step info in multi md service
      _multiMDCellService->plotEveryMacroscopicTimestep(cycle);

      _mdStepCounter += _simpleMDConfig.getSimulationConfiguration().getNumberOfTimesteps();

      if (_isRootRank) {
        gettimeofday(&_tv.end, NULL);
        _tv.micro += (_tv.end.tv_sec - _tv.start.tv_sec) * 1000000 + (_tv.end.tv_usec - _tv.start.tv_usec);
      }

      // send back data from MD instances and merge it
      if (_cfg.md2Macro) {
#ifdef USE_COLLECTIVE_MPI
        _tv.filter += _multiMDCellService->reduceFromMD2Macro(_couplingBuffer.md2macroBuffer);
#else
        _tv.filter += _multiMDCellService->sendFromMD2Macro(_couplingBuffer.md2macroBuffer);
#endif
        // std::cout << "Finish _multiMDCellService->sendFromMD2Macro " <<
        // std::endl;
      }
    }

    if (_cfg.miSolverType == coupling::configurations::CouetteConfig::SYNTHETIC) {
      if (_cfg.md2Macro) {
        //_couplingBuffer does not get used here: Instead, the synthetic MD in the
        // SYNTHETICMD_SEQUENCE generates values. To prevent segfaults, it has
        // to be nonempty, though.
#ifdef USE_COLLECTIVE_MPI
        _tv.filter += _multiMDCellService->reduceFromMD2Macro(_couplingBuffer.md2macroBuffer);
#else
        _tv.filter += _multiMDCellService->sendFromMD2Macro(_couplingBuffer.md2macroBuffer);
#endif
      }
    }
  }

  /** @brief used to compute signal-to-noise ratio, stores values in _sum_noise and _sum_signal */
  void computeSNR(int cycle) {
    if (_cfg.computeSNR && cycle >= _cfg.filterInitCycles) {
      std::cout << cycle - _cfg.filterInitCycles << ", ";
      using namespace coupling::indexing;
      const tarch::la::Vector<3, double> dx(IndexingService<3>::getInstance().getCouplingCellSize());
      const double mass = _cfg.density * dx[0] * dx[1] * dx[2];
      unsigned int i = 0;
      for (auto pair : _couplingBuffer.md2macroBuffer) {
        /// todo@ use more cells
        if (i == 87) {
          I01 idx;
          coupling::datastructures::CouplingCell<3>* couplingCell;
          std::tie(couplingCell, idx) = pair;
          auto midPoint = idx.getCellMidPoint();
          double vx_macro = _couetteSolver->getVelocity(midPoint)[0];
          double vx_filter = (1 / mass * couplingCell->getMacroscopicMomentum())[0];
          _sum_noise += (vx_macro - vx_filter) * (vx_macro - vx_filter);
          _sum_signal += vx_macro * vx_macro;
          std::cout << vx_macro << ", " << vx_filter << std::endl;
        }
        i++;
      }
    }
  }

  /** @brief sets up the boundaries in the macro solver for the coupling and
   * applies the values from the md in the macro solver
   *  @param cycle current time step of the coupled simulation */
  void twoWayCoupling(int cycle) {
    using coupling::configurations::CouetteConfig;
    if (_cfg.twoWayCoupling) {
      if ((_cfg.maSolverType == CouetteConfig::COUETTE_LB || _cfg.maSolverType == CouetteConfig::COUETTE_FD) && cycle == _cfg.filterInitCycles) {
        static_cast<coupling::solvers::LBCouetteSolver*>(_couetteSolver)
            ->setMDBoundary(_simpleMDConfig.getDomainConfiguration().getGlobalDomainOffset(), _simpleMDConfig.getDomainConfiguration().getGlobalDomainSize(),
                            _mamicoConfig.getMomentumInsertionConfiguration().getInnerOverlap());
      }
#if (BUILD_WITH_OPENFOAM)
      else if ((_cfg.maSolverType == CouetteConfig::COUETTE_FOAM) && cycle == _cfg.filterInitCycles && _couetteSolver != NULL) {
        static_cast<coupling::solvers::IcoFoamInterface*>(_couetteSolver)->setupMDBoundary();
      }
#endif
      if ((_cfg.maSolverType == CouetteConfig::COUETTE_LB || _cfg.maSolverType == CouetteConfig::COUETTE_FD) && cycle >= _cfg.filterInitCycles) {
        static_cast<coupling::solvers::LBCouetteSolver*>(_couetteSolver)->setMDBoundaryValues(_couplingBuffer.md2macroBuffer);
      }
#if (BUILD_WITH_OPENFOAM)
      else if (_cfg.maSolverType == CouetteConfig::COUETTE_FOAM && cycle >= _cfg.filterInitCycles && _couetteSolver != NULL) {
        static_cast<coupling::solvers::IcoFoamInterface*>(_couetteSolver)->setMDBoundaryValues(_couplingBuffer.md2macroBuffer);
      }
#endif
    }
    // write data to csv-compatible file for evaluation
    write2CSV(_couplingBuffer.md2macroBuffer, cycle + 1);
  }

  /** @brief finalize the time measurement, and cleans up at the end of the
   * simulation */
  void shutdown() {
    if (_cfg.computeSNR) {
      std::cout << "SNR = " << 10 * log10(_sum_signal / _sum_noise) << std::endl;
    }

    // finish time measurement for coupled simulation
    if (_isRootRank) {
      gettimeofday(&_tv.end, NULL);
      double time_total = (_tv.end.tv_sec - _tv.start_total.tv_sec) * 1000000 + (_tv.end.tv_usec - _tv.start_total.tv_usec);
      std::cout << "Finished all coupling cycles after " << time_total / 1000000 << " s" << std::endl;
      if (_cfg.twsLoop)
        std::cout << "TWS = " << _tws << std::endl;
      std::cout << "Time percentages Micro, Macro, Filter, Other: " << std::endl;
      std::cout << _tv.micro / time_total * 100 << ", " << _tv.macro / time_total * 100 << ",  " << _tv.filter / time_total * 100 << ", "
                << (1 - (_tv.micro + _tv.macro + _tv.filter) / time_total) * 100 << std::endl;
    }

    if (_instanceHandling != nullptr) {
      delete _instanceHandling;
    }

    coupling::interface::MacroscopicSolverInterface<3>* couetteSolverInterface =
        coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL, 3>::getInstance().getMacroscopicSolverInterface();

    if (_multiMDService != NULL) {
      delete _multiMDService;
      _multiMDService = NULL;
    }
    if (couetteSolverInterface != NULL) {
      delete couetteSolverInterface;
      couetteSolverInterface = NULL;
    }
    if (_couetteSolver != NULL) {
      delete _couetteSolver;
      _couetteSolver = NULL;
    }
    if (_multiMDCellService != NULL) {
      delete _multiMDCellService;
      _multiMDCellService = NULL;
    }
    if (_multiMDMediator != nullptr) {
      delete _multiMDMediator;
      _multiMDMediator = nullptr;
    }

    std::cout << "Finish CouetteScenario::shutdown() " << std::endl;
  }

  /** computes global number of coupling cells from configs. Required by couette solver interface before CouplingCellService is initialised! */
  tarch::la::Vector<3, unsigned int> getGlobalNumberCouplingCells(const simplemd::configurations::MolecularDynamicsConfiguration& simpleMDConfig,
                                                                  const coupling::configurations::MaMiCoConfiguration<3>& mamicoConfig) const {
    tarch::la::Vector<3, double> domainSize(simpleMDConfig.getDomainConfiguration().getGlobalDomainSize());
    tarch::la::Vector<3, double> dx(mamicoConfig.getCouplingCellConfiguration().getCouplingCellSize());
    tarch::la::Vector<3, unsigned int> globalNumberCouplingCells(0);
    for (unsigned int d = 0; d < 3; d++) {
      int buf = floor(domainSize[d] / dx[d] + 0.5);
      globalNumberCouplingCells[d] = (unsigned int)buf;
    }
    return globalNumberCouplingCells;
  }

  /**
   *  @brief allocates the send buffer (with values for all coupling cells).
   *  @param msi macroscopic solver interface for the continuum solver */
  void allocateMacro2mdBuffer(coupling::interface::MacroscopicSolverInterface<3>& msi) {
    for (auto idx : I08()) {
      if (!I12::contains(idx)) {
        if (tarch::utils::contains(msi.getSourceRanks(idx), (unsigned int)_rank)) {
          coupling::datastructures::CouplingCell<3>* couplingCell = new coupling::datastructures::CouplingCell<3>();
          _couplingBuffer.macro2MDBuffer << std::make_pair(couplingCell, idx);
          if (couplingCell == nullptr)
            throw std::runtime_error(std::string("ERROR CouetteScenario::allocateMacro2mdBuffer: couplingCell==NULL!"));
        }
      }
    }
  }

  /** allocates the recv-buffer. This buffer contains all global inner coupling cells, but only on rank 0. On all other ranks, no cells are stored and a NULL
   * ptr is returned */
  void allocateMd2macroBuffer(coupling::interface::MacroscopicSolverInterface<3>& msi) {
    for (auto idx : I12()) {
      if (tarch::utils::contains(msi.getTargetRanks(idx), (unsigned int)_rank)) {
        coupling::datastructures::CouplingCell<3>* couplingCell = new coupling::datastructures::CouplingCell<3>();
        _couplingBuffer.md2macroBuffer << std::make_pair(couplingCell, idx);
        if (couplingCell == nullptr)
          throw std::runtime_error(std::string("ERROR CouetteScenario::allocateMacro2mdBuffer: couplingCell==NULL!"));
      }
    }
  }

  /** @brief write coupling cells that have been received from MD to csv file
   *  @param recvBuffer the buffer for the data, which comes from md
   *  @param recvIndices the indices for the macr cells in the buffer
   *  @param couplingCycle the current number of coupling cycle */
  template <class Container_T> void write2CSV(Container_T& md2macroBuffer, int couplingCycle) const {
    if (md2macroBuffer.size() == 0)
      return;
    if (_cfg.csvEveryTimestep < 1 || couplingCycle % _cfg.csvEveryTimestep > 0)
      return;
    // form file name and open file
    std::stringstream ss;
    ss << "CouetteAvgMultiMDCells_" << _timeIntegrationService->getPintDomain() << "_" << _rank << "_" << couplingCycle << ".csv";
    std::ofstream file(ss.str().c_str());
    if (!file.is_open()) {
      std::cout << "ERROR CouetteScenario::write2CSV(): Could not open file " << ss.str() << "!" << std::endl;
      exit(EXIT_FAILURE);
    }

    // loop over md2macro cells; read macroscopic mass+momentum buffers and
    // write cell index, mass and velocity to one line in the csv-file
    file << "I01_x,I01_y,I01_z,vel_x,vel_y,vel_z,mass" << std::endl;
    for (auto pair : md2macroBuffer) {
      I01 idx;
      coupling::datastructures::CouplingCell<3>* couplingCell;
      std::tie(couplingCell, idx) = pair;
      tarch::la::Vector<3, double> vel(couplingCell->getMacroscopicMomentum());
      if (couplingCell->getMacroscopicMass() != 0.0) {
        vel = (1.0 / couplingCell->getMacroscopicMass()) * vel;
      }
      file << idx << "," << vel[0] << "," << vel[1] << "," << vel[2] << "," << couplingCell->getMacroscopicMass() << std::endl;
    }

    // close file
    file.close();
  }

  /** @brief fills send buffer with data from macro/continuum solver
   *  @param density the general density of the fluid
   *  @param couetteSolver the continuum solver
   *  @param macro2MDBuffer the bufffer to send data from macro to micro
   *  @param globalCellIndices4SendBuffer the global linearized indices of the
   * coupling cells in the buffer  */
  void fillSendBuffer(const double density, const coupling::solvers::AbstractCouetteSolver<3>& couetteSolver,
                      coupling::datastructures::FlexibleCellContainer<3>& macro2MDBuffer) const {
    using coupling::configurations::CouetteConfig;
    using namespace coupling::indexing;
    const tarch::la::Vector<3, double> dx(IndexingService<3>::getInstance().getCouplingCellSize());
    double mass = density * dx[0] * dx[1] * dx[2];
    for (auto pair : macro2MDBuffer) {
      I01 idx;
      coupling::datastructures::CouplingCell<3>* couplingCell;
      std::tie(couplingCell, idx) = pair;
      auto midPoint = idx.getCellMidPoint();
      if (_cfg.maSolverType == CouetteConfig::COUETTE_LB || _cfg.maSolverType == CouetteConfig::COUETTE_FD)
        mass *= static_cast<const coupling::solvers::LBCouetteSolver*>(&couetteSolver)->getDensity(midPoint);
      // compute momentum
      tarch::la::Vector<3, double> momentum(mass * couetteSolver.getVelocity(midPoint));
      couplingCell->setMicroscopicMass(mass);
      couplingCell->setMicroscopicMomentum(momentum);
    }
  }

  /** @returns the correct macro/continuum solver for the given setup
   *  @param dx the grid size (equidistant mesh)
   *  @param dt the time step  */
  coupling::solvers::AbstractCouetteSolver<3>* getCouetteSolver(const double dx, const double dt) {
    using coupling::configurations::CouetteConfig;
    coupling::solvers::AbstractCouetteSolver<3>* solver = NULL;
    tarch::la::Vector<3, double> vel = _cfg.wallInitCycles > 0 ? _cfg.wallInitVelocity : _cfg.wallVelocity;
    // analytical solver: is only active on rank 0
    if (_cfg.maSolverType == CouetteConfig::COUETTE_ANALYTICAL) {
      if (_rank == 0 || _cfg.miSolverType == CouetteConfig::SYNTHETIC) {
        solver = new coupling::solvers::CouetteSolver<3>(_cfg.channelheight, vel[0], _cfg.kinVisc);
        if (solver == NULL) {
          std::cout << "ERROR CouetteScenario::getCouetteSolver(): Analytic solver==NULL!" << std::endl;
          exit(EXIT_FAILURE);
        }
      }
    }
#if (BUILD_WITH_OPENFOAM)
    else if (_cfg.maSolverType == CouetteConfig::COUETTE_FOAM) {
      solver = new coupling::solvers::IcoFoamInterface(_rank, _cfg.plotEveryTimestep, _cfg.channelheight, _cfg.foam.directory, _cfg.foam.folder,
                                                       _cfg.foam.boundariesWithMD, _cfg.wallVelocity);
      if (solver == NULL) {
        std::cout << "ERROR CouetteScenario::getCouetteSolver(): IcoFoamInterface solver==NULL!" << std::endl;
        exit(EXIT_FAILURE);
      }
    }
#endif
    // LB solver: active on lbNumberProcesses
    else if (_cfg.maSolverType == CouetteConfig::COUETTE_LB) {
      solver = new coupling::solvers::LBCouetteSolver(_cfg.channelheight, vel, _cfg.kinVisc, dx, dt, _cfg.plotEveryTimestep, "LBCouette",
                                                      _cfg.lbNumberProcesses, 1, this);
      if (solver == NULL) {
        std::cout << "ERROR CouetteScenario::getCouetteSolver(): LB solver==NULL!" << std::endl;
        exit(EXIT_FAILURE);
      }
    } else if (_cfg.maSolverType == CouetteConfig::COUETTE_FD) {
      solver = new coupling::solvers::FiniteDifferenceSolver(_cfg.channelheight, vel, _cfg.kinVisc, dx, dt, _cfg.plotEveryTimestep, "FDCouette",
                                                             _cfg.lbNumberProcesses, 1);
      if (solver == NULL) {
        std::cout << "ERROR CouetteScenario::getCouetteSolver(): FD solver==NULL!" << std::endl;
        exit(EXIT_FAILURE);
      }
    } else {
      std::cout << "ERROR CouetteScenario::getCouetteSolver(): Unknown solver type!" << std::endl;
      exit(EXIT_FAILURE);
    }
    return solver;
  }

  /** @returns the interface for the macro/continuum solver
   *  @param couetteSolver the macro/continuum solver
   *  @param mdOffset the offset of the md domain from (0.0.0)
   *  @param mamicoMeshsize
   *  @param globalNumberCouplingCells the total number coupling cells for
   * the whole domain
   *  @param outerRegion
   *  @todo piet, what is the mamicoMeshsize & the outer layer  */
  coupling::interface::MacroscopicSolverInterface<3>* getCouetteSolverInterface(const double dx, tarch::la::Vector<3, double> mdOffset,
                                                                                tarch::la::Vector<3, double> mamicoMeshsize,
                                                                                tarch::la::Vector<3, unsigned int> globalNumberCouplingCells,
                                                                                unsigned int outerRegion) {
    using coupling::configurations::CouetteConfig;
    coupling::interface::MacroscopicSolverInterface<3>* interface = NULL;
    if (_cfg.maSolverType == CouetteConfig::COUETTE_ANALYTICAL) {
      interface = new coupling::solvers::CouetteSolverInterface<3>(globalNumberCouplingCells, outerRegion);
    } else if (_cfg.maSolverType == CouetteConfig::COUETTE_LB) {
      // compute number of cells of MD offset; detect any mismatches!
      tarch::la::Vector<3, unsigned int> offsetMDDomain(0);
      for (unsigned int d = 0; d < 3; d++) {
        if (mdOffset[d] < 0.0) {
          std::cout << "ERROR CouetteScenario::getCouetteSolverInterface(...): mdOffset[" << d << "]<0.0!" << std::endl;
          exit(EXIT_FAILURE);
        }
        offsetMDDomain[d] = floor(mdOffset[d] / mamicoMeshsize[d] + 0.5);
        if (fabs((offsetMDDomain[d] * mamicoMeshsize[d] - mdOffset[d]) / mamicoMeshsize[d]) > 1.0e-8) {
          std::cout << "ERROR CouetteScenario::getCouetteSolverInterface: MD offset and mesh size mismatch!" << std::endl;
          exit(EXIT_FAILURE);
        }
      }
      tarch::la::Vector<3, unsigned int> cells_per_process{
          tarch::la::Vector<3, int>{coupling::solvers::NumericalSolver::getAvgDomainSize(_cfg.channelheight, dx, _cfg.lbNumberProcesses, 0),
                                    coupling::solvers::NumericalSolver::getAvgDomainSize(_cfg.channelheight, dx, _cfg.lbNumberProcesses, 1),
                                    coupling::solvers::NumericalSolver::getAvgDomainSize(_cfg.channelheight, dx, _cfg.lbNumberProcesses, 2)}};
      interface =
          new coupling::solvers::LBCouetteSolverInterface(cells_per_process, _cfg.lbNumberProcesses, offsetMDDomain, globalNumberCouplingCells, outerRegion);
    }
#if (BUILD_WITH_OPENFOAM)
    else if (_cfg.maSolverType == CouetteConfig::COUETTE_FOAM) {
      interface = new coupling::solvers::CouetteSolverInterface<3>(globalNumberCouplingCells, outerRegion);
    }
#endif
    else if (_cfg.maSolverType == CouetteConfig::COUETTE_FD) {
      // compute number of cells of MD offset; detect any mismatches!
      tarch::la::Vector<3, unsigned int> offsetMDDomain(0);
      for (unsigned int d = 0; d < 3; d++) {
        if (mdOffset[d] < 0.0) {
          std::cout << "ERROR CouetteScenario::getCouetteSolverInterface(...): mdOffset[" << d << "]<0.0!" << std::endl;
          exit(EXIT_FAILURE);
        }
        offsetMDDomain[d] = floor(mdOffset[d] / mamicoMeshsize[d] + 0.5);
        if (fabs((offsetMDDomain[d] * mamicoMeshsize[d] - mdOffset[d]) / mamicoMeshsize[d]) > 1.0e-8) {
          std::cout << "ERROR CouetteScenario::getCouetteSolverInterface: MD offset and mesh size mismatch!" << std::endl;
          exit(EXIT_FAILURE);
        }
      }
      tarch::la::Vector<3, unsigned int> cells_per_process{
          tarch::la::Vector<3, int>{coupling::solvers::NumericalSolver::getAvgDomainSize(_cfg.channelheight, dx, _cfg.lbNumberProcesses, 0),
                                    coupling::solvers::NumericalSolver::getAvgDomainSize(_cfg.channelheight, dx, _cfg.lbNumberProcesses, 1),
                                    coupling::solvers::NumericalSolver::getAvgDomainSize(_cfg.channelheight, dx, _cfg.lbNumberProcesses, 2)}};
      interface =
          new coupling::solvers::LBCouetteSolverInterface(cells_per_process, _cfg.lbNumberProcesses, offsetMDDomain, globalNumberCouplingCells, outerRegion);
    }

    if (interface == NULL) {
      std::cout << "ERROR CouetteScenario::getCouetteSolverInterface(...), rank=" << _rank << ": interface==NULL!" << std::endl;
      exit(EXIT_FAILURE);
    }
    return interface;
  }

  /** the buffers store coupling cells, so momentum and density will be
   * transferred
   *  @brief holds the buffers for the data transfer */
  struct CouplingBuffer {
    /** @brief the buffer for data transfer from macro to md */
    coupling::datastructures::FlexibleCellContainer<3> macro2MDBuffer;
    /** @brief the buffer for data transfer from md to macro */
    coupling::datastructures::FlexibleCellContainer<3> md2macroBuffer;
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
    timeval output;
    /** @brief */
    double micro;
    /** @brief */
    double macro;
    /** @brief */
    double filter;
  };

  /** @brief the rank of the current MPI process in the local time domain*/
  int _rank;
  /** @brief if this is the world global root process */
  bool _isRootRank;
  /** @todo Piet */
  int _tws;
  /** @brief the config data and information for SimpleMD */
  simplemd::configurations::MolecularDynamicsConfiguration _simpleMDConfig;
  /** @brief the config data and information for MaMiCo*/
  coupling::configurations::MaMiCoConfiguration<3> _mamicoConfig;
  /** @brief the CouetteConfig for the current setup */
  coupling::configurations::CouetteConfig _cfg;
  /** @brief the counter for the time steps, which are done within the md */
  unsigned int _mdStepCounter;
  /** @brief the current macro/continuum solver */
  coupling::solvers::AbstractCouetteSolver<3>* _couetteSolver;
  /** @todo piet */
  tarch::utils::MultiMDService<3>* _multiMDService;
  /** @todo piet */
  coupling::services::MultiMDCellService<MY_LINKEDCELL, 3>* _multiMDCellService;
  /** @brief the current buffer for the coupling */
  CouplingBuffer _couplingBuffer;
  /** @brief the interface to the md solver */
  std::vector<coupling::interface::MDSolverInterface<MY_LINKEDCELL, 3>*> _mdSolverInterface;
  coupling::InstanceHandling<MY_LINKEDCELL, 3>* _instanceHandling;
  /** @todo piet */
  std::default_random_engine _generator;
  /** @todo piet */
  double _sum_signal;
  /** @todo piet */
  double _sum_noise;
  /** @brief a instance of the timingValues */
  TimingValues _tv;
  coupling::MultiMDMediator<MY_LINKEDCELL, 3>* _multiMDMediator;
};

#endif // _COUPLING_SCENARIO_COUETTESCENARIO_H_
