#pragma once

#include "coupling/CouplingMDDefinitions.h"
#include "coupling/configurations/CouetteConfiguration.h"
#include "coupling/configurations/MaMiCoConfiguration.h"
#include "coupling/indexing/IndexingService.h"
#include "coupling/interface/MDSimulationFactory.h"
#include "coupling/interface/impl/SimpleMD/SimpleMDSolverInterface.h"
#include "coupling/services/ParallelTimeIntegrationService.h"
#include "coupling/scenario/Scenario.h"
#include "coupling/services/MultiMDCellService.h"
#include "coupling/solvers/CouetteSolver.h"
#include "coupling/solvers/CouetteSolverInterface.h"
#include "coupling/scenario/precice/PreciceAdapter.h"
#include "simplemd/configurations/MolecularDynamicsConfiguration.h"
#include "tarch/configuration/Configuration.h"
#include "tarch/configuration/ParseConfiguration.h"
#include "tarch/utils/MultiMDService.h"
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
#include <mpi.h>
#endif
#include <random>
#include <string>
#include <sys/time.h>

#if defined(LS1_MARDYN)
#include "coupling/interface/impl/ls1/LS1MDSolverInterface.h"
#include "coupling/interface/impl/ls1/LS1StaticCommData.h"
#include "utils/Logger.h"
using Log::global_log;
#endif

namespace coupling {
namespace scenario {
namespace precice {
  class CouetteScenario;
}
}
}


class coupling::scenario::precice::CouetteScenario : public Scenario {
public:
  /** @param configFilename couette.xml in the standard <couette-test> schema */
  CouetteScenario(const std::string& configFilename = "couette.xml") : Scenario("couette_precice"), _configFilename(configFilename) {}
  ~CouetteScenario() {
    if (_instanceHandling != nullptr) {
      delete _instanceHandling;
      _instanceHandling = nullptr;
    }
    if (_multiMDService != nullptr) {
      delete _multiMDService;
      _multiMDService = nullptr;
    }
    if (_preciceAdapter != nullptr) {
      delete _preciceAdapter;
      _preciceAdapter = nullptr;
    }
    if (_preciceInterface != nullptr) {
      delete _preciceInterface;
      _preciceInterface = nullptr;
    }
    if (_multiMDCellService != nullptr) {
      delete _multiMDCellService;
      _multiMDCellService = nullptr;
    }
  }

  /**
   * Run the Couette scenario for a certain amount of timesteps given in the configuration file.
   * PreCICE is steering the solving between the CFD solver and MaMiCo/MD solver.
   */
  void run () override {
    const unsigned int dim = 3;
    int rank = 0;
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
#if defined(LS1_MARDYN)
    Log::global_log = std::make_unique<Log::Logger>(Log::Error); //Log::Info
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
    global_log->set_mpi_output_root(0);
#endif
#endif

    // Same three-section couette.xml the in-process scenarios read; <macroscopic-solver type="precice">
    simplemd::configurations::MolecularDynamicsConfiguration mdConfig;
    tarch::configuration::ParseConfiguration::parseConfiguration<simplemd::configurations::MolecularDynamicsConfiguration>(_configFilename,
                                                                                                                           "molecular-dynamics", mdConfig);
    if (!mdConfig.isValid()) {
      std::cout << "ERROR precice::CouetteScenario: Invalid SimpleMD config!" << std::endl;
      exit(EXIT_FAILURE);
    }
    coupling::configurations::MaMiCoConfiguration<3> mamicoConfig;
    tarch::configuration::ParseConfiguration::parseConfiguration<coupling::configurations::MaMiCoConfiguration<dim>>(_configFilename, "mamico",
                                                                                                                     mamicoConfig);
    if (!mamicoConfig.isValid()) {
      std::cout << "ERROR precice::CouetteScenario: Invalid MaMiCo config!" << std::endl;
      exit(EXIT_FAILURE);
    }
    const coupling::configurations::CouetteConfig cfg = coupling::configurations::CouetteConfig::parseCouetteConfiguration(_configFilename);
    if (cfg.maSolverType != coupling::configurations::CouetteConfig::COUETTE_PRECICE) {
      std::cout << "ERROR precice::CouetteScenario: <macroscopic-solver type> must be \"precice\"" << std::endl;
      exit(EXIT_FAILURE);
    }

#if defined(LS1_MARDYN)
    auto offset = mdConfig.getDomainConfiguration().getGlobalDomainOffset();
    coupling::interface::LS1StaticCommData::getInstance().setConfigFilename("ls1config.xml");
    coupling::interface::LS1StaticCommData::getInstance().setBoxOffsetAtDim(0, offset[0]); // temporary till ls1 offset is natively supported
    coupling::interface::LS1StaticCommData::getInstance().setBoxOffsetAtDim(1, offset[1]);
    coupling::interface::LS1StaticCommData::getInstance().setBoxOffsetAtDim(2, offset[2]);
#endif

    _multiMDService = new tarch::utils::MultiMDService<dim>(mdConfig.getMPIConfiguration().getNumberOfProcesses(), cfg.totalNumberMDSimulations);

    _instanceHandling = new coupling::InstanceHandling<MY_LINKEDCELL, 3>(mdConfig, mamicoConfig, *_multiMDService);

    // MD equilibration steps
    unsigned int mdStepCounter = 0;
    _instanceHandling->switchOffCoupling();
    _instanceHandling->equilibrate(cfg.equSteps, mdStepCounter);
    _instanceHandling->switchOnCoupling();
    mdStepCounter += cfg.equSteps;

    _instanceHandling->setMDSolverInterface();

    // Initialize domain and physical parameters
    const tarch::la::Vector<3, double> domainOffset{mdConfig.getDomainConfiguration().getGlobalDomainOffset()};
    const tarch::la::Vector<3, double> cellSize{mamicoConfig.getCouplingCellConfiguration().getCouplingCellSize()};
    const tarch::la::Vector<dim, double> domainSize{mdConfig.getDomainConfiguration().getGlobalDomainSize()};
    const double massCell = cfg.density * cellSize[0] * cellSize[1] * cellSize[2];

    // Create the preCICE interface to enable communication between MaMiCo and preCICE
    _preciceInterface = new CouettePreciceInterface(massCell, cfg.twoWayCoupling);

    // init indexing
    coupling::indexing::IndexingService<3>::getInstance().initWithMDSize(
        mdConfig.getDomainConfiguration().getGlobalDomainSize(), mdConfig.getDomainConfiguration().getGlobalDomainOffset(),
        mdConfig.getMPIConfiguration().getNumberOfProcesses(), mamicoConfig.getCouplingCellConfiguration().getCouplingCellSize(),
        mamicoConfig.getParallelTopologyConfiguration().getParallelTopologyType(), mamicoConfig.getMomentumInsertionConfiguration().getInnerOverlap(),
        (unsigned int)rank);

    _multiMDCellService = new coupling::services::MultiMDCellService<MY_LINKEDCELL, dim>(
        _instanceHandling->getMDSolverInterface(), _preciceInterface, mdConfig, mamicoConfig, _configFilename.c_str(), *_multiMDService);
    // init filtering for all md instances
    _multiMDCellService->constructFilterPipelines();

    coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL, dim>::getInstance().setMacroscopicSolverInterface(_preciceInterface);

    _instanceHandling->setCouplingCellServices(*_multiMDCellService);

    _multiMDCellService->computeAndStoreTemperature(cfg.temp);

    // Create the MaMiCo/preCICE adapter and initialie it, i.e.:
    //  - create and communicate the preCICE vertices corresponding to the MaMiCo cartesian grid
    //  - create the CouplingCellContainers used as buffer for communication between MaMiCo and preCICE
    _preciceAdapter = new coupling::preciceadapter::PreciceAdapter<MY_LINKEDCELL, dim>(_preciceInterface, cfg.preciceConfigFile);
    _preciceAdapter->initialize();

    // Start the preCICE solving until the coupling is on going, i.e. until the final time is reached
    double mamico_dt = mdConfig.getSimulationConfiguration().getNumberOfTimesteps() * mdConfig.getSimulationConfiguration().getDt();
    int cycle = 0;
    while (_preciceAdapter->isCouplingOngoing()) {
      // In case of implicit coupling scheme, it might be required to save the current state
      if (_preciceAdapter->requiresWritingCheckpoint()) {
        std::cout << "Writing checkpoint" << std::endl;
        _instanceHandling->writeCheckpoint("checkpoint_restart_" + std::to_string(rank), 0);
      }

      _preciceAdapter->readData();
      _preciceAdapter->sendFromMacro2MD(_multiMDCellService);

      // Subcycle MD solving by a number of timesteps given in the config file
      _instanceHandling->simulateTimesteps(mdConfig.getSimulationConfiguration().getNumberOfTimesteps(), mdStepCounter, *_multiMDCellService);
      mdStepCounter += mdConfig.getSimulationConfiguration().getNumberOfTimesteps();
      _preciceAdapter->sendFromMD2Macro(_multiMDCellService);
      if (_preciceInterface->twoWayCoupling()) _preciceAdapter->writeData();

      // Write csv file if needed
      _multiMDCellService->plotEveryMacroscopicTimestep(cycle);

      // Send all data to CFD and advance the preCICE/MaMiCo adapter
      _preciceAdapter->advance(mamico_dt);

      // In case of implicit coupling scheme, it might be required to restart the coupling cycle
      if (_preciceAdapter->requiresReadingCheckpoint()) {
        std::cout << "Reading checkpoint" << std::endl;
      } else {
        cycle++;
        if (cfg.csvEveryTimestep >= 1 && cycle % cfg.csvEveryTimestep == 0)
          _preciceAdapter->write2csv(cycle);
      }
    }
  }

  void init() override { throw std::runtime_error("not supported yet"); }
  void runOneCouplingCycle(int cycle) override { throw std::runtime_error("not supported yet"); }
  coupling::solvers::AbstractCouetteSolver<3>* getSolver() override { throw std::runtime_error("not supported yet"); }

private:
  class CouettePreciceInterface : public coupling::preciceadapter::PreciceInterface<3> {
  private:
    const std::string _macro2MDMeshName;
    const std::string _MD2macroMeshName;
    const double _massCell;

  public:
    CouettePreciceInterface(const double massCell, const bool twoWayCoupling)
        : coupling::preciceadapter::PreciceInterface<3>(twoWayCoupling), _macro2MDMeshName("mamico-macro2MD-mesh"), _MD2macroMeshName("mamico-MD2macro-mesh"), _massCell{massCell} {
      coupling::preciceadapter::PreciceInterface<3>::addDataDescription(
        _macro2MDMeshName, coupling::preciceadapter::DataDescription{"VelocityMacro", coupling::preciceadapter::DataType::vector});
      coupling::preciceadapter::PreciceInterface<3>::addDataDescription(
        _MD2macroMeshName, coupling::preciceadapter::DataDescription{"VelocityMicro", coupling::preciceadapter::DataType::vector});
    }

    unsigned int getOuterRegion() override {
      return 3;
    }
#ifdef MAMICO_MD2MACRO_LAYER3
    // At hole inset 0 the CFD face lies between the innermost macro2MD layer and
    // the first md2macro cell. Putting that layer on the MD2macro mesh, carrying
    // the continuum velocity it was itself given, turns linear-cell-interpolation
    // on the face into (CFD_own + MD)/2: the IcoFoamInterface operator.
    bool isMD2Macro(I01 idx) const override { return I12::contains(idx) || isInnermostMacro2MDLayer(idx); }
    static bool isInnermostMacro2MDLayer(I01 idx) {
      if (!I08::contains(idx) || I12::contains(idx)) return false;
      const auto lo = I12::lowerBoundary.get(), hi = I12::upperBoundary.get(), v = idx.get();
      for (unsigned int d = 0; d < 3; d++)
        if (v[d] < lo[d] - 1 || v[d] > hi[d] + 1) return false;
      return true;
    }
    mutable std::map<unsigned int, tarch::la::Vector<3, double>> _continuumVelocity;
#endif

    std::vector<unsigned int> getRanks(I01 idx) override { return {0}; }

    std::string getMacro2MDSolverMeshName(I01 idx) const override {
      return _macro2MDMeshName;
    }

    std::string getMD2MacroSolverMeshName(I01 idx) const override {
      return _MD2macroMeshName;
    }

    bool contains(std::string meshName, I01 idx) const override {
      if (meshName == _macro2MDMeshName) {
        return isMacro2MD(idx);
      } else if (meshName == _MD2macroMeshName) {
        return isMD2Macro(idx);
      }
      return false;
    }

    void readVectorData(const std::string& meshName, const std::string& dataName, coupling::datastructures::CouplingCell<3>* const cell,
    const I01& idx, const double& vx, const double& vy, const double& vz) const override {
#ifdef MAMICO_MD2MACRO_LAYER3
      _continuumVelocity[I00{idx}.get()] = tarch::la::Vector<3, double>{vx, vy, vz};
#endif
      tarch::la::Vector<3, double> momentum{vx, vy, vz};
      momentum=momentum*_massCell;
      cell->setMicroscopicMass(_massCell);
      cell->setMicroscopicMomentum(momentum);
    }

    void writeVectorData(const std::string& meshName, const std::string& dataName, const coupling::datastructures::CouplingCell<3>* const cell,
    const I01& idx, double& vx, double& vy, double& vz) const override {
#ifdef MAMICO_MD2MACRO_LAYER3
      if (!I12::contains(idx)) {
        const auto& u = _continuumVelocity.at(I00{idx}.get());
        vx = u[0]; vy = u[1]; vz = u[2];
        return;
      }
#endif
      tarch::la::Vector<3, double> velocity(0.0);
      if (cell->getMacroscopicMass() != 0.0) {
        velocity = (1.0 / cell->getMacroscopicMass()) * cell->getMacroscopicMomentum();
      }
      vx=velocity[0];
      vy=velocity[1];
      vz=velocity[2];
    }
  };

  const std::string _configFilename;
  coupling::preciceadapter::PreciceAdapter<MY_LINKEDCELL, 3>* _preciceAdapter = nullptr;
  coupling::preciceadapter::PreciceInterface<3>* _preciceInterface = nullptr;
  coupling::InstanceHandling<MY_LINKEDCELL, 3>* _instanceHandling = nullptr;
  tarch::utils::MultiMDService<3>* _multiMDService = nullptr;
  coupling::services::MultiMDCellService<MY_LINKEDCELL, 3>* _multiMDCellService = nullptr;
};
