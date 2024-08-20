#pragma once

#include "coupling/CouplingMDDefinitions.h"
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
  CouetteScenario() : Scenario("couette_precice") {}
  ~CouetteScenario() {
    if (_instanceHandling != nullptr) {
      delete _instanceHandling;
    }
    if (_multiMDService != NULL) {
      delete _multiMDService;
      _multiMDService = NULL;
    }
    if (_preciceInterface != NULL) {
      delete _preciceInterface;
      _preciceInterface = NULL;
    }
    if (_preciceAdapter != NULL) {
      delete _preciceAdapter;
      _preciceAdapter = NULL;
    }
    if (_multiMDCellService != NULL) {
      delete _multiMDCellService;
      _multiMDCellService = NULL;
    }
  }

  /**
   * Run the Couette scenario for a certain amount of timesteps given in the configuration file.
   * PreCICE is steering the solving between the CFD solver and MaMiCo/MD solver.
   */
  void run () override {
    const unsigned int dim = 3;
    int rank;
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
#if defined(LS1_MARDYN)
    global_log = new Log::Logger(Log::Info); // Info
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
    global_log->set_mpi_output_root(0);
#endif
#endif
    std::string xmlConfigurationFilename("config.xml");
    // Read the MD config
    simplemd::configurations::MolecularDynamicsConfiguration mdConfig;
    tarch::configuration::ParseConfiguration::parseConfiguration<simplemd::configurations::MolecularDynamicsConfiguration>(xmlConfigurationFilename,
                                                                                                                           "molecular-dynamics", mdConfig);
    // Read the MaMiCo config
    coupling::configurations::MaMiCoConfiguration<3> mamicoConfig;
    tarch::configuration::ParseConfiguration::parseConfiguration<coupling::configurations::MaMiCoConfiguration<dim>>(xmlConfigurationFilename, "mamico",
                                                                                                                     mamicoConfig);
    // Read specific Couette config
    ScenarioConfig scenarioConfig;
    tarch::configuration::ParseConfiguration::parseConfiguration<ScenarioConfig>(xmlConfigurationFilename, "scenario", scenarioConfig);
#if defined(LS1_MARDYN)
    auto offset = mdConfig.getDomainConfiguration().getGlobalDomainOffset();
    coupling::interface::LS1StaticCommData::getInstance().setConfigFilename("ls1config.xml");
    coupling::interface::LS1StaticCommData::getInstance().setBoxOffsetAtDim(0, offset[0]); // temporary till ls1 offset is natively supported
    coupling::interface::LS1StaticCommData::getInstance().setBoxOffsetAtDim(1, offset[1]);
    coupling::interface::LS1StaticCommData::getInstance().setBoxOffsetAtDim(2, offset[2]);
#endif
    _multiMDService = new tarch::utils::MultiMDService<dim>(mdConfig.getMPIConfiguration().getNumberOfProcesses(), scenarioConfig.totalNumberMDSimulations);

    _instanceHandling = new coupling::InstanceHandling<MY_LINKEDCELL, 3>(mdConfig, mamicoConfig, *_multiMDService);

    // MD equilibration steps
    unsigned int mdStepCounter = 0;
    _instanceHandling->switchOffCoupling();
    _instanceHandling->equilibrate(scenarioConfig.equSteps, mdStepCounter);
    _instanceHandling->switchOnCoupling();
    mdStepCounter += scenarioConfig.equSteps;

    _instanceHandling->setMDSolverInterface();

    // Initialize domain and physical parameters
    const tarch::la::Vector<3, double> domainOffset{mdConfig.getDomainConfiguration().getGlobalDomainOffset()};
    const tarch::la::Vector<3, double> cellSize{mamicoConfig.getMacroscopicCellConfiguration().getMacroscopicCellSize()};
    const tarch::la::Vector<dim, double> domainSize{mdConfig.getDomainConfiguration().getGlobalDomainSize()};
    tarch::la::Vector<3, int> numberCells;
    for (unsigned int d = 0; d < 3; d++) {
      numberCells[d] = floor(domainSize[d] / cellSize[d] + 0.5);
    }
    const unsigned int overLap = mamicoConfig.getMomentumInsertionConfiguration().getInnerOverlap();
    const double densityCell = mdConfig.getDomainConfiguration().getMoleculesPerDirection()[0] *
                               mdConfig.getDomainConfiguration().getMoleculesPerDirection()[1] *
                               mdConfig.getDomainConfiguration().getMoleculesPerDirection()[2] / (domainSize[0] * domainSize[1] * domainSize[2]);
    const double massCell = densityCell * cellSize[0] * cellSize[1] * cellSize[2];

    // Create the preCICE interface to enable communication between MaMiCo and preCICE
    _preciceInterface = new PreciceInterface(numberCells, overLap, massCell, scenarioConfig.twoWayCoupling);

    // Initialize the indexing service
    coupling::indexing::IndexingService<dim>::getInstance().init(mdConfig, mamicoConfig, _preciceInterface, rank);

    _multiMDCellService = new coupling::services::MultiMDCellService<MY_LINKEDCELL, dim>(
        _instanceHandling->getMDSolverInterface(), _preciceInterface, mdConfig, mamicoConfig, xmlConfigurationFilename.c_str(), *_multiMDService);
    
    coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL, dim>::getInstance().setMacroscopicSolverInterface(_preciceInterface);

    _instanceHandling->setMacroscopicCellServices(*_multiMDCellService);

    _multiMDCellService->computeAndStoreTemperature(scenarioConfig.temp);

    // Create the MaMiCo/preCICE adapter and initialie it, i.e.:
    //  - create and communicate the preCICE vertices corresponding to the MaMiCo cartesian grid
    //  - create the CouplingCellContainers used as buffer for communication between MaMiCo and preCICE
    _preciceAdapter = new coupling::preciceadapter::PreciceAdapter<dim>();
    _preciceAdapter->setMeshes(_preciceInterface, domainOffset, cellSize);
    _preciceAdapter->initialize();

    // Start the preCICE solving until the coupling is on going, i.e. until the final time is reached
    double mamico_dt = mdConfig.getSimulationConfiguration().getNumberOfTimesteps() * mdConfig.getSimulationConfiguration().getDt();
    int cycle = 0;
    while (_preciceAdapter->isCouplingOngoing()) {
      // In case of implicit coupling scheme, it might be required to save the current state
      if (_preciceAdapter->requiresWritingCheckpoint()) {
        std::cout << "Writing checkpoint" << std::endl;
        _instanceHandling->writeCheckpoint("checkpoint_restart_" + rank, 0);
      }

      // Read the data send from CFD to the preCICE M2m vertices and store it in the preCICE adapter M2m coupling cell container
      _preciceAdapter->readData(_preciceInterface);

      // Send the data stored in the preCICE M2m coupling cell container to the MaMiCo M2m coupling cell container
      _preciceAdapter->sendFromMacro2MD(_multiMDCellService);


      if (!scenarioConfig.couetteAnalytical) { // Non-analytical MD solver
        // Subcycle MD solving by a number of timesteps given in the config file
        _instanceHandling->simulateTimesteps(mdConfig.getSimulationConfiguration().getNumberOfTimesteps(), mdStepCounter, *_multiMDCellService);
        mdStepCounter += mdConfig.getSimulationConfiguration().getNumberOfTimesteps();

        // Send the data contained in the MaMiCo m2M coupling cell container to the preCICE m2M coupling cell container
        _preciceAdapter->sendFromMD2Macro(_multiMDCellService);

        // Write csv file if needed
        _multiMDCellService->plotEveryMacroscopicTimestep(cycle);
      } else { // Analytical Couette solver, usefull for debugging
        coupling::solvers::CouetteSolver<3>* couetteSolver =
            new coupling::solvers::CouetteSolver<3>(scenarioConfig.channelHeight, scenarioConfig.wallVelocity, scenarioConfig.kinematicViscosity);
        couetteSolver->advance(mamico_dt * (cycle + 1));
        for (auto pair : _preciceAdapter->getm2MCells()) {
          I01 idx;
          coupling::datastructures::CouplingCell<3>* couplingCell;
          std::tie(couplingCell, idx) = pair;
          auto midPoint = idx.getCellMidPoint();
          const auto momentum{massCell * (couetteSolver->getVelocity(cellMidPoint))};
          _preciceAdapter->getm2MCells()[i]->setMacroscopicMass(massCell);
          _preciceAdapter->getm2MCells()[i]->setMacroscopicMomentum(momentum);
        }
      }

      // Write the data stored in the preCICE m2M coupling cell container to the preCICE m2M vertices 
      _preciceAdapter->writeData(_preciceInterface);

      // Send all data to CFD and advance the preCICE/MaMiCo adapter
      _preciceAdapter->advance(mamico_dt);

      // In case of implicit coupling scheme, it might be required to restart the coupling cycle
      if (_preciceAdapter->requiresReadingCheckpoint()) {
        std::cout << "Reading checkpoint" << std::endl;
      } else {
        cycle++;
        if (!_preciceAdapter->getm2MCells().empty() && scenarioConfig.csvEveryTimestep >= 1 && cycle % scenarioConfig.csvEveryTimestep == 0)
          write2CSV(_preciceAdapter->getm2MCells(), cycle, rank);
      } 
    }
  }

  void init() override { throw std::runtime_error("not supported yet"); }
  void runOneCouplingCycle(int cycle) override { throw std::runtime_error("not supported yet"); }
  coupling::solvers::AbstractCouetteSolver<3>* getSolver() override { throw std::runtime_error("not supported yet"); }

private:

  void write2CSV(coupling::datastructures::FlexibleCellContainer<3>& m2MCells, int couplingCycle, unsigned int rank) {
    std::stringstream ss;
    ss << "results_" << rank << "_" << couplingCycle << ".csv";
    std::ofstream file(ss.str().c_str());
    if (!file.is_open()) {
      exit(EXIT_FAILURE);
    }
    file << "i;j;k;x;y;z;v_x;v_y;v_z;T;m" << std::endl;
    for (auto pair : m2MCells) {
      I01 idx;
      if (I12.contains(idx)) {
        coupling::datastructures::CouplingCell<3>* couplingCell;
        std::tie(couplingCell, idx) = pair;
        auto cellMidPoint = idx.getCellMidPoint();
        tarch::la::Vector<3, double> vel(couplingCell->getMacroscopicMomentum());
        if (couplingCell->getMacroscopicMass() != 0.0) {
          vel = (1.0 / couplingCell->getMacroscopicMass()) * vel;
        }
        file   << idx.get()[0] << ";" << idx.get()[1] << ";" << idx.get()[2] << ";" 
               << cellMidPoint[0] << ";" << cellMidPoint[1] << ";" << cellMidPoint[2] << ";" 
               << vel[0] << ";" << vel[1] << ";" << vel[2] << ";" 
               << couplingCell->getTemperature() << ";" << couplingCellgetMacroscopicMass();
        file   << std::endl;
      }
    }
    file.close();
  }

  struct ScenarioConfig : public tarch::configuration::Configuration {
    ~ScenarioConfig() {}

    void parseSubtag(tinyxml2::XMLElement* node) override {
      tinyxml2::XMLElement* subtag = node->FirstChildElement("coupling");
      tarch::configuration::ParseConfiguration::readIntMandatory(csvEveryTimestep, subtag, "write-csv-every-timestep");
      tarch::configuration::ParseConfiguration::readBoolMandatory(twoWayCoupling, subtag, "two-way-coupling");
      subtag = node->FirstChildElement("microscopic-solver");
      tarch::configuration::ParseConfiguration::readDoubleMandatory(temp, subtag, "temperature");
      tarch::configuration::ParseConfiguration::readIntMandatory(equSteps, subtag, "equilibration-steps");
      tarch::configuration::ParseConfiguration::readIntMandatory(totalNumberMDSimulations, subtag, "number-md-simulations");
      tarch::configuration::ParseConfiguration::readBoolOptional(couetteAnalytical, subtag, "couette-analytical");
      tarch::configuration::ParseConfiguration::readDoubleOptional(channelHeight, subtag, "channel-height");
      tarch::configuration::ParseConfiguration::readDoubleOptional(wallVelocity, subtag, "wall-velocity");
      tarch::configuration::ParseConfiguration::readDoubleOptional(kinematicViscosity, subtag, "kinematic-viscosity");
    };

    std::string getTag() const override { return "scenario"; };

    bool isValid() const override { return true; };

    int csvEveryTimestep;
    int equSteps;
    double temp;
    int totalNumberMDSimulations;
    bool couetteAnalytical = false;
    double channelHeight;
    double wallVelocity;
    double kinematicViscosity;
    bool twoWayCoupling;
  };

  class PreciceInterface : public coupling::preciceadapter::PreciceInterface<3> {
  private:
    const tarch::la::Vector<3, unsigned int> _globalNumberMacroscopicCells;
    const unsigned int _overlap;
    const double _massCell;
    const std::string _M2mMeshName;
    const std::string _m2MMeshName;

  public:
    PreciceInterface(const tarch::la::Vector<3, int> globalNumberMacroscopicCells, const unsigned int overlap, const double massCell, const bool twoWayCoupling)
        : coupling::preciceadapter::PreciceInterface<3>(twoWayCoupling), 
        _globalNumberMacroscopicCells(globalNumberMacroscopicCells), _overlap(overlap), _massCell{massCell},
        _M2mMeshName("mamico-M2m-mesh"), _m2MMeshName("mamico-m2M-mesh") {
      coupling::preciceadapter::PreciceInterface<3>::addData(
        _M2mMeshName, coupling::preciceadapter::Data{"VelocityMacro", coupling::preciceadapter::DataType::vector});
      coupling::preciceadapter::PreciceInterface<3>::addData(
        _m2MMeshName, coupling::preciceadapter::Data{"VelocityMicro", coupling::preciceadapter::DataType::vector});      
    }

    virtual bool receiveMacroscopicQuantityFromMDSolver(tarch::la::Vector<3, unsigned int> globalCellIndex) override {
      bool rcv = true;
      for (unsigned int currentDim = 0; currentDim < 3; currentDim++) {
        rcv &= globalCellIndex[currentDim] >= 1 + (_overlap - 1);
        rcv &= globalCellIndex[currentDim] < _globalNumberMacroscopicCells[currentDim] + 1 - (_overlap - 1);
      }
      return rcv;
    }

    virtual bool sendMacroscopicQuantityToMDSolver(tarch::la::Vector<3, unsigned int> globalCellIndex) override {
      bool isGhostCell = false;
      bool isInner = true;
      for (unsigned int currentDim = 0; currentDim < 3; currentDim++) {
        isGhostCell |= globalCellIndex[currentDim] > _globalNumberMacroscopicCells[currentDim];
        isGhostCell |= globalCellIndex[currentDim] < 1;
        isInner &= globalCellIndex[currentDim] >= 1 + _overlap;
        isInner &= globalCellIndex[currentDim] < _globalNumberMacroscopicCells[currentDim] + 1 - _overlap;
      }
      return (!isGhostCell) && (!isInner);
    }

    std::vector<unsigned int> getRanks(tarch::la::Vector<3, unsigned int> globalCellIndex) override { return {0}; }

    // std::vector<unsigned int> getSourceRanks(tarch::la::Vector<3, unsigned int> globalCellIndex) override { 
    //   coupling::indexing::CellIndex<3, coupling::indexing::IndexTrait::vector> cellIndex_v{static_cast<tarch::la::Vector<3,int>>(globalCellIndex)};
    //   std::vector<unsigned int> ranks = coupling::indexing::IndexingService<3>::getInstance().getRanksForGlobalIndex(cellIndex_v);
    //   return ranks;
    // }

    // std::vector<unsigned int> getTargetRanks(tarch::la::Vector<3, unsigned int> globalCellIndex) override {
    //   const unsigned int rank = coupling::indexing::IndexingService<3>::getInstance().getUniqueRankForGlobalIndex(globalCellIndex);
    //   std::vector<unsigned int> ranks;
    //   ranks.push_back(rank);
    //   return ranks; 
    // }

    std::string getMacroscopicToMDSolverMeshName(tarch::la::Vector<3, unsigned int> globalCellIndex) override {
      return _M2mMeshName;
    }

    std::string getMDToMacroscopicSolverMeshName(tarch::la::Vector<3, unsigned int> globalCellIndex) override {
      return _m2MMeshName;
    }

    void readVectorData(std::string meshName, std::string dataName, coupling::datastructures::MacroscopicCell<3>* const cell, const double vx, const double vy, const double vz) override {
      tarch::la::Vector<3, double> momentum{vx, vy, vz};
      momentum=momentum*_massCell;
      cell->setMicroscopicMass(_massCell);
      cell->setMicroscopicMomentum(momentum);
    }

    void writeVectorData(std::string meshName, std::string dataName, const coupling::datastructures::MacroscopicCell<3>* const cell, double& vx, double& vy, double& vz) {
      tarch::la::Vector<3, double> velocity;
      if (cell->getMacroscopicMass() != 0.0) {
        velocity = (1.0 / cell->getMacroscopicMass()) * cell->getMacroscopicMomentum();
      }
      vx=velocity[0];
      vy=velocity[1];
      vz=velocity[2];
    }
  };

  coupling::preciceadapter::PreciceAdapter<3>* _preciceAdapter;
  coupling::preciceadapter::PreciceInterface<3>* _preciceInterface;
  coupling::InstanceHandling<MY_LINKEDCELL, 3>* _instanceHandling;
  tarch::utils::MultiMDService<3>* _multiMDService;
  coupling::services::MultiMDCellService<MY_LINKEDCELL, 3>* _multiMDCellService;
};
