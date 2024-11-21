#pragma once

#include "coupling/CouplingMDDefinitions.h"
#include "coupling/configurations/MaMiCoConfiguration.h"
#include "coupling/indexing/IndexingService.h"
#include "coupling/interface/MDSimulationFactory.h"
#include "coupling/scenario/Scenario.h"
#include "coupling/services/MultiMDCellService.h"
#include "coupling/solvers/CouetteSolverInterface.h"
#include "coupling/scenario/precice/PreciceAdapter.h"
#include "simplemd/configurations/MolecularDynamicsConfiguration.h"
#include "tarch/configuration/Configuration.h"
#include "tarch/configuration/ParseConfiguration.h"
#include "tarch/utils/MultiMDService.h"
#include "coupling/scenario/precice/PreciceInterface.h"
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
#include <mpi.h>
#endif
#include <random>
#include <sys/time.h>

#if defined(LS1_MARDYN)
#include "coupling/interface/impl/ls1/LS1MDSolverInterface.h"
#include "coupling/interface/impl/ls1/LS1StaticCommData.h"
#include "utils/Logger.h"
#include "plugins/NEMD/MettDeamonFeedrateDirector.h"
#include "plugins/PluginBase.h"
using Log::global_log;
#endif

namespace coupling {
namespace scenario {
namespace precice {
  class EvaporationScenario;
}
}
}


class coupling::scenario::precice::EvaporationScenario : public Scenario {
public:
  EvaporationScenario() : Scenario("evaporation") {}
  ~EvaporationScenario() {
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
    if (_preciceInterface != NULL) {
      delete _preciceInterface;
      _preciceInterface = NULL;
    }
    if (_multiMDCellService != NULL) {
      delete _multiMDCellService;
      _multiMDCellService = NULL;
    }
  }

  void run() override {
    #if !defined(LS1_MARDYN)
    throw new std::runtime_error("EvaporationScenario requires LS1 MD solver");
    #endif

    const unsigned int dim = 3;
    int rank;
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
    Log::global_log = std::make_unique<Log::Logger>(Log::Error); //Log::Info
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
    global_log->set_mpi_output_root(0);
#endif

    // Read the configurations
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

    auto offset = mdConfig.getDomainConfiguration().getGlobalDomainOffset();
    coupling::interface::LS1StaticCommData::getInstance().setConfigFilename("ls1config.xml");
    coupling::interface::LS1StaticCommData::getInstance().setBoxOffsetAtDim(0, offset[0]); // temporary till ls1 offset is natively supported
    coupling::interface::LS1StaticCommData::getInstance().setBoxOffsetAtDim(1, offset[1]);
    coupling::interface::LS1StaticCommData::getInstance().setBoxOffsetAtDim(2, offset[2]);

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
    const tarch::la::Vector<3, double> cellSize{mamicoConfig.getCouplingCellConfiguration().getCouplingCellSize()};
    const tarch::la::Vector<dim, double> domainSize{mdConfig.getDomainConfiguration().getGlobalDomainSize()};
    tarch::la::Vector<3, int> numberCells;
    for (unsigned int d = 0; d < 3; d++) {
      numberCells[d] = floor(domainSize[d] / cellSize[d] + 0.5);
    }
    const double densityCell = mdConfig.getDomainConfiguration().getMoleculesPerDirection()[0] *
                               mdConfig.getDomainConfiguration().getMoleculesPerDirection()[1] *
                               mdConfig.getDomainConfiguration().getMoleculesPerDirection()[2] / (domainSize[0] * domainSize[1] * domainSize[2]);
    const double massCell = densityCell * cellSize[0] * cellSize[1] * cellSize[2];

    // Create the preCICE interface to enable communication between MaMiCo and preCICE
    _preciceInterface = new EvaporationPreciceInterface(massCell); 

    // init indexing
    coupling::indexing::IndexingService<3>::getInstance().initWithMDSize(
        mdConfig.getDomainConfiguration().getGlobalDomainSize(), mdConfig.getDomainConfiguration().getGlobalDomainOffset(),
        mdConfig.getMPIConfiguration().getNumberOfProcesses(), mamicoConfig.getCouplingCellConfiguration().getCouplingCellSize(),
        mamicoConfig.getParallelTopologyConfiguration().getParallelTopologyType(), mamicoConfig.getMomentumInsertionConfiguration().getInnerOverlap(),
        (unsigned int)rank);

    _multiMDCellService = new coupling::services::MultiMDCellService<MY_LINKEDCELL, dim>(
        _instanceHandling->getMDSolverInterface(), _preciceInterface, mdConfig, mamicoConfig, xmlConfigurationFilename.c_str(), *_multiMDService);
    // init filtering for all md instances
    _multiMDCellService->constructFilterPipelines();
    
    coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL, dim>::getInstance().setMacroscopicSolverInterface(_preciceInterface);

    _instanceHandling->setCouplingCellServices(*_multiMDCellService);

    _multiMDCellService->computeAndStoreTemperature(scenarioConfig.temp);

    // Create the MaMiCo/preCICE adapter and initialie it, i.e.:
    //  - create and communicate the preCICE vertices corresponding to the MaMiCo cartesian grid
    //  - create the CouplingCellContainers used as buffer for communication between MaMiCo and preCICE
    _preciceAdapter = new coupling::preciceadapter::PreciceAdapter<MY_LINKEDCELL, dim>(_preciceInterface);
    _preciceAdapter->initialize();

    // Start the preCICE solving until the coupling is on going, i.e. until the final time is reached
    double mamico_dt = mdConfig.getSimulationConfiguration().getNumberOfTimesteps() * mdConfig.getSimulationConfiguration().getDt();
    int cycle = 0;

    MettDeamonFeedrateDirector* mettDeamonFeedrateDirector = nullptr;
    std::list<PluginBase*>& plugins = *(global_simulation->getPluginList() );
    for (auto&& pit:plugins) {
      std::string name = pit->getPluginName();
      if(name == "MettDeamonFeedrateDirector")
        mettDeamonFeedrateDirector = dynamic_cast<MettDeamonFeedrateDirector*>(pit);
    }
    int updateFrequency = mettDeamonFeedrateDirector->getUpdateFreq();
    dynamic_cast<EvaporationPreciceInterface*>(_preciceInterface)->updateFeedrate(mettDeamonFeedrateDirector->getInitFeedrate());
    dynamic_cast<EvaporationPreciceInterface*>(_preciceInterface)->updateFeedrate(0.0);
   
    while (_preciceAdapter->isCouplingOngoing()) {
      // In case of implicit coupling scheme, it might be required to save the current state
      if (_preciceAdapter->requiresWritingCheckpoint()) {
        std::cout << "Writing checkpoint" << std::endl;
        _instanceHandling->writeCheckpoint("checkpoint_restart_" + rank, 0);
      }

      _preciceAdapter->readData();
      _preciceAdapter->sendFromMacro2MD(_multiMDCellService);

      // Subcycle MD solving by a number of timesteps given in the config file
      _instanceHandling->simulateTimesteps(mdConfig.getSimulationConfiguration().getNumberOfTimesteps(), mdStepCounter, *_multiMDCellService);
      mdStepCounter += mdConfig.getSimulationConfiguration().getNumberOfTimesteps();
      if (mdStepCounter % updateFrequency == 0) {
        double feedrate = mettDeamonFeedrateDirector->getFeedrate();
        dynamic_cast<EvaporationPreciceInterface*>(_preciceInterface)->updateFeedrate(mettDeamonFeedrateDirector->getFeedrate());
        if (rank==0) std::cout << "updating CFD feedrate given by MDFD :" << feedrate << std::endl;
      }

      _preciceAdapter->sendFromMD2Macro(_multiMDCellService);
      _preciceAdapter->writeData();

      // Write csv file if needed
      _multiMDCellService->plotEveryMacroscopicTimestep(cycle);

      // Send all data to CFD and advance the preCICE/MaMiCo adapter
      _preciceAdapter->advance(mamico_dt);

      // In case of implicit coupling scheme, it might be required to restart the coupling cycle
      if (_preciceAdapter->requiresReadingCheckpoint()) {
        std::cout << "Reading checkpoint" << std::endl;
      } else {
        cycle++;
        if (scenarioConfig.csvEveryTimestep >= 1 && cycle % scenarioConfig.csvEveryTimestep == 0)
          _preciceAdapter->write2csv(cycle);
      } 
    }
  }

  void init() override { throw std::runtime_error("not supported yet"); }
  void runOneCouplingCycle(int cycle) override { throw std::runtime_error("not supported yet"); }
  coupling::solvers::AbstractCouetteSolver<3>* getSolver() override { throw std::runtime_error("not supported yet"); }

private:
  class EvaporationPreciceInterface : public coupling::preciceadapter::PreciceInterface<3> {
  private:
    const double _massCell;
    const std::string _M2mMeshName;
    const std::string _m2MLMeshName;
    const std::string _m2MVMeshName;
    double _feedrate;

  public:
    EvaporationPreciceInterface(const double massCell)
        : coupling::preciceadapter::PreciceInterface<3>(true), 
        _massCell{massCell}, 
        _M2mMeshName("mamico-M2m-mesh"), _m2MLMeshName("mamico-m2ML-mesh"), _m2MVMeshName("mamico-m2MV-mesh"), _feedrate(0.0) {
        // macro to micro
        coupling::preciceadapter::PreciceInterface<3>::addDataDescription(
      _M2mMeshName, coupling::preciceadapter::DataDescription{"VelocityMacro", coupling::preciceadapter::DataType::vector});
        // micro to macro in the liquid
          coupling::preciceadapter::PreciceInterface<3>::addDataDescription(
      _m2MLMeshName, coupling::preciceadapter::DataDescription{"VelocityMicro", coupling::preciceadapter::DataType::vector});
          // micro to macro vapor to liquid
            coupling::preciceadapter::PreciceInterface<3>::addDataDescription(
      _m2MVMeshName, coupling::preciceadapter::DataDescription{"VelocityMicro", coupling::preciceadapter::DataType::vector});      
    }

    void updateFeedrate(double feedrate) { _feedrate = feedrate; }

    unsigned int getOuterRegion() override {
      return 3;
    }

    std::vector<unsigned int> getRanks(I01 idx) override { return {0}; }

    std::string getMacro2MDSolverMeshName(I01 idx) const override {
      return _M2mMeshName;
    }

    bool isMD2MacroLiquid(I01 idx) const {
      return idx.get()[1] == 4;
    }

    bool isMD2MacroVapor(I01 idx) const {
      return idx.get()[1] - 1 == coupling::indexing::IndexingService<3>::getInstance().getGlobalMDDomainSize()[1];
    }

    bool isMD2Macro(I01 idx) const override {
      return (isMD2MacroLiquid(idx) || isMD2MacroVapor(idx));
    }

    std::string getMD2MacroSolverMeshName(I01 idx) const override {
      std::string meshName;
      if (isMD2MacroLiquid(idx)) {
        meshName = _m2MLMeshName;
      } else if (isMD2MacroVapor(idx)) {
        meshName = _m2MVMeshName;
      } 
      return meshName;
    }

    tarch::la::Vector<3, double> getMD2MacroSolverMeshOffset(I01 idx) const override {
      tarch::la::Vector<3, double> offset;
      if (isMD2MacroVapor(idx)) {
        offset = tarch::la::Vector<3,double>{0, /*MD size*/-115/*CFD size*/-15/*half a mamico cell*/+1.25, 0};
      }
      return offset;
    }

    bool contains(std::string meshName, I01 idx) const override {
      bool contained = false;
      if (meshName == _M2mMeshName) {
        contained = isMacro2MD(idx);
      } else if (meshName == _m2MLMeshName) {
        contained = isMD2MacroLiquid(idx);
      } else if (meshName == _m2MVMeshName) {
        contained = isMD2MacroVapor(idx);
      }
      return contained;
    }

    void readVectorData(const std::string& meshName, const std::string& dataName, coupling::datastructures::CouplingCell<3>* const cell, 
    const I01& idx, const double& vx, const double& vy, const double& vz) const override {
      tarch::la::Vector<3, double> momentum{vx, vy, vz};
      momentum=momentum*_massCell;
      cell->setMicroscopicMass(_massCell);
      cell->setMicroscopicMomentum(momentum);
    }

    virtual void writeVectorData(const std::string& meshName, const std::string& dataName, const coupling::datastructures::CouplingCell<3>* const cell, 
    const I01& idx, double& vx, double& vy, double& vz) const override {
      if (meshName == _m2MLMeshName) {
        tarch::la::Vector<3, double> velocity;
        if (cell->getMacroscopicMass() != 0.0) {
          velocity = (1.0 / cell->getMacroscopicMass()) * cell->getMacroscopicMomentum();
        }
        vx=velocity[0];
        vy=velocity[1];
        vz=velocity[2];
      } else if (meshName == _m2MVMeshName) {
        vx = 0.0;
        vy = _feedrate;
        vz = 0.0;
      }
    }
  };

  struct ScenarioConfig : public tarch::configuration::Configuration {
    ~ScenarioConfig() {}

    void parseSubtag(tinyxml2::XMLElement* node) override {
      tinyxml2::XMLElement* subtag = node->FirstChildElement("coupling");
      tarch::configuration::ParseConfiguration::readIntMandatory(csvEveryTimestep, subtag, "write-csv-every-timestep");
      subtag = node->FirstChildElement("microscopic-solver");
      tarch::configuration::ParseConfiguration::readDoubleMandatory(temp, subtag, "temperature");
      tarch::configuration::ParseConfiguration::readIntMandatory(equSteps, subtag, "equilibration-steps");
      tarch::configuration::ParseConfiguration::readIntMandatory(totalNumberMDSimulations, subtag, "number-md-simulations");
    };

    std::string getTag() const override { return "scenario"; };

    bool isValid() const override { return true; };

    int csvEveryTimestep;
    bool twoWayCoupling;
    int equSteps;
    double temp;
    int totalNumberMDSimulations;
  };

  coupling::preciceadapter::PreciceAdapter<MY_LINKEDCELL, 3>* _preciceAdapter;
  coupling::preciceadapter::PreciceInterface<3>* _preciceInterface;
  coupling::InstanceHandling<MY_LINKEDCELL, 3>* _instanceHandling;
  tarch::utils::MultiMDService<3>* _multiMDService;
  coupling::services::MultiMDCellService<MY_LINKEDCELL, 3>* _multiMDCellService;
};
