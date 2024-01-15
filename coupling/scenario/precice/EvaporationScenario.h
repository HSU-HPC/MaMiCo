#pragma once

#include "coupling/CouplingMDDefinitions.h"
#include "coupling/configurations/MaMiCoConfiguration.h"
#include "coupling/indexing/IndexingService.h"
#include "coupling/interface/MDSimulationFactory.h"
#include "coupling/interface/impl/SimpleMD/SimpleMDSolverInterface.h"
#include "coupling/scenario/Scenario.h"
#include "coupling/services/MultiMDCellService.h"
#include "coupling/solvers/CouetteSolver.h"
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
    if (_multiMDCellService != NULL) {
      delete _multiMDCellService;
      _multiMDCellService = NULL;
    }
  }

  void run() override {
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
    simplemd::configurations::MolecularDynamicsConfiguration mdConfig;
    tarch::configuration::ParseConfiguration::parseConfiguration<simplemd::configurations::MolecularDynamicsConfiguration>(xmlConfigurationFilename,
                                                                                                                           "molecular-dynamics", mdConfig);
    coupling::configurations::MaMiCoConfiguration<3> mamicoConfig;
    tarch::configuration::ParseConfiguration::parseConfiguration<coupling::configurations::MaMiCoConfiguration<dim>>(xmlConfigurationFilename, "mamico",
                                                                                                                     mamicoConfig);
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
    unsigned int mdStepCounter = 0;
    _instanceHandling->switchOffCoupling();
    _instanceHandling->equilibrate(scenarioConfig.equSteps, mdStepCounter);
    _instanceHandling->switchOnCoupling();
    mdStepCounter += scenarioConfig.equSteps;
    _instanceHandling->setMDSolverInterface();
    const tarch::la::Vector<3, double> domainOffset{mdConfig.getDomainConfiguration().getGlobalDomainOffset()};
    const tarch::la::Vector<3, double> cellSize{mamicoConfig.getMacroscopicCellConfiguration().getMacroscopicCellSize()};
    const tarch::la::Vector<dim, double> domainSize{mdConfig.getDomainConfiguration().getGlobalDomainSize()};
    tarch::la::Vector<3, int> numberCells;
    for (unsigned int d = 0; d < 3; d++) {
      numberCells[d] = floor(domainSize[d] / cellSize[d] + 0.5);
    }
    const unsigned int overLap = mamicoConfig.getMomentumInsertionConfiguration().getInnerOverlap();
    const double densityCell = 0.81;
    const double massCell = densityCell * cellSize[0] * cellSize[1] * cellSize[2];
    _preciceInterface = new PreciceInterface(numberCells, overLap, massCell, scenarioConfig.twoWayCoupling);
    coupling::indexing::IndexingService<dim>::getInstance().init(mdConfig, mamicoConfig, _preciceInterface, rank);
    _multiMDCellService = new coupling::services::MultiMDCellService<MY_LINKEDCELL, dim>(
        _instanceHandling->getMDSolverInterface(), _preciceInterface, mdConfig, mamicoConfig, xmlConfigurationFilename.c_str(), *_multiMDService);
    coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL, dim>::getInstance().setMacroscopicSolverInterface(_preciceInterface);
    _instanceHandling->setMacroscopicCellServices(*_multiMDCellService);
    _multiMDCellService->computeAndStoreTemperature(scenarioConfig.temp);
    _preciceAdapter = new coupling::preciceadapter::PreciceAdapter<dim>();
    _preciceAdapter->setMeshes(_preciceInterface, domainOffset, cellSize);
    _preciceAdapter->initialize();
    int cycle = 0;
#if defined(LS1_MARDYN)
    MettDeamonFeedrateDirector* mettDeamonFeedrateDirector = nullptr;
    std::list<PluginBase*>& plugins = *(global_simulation->getPluginList() );
    for (auto&& pit:plugins) {
      std::string name = pit->getPluginName();
      if(name == "MettDeamonFeedrateDirector")
        mettDeamonFeedrateDirector = dynamic_cast<MettDeamonFeedrateDirector*>(pit);
    }
    int updateFrequency = mettDeamonFeedrateDirector->getUpdateFreq();
    _preciceInterface->updateFeedrate(mettDeamonFeedrateDirector->getInitFeedrate());
#endif
    while (_preciceAdapter->isCouplingOngoing()) {
      _preciceAdapter->readData(_preciceInterface);
      _multiMDCellService->sendFromMacro2MD(_preciceAdapter->getM2mCells(), _preciceAdapter->getM2mCellIndices());
      int numberOfMDTimesteps = _preciceAdapter->getMaxTimeStepSize()/mdConfig.getSimulationConfiguration().getDt();
      _instanceHandling->simulateTimesteps(numberOfMDTimesteps, mdStepCounter, *_multiMDCellService);
      mdStepCounter+=mdConfig.getSimulationConfiguration().getNumberOfTimesteps();
#if defined(LS1_MARDYN)
      if (mdStepCounter % updateFrequency == 0) {
        double feedrate = mettDeamonFeedrateDirector->getFeedrate();
        _preciceInterface->updateFeedrate(mettDeamonFeedrateDirector->getFeedrate());
        if (rank==0) std::cout << "updating CFD feedrate given by MDFD :" << feedrate << std::endl;
      }
#endif
      _multiMDCellService->sendFromMD2Macro(_preciceAdapter->getm2MCells(), _preciceAdapter->getm2MCellIndices());
      _multiMDCellService->plotEveryMacroscopicTimestep(cycle);
      _preciceAdapter->writeData(_preciceInterface);
      _preciceAdapter->advance(numberOfMDTimesteps * mdConfig.getSimulationConfiguration().getDt());
      cycle++;
    }
  }

  void init() override { throw std::runtime_error("not supported yet"); }
  void runOneCouplingCycle(int cycle) override { throw std::runtime_error("not supported yet"); }
  coupling::solvers::AbstractCouetteSolver<3>* getSolver() override { throw std::runtime_error("not supported yet"); }

private:
  tarch::la::Vector<3, double> getCellMidPoint(const tarch::la::Vector<3, int> cellIndex, const tarch::la::Vector<3, double> domainOffset,
                                               const tarch::la::Vector<3, double> cellSize) const {
    tarch::la::Vector<3, double> cellMidPoint(domainOffset - 0.5 * cellSize);
    for (unsigned int d = 0; d < 3; d++) {
      cellMidPoint[d] = cellMidPoint[d] + cellIndex[d] * cellSize[d];
    }
    return cellMidPoint;
  }

  class PreciceInterface : public coupling::preciceadapter::PreciceInterface<3> {
  private:
    const tarch::la::Vector<3, unsigned int> _globalNumberMacroscopicCells;
    const unsigned int _overlap;
    const double _massCell;
    const std::string _M2mMeshName;
    const std::string _m2MLMeshName;
    const std::string _m2MVMeshName;
    double _feedrate;

  public:
    PreciceInterface(const tarch::la::Vector<3, int> globalNumberMacroscopicCells, const unsigned int overlap, const double massCell, const bool twoWayCoupling)
        : coupling::preciceadapter::PreciceInterface<3>(twoWayCoupling), 
        _globalNumberMacroscopicCells(globalNumberMacroscopicCells), _overlap(overlap), _massCell{massCell}, 
        _M2mMeshName("mamico-M2m-mesh"), _m2MLMeshName("mamico-m2ML-mesh"), _m2MVMeshName("mamico-m2MV-mesh") {
          coupling::preciceadapter::PreciceInterface<3>::addData(
        _M2mMeshName, coupling::preciceadapter::Data{"VelocityMacro", coupling::preciceadapter::DataType::vector});
            coupling::preciceadapter::PreciceInterface<3>::addData(
        _m2MLMeshName, coupling::preciceadapter::Data{"VelocityMicro", coupling::preciceadapter::DataType::vector});      
              coupling::preciceadapter::PreciceInterface<3>::addData(
        _m2MVMeshName, coupling::preciceadapter::Data{"VelocityMicro", coupling::preciceadapter::DataType::vector});      

        }

    bool isLiquid(tarch::la::Vector<3, unsigned int> globalCellIndex) {
      return (globalCellIndex[1] == 3 || globalCellIndex[1] == 4);
    }

    bool isVapor(tarch::la::Vector<3, unsigned int> globalCellIndex) {
      return (globalCellIndex[1] == this->_globalNumberMacroscopicCells[1]);
    }

    void updateFeedrate(double feedrate) { _feedrate = feedrate; }

    bool receiveMacroscopicQuantityFromMDSolver(tarch::la::Vector<3, unsigned int> globalCellIndex) override {
      bool isGhostCell = false;
      for (unsigned int currentDim = 0; currentDim < 3; currentDim++) {
        isGhostCell |= globalCellIndex[currentDim] > this->_globalNumberMacroscopicCells[currentDim];
        isGhostCell |= globalCellIndex[currentDim] < 1;
      }
      return !isGhostCell && (isVapor(globalCellIndex) || isLiquid(globalCellIndex));
    }

    bool sendMacroscopicQuantityToMDSolver(tarch::la::Vector<3, unsigned int> globalCellIndex) override {
      bool isGhostCell = false;
      for (unsigned int currentDim = 0; currentDim < 3; currentDim++) {
        isGhostCell |= globalCellIndex[currentDim] > this->_globalNumberMacroscopicCells[currentDim];
        isGhostCell |= globalCellIndex[currentDim] < 1;
      }
      return !isGhostCell && (globalCellIndex[1] == 1 || globalCellIndex[1] == 2 || globalCellIndex[1] == 3);
    }

    std::vector<unsigned int> getRanks(tarch::la::Vector<3, unsigned int> globalCellIndex) override { return {0}; }

    std::vector<unsigned int> getSourceRanks(tarch::la::Vector<3, unsigned int> globalCellIndex) override { 
      std::vector<unsigned int> ranks;
      int size = 1;
  #if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
      MPI_Comm_size(MPI_COMM_WORLD, &size);
  #endif
      for (int i = 0; i < size; i++)
      {
        ranks.push_back(i);
      }
      
      return ranks;
    }

    std::string getMacroscopicToMDSolverMeshName(tarch::la::Vector<3, unsigned int> globalCellIndex) override {
      return _M2mMeshName;
    }

    std::string getMDToMacroscopicSolverMeshName(tarch::la::Vector<3, unsigned int> globalCellIndex) override {
      std::string meshName;
      if (isLiquid(globalCellIndex)) {
        meshName = _m2MLMeshName;
      } else if (isVapor(globalCellIndex)) {
        meshName = _m2MVMeshName;
      } else {
        std::cout << "PreciceInterface::getMDToMacroscopicSolverMeshName: wrong cell index " << globalCellIndex << std::endl;
        exit(EXIT_FAILURE);
      }
      return meshName;
    }

    tarch::la::Vector<3, double> getMDToMacroscopicSolverMeshOffset(tarch::la::Vector<3, unsigned int> globalCellIndex) override {
      tarch::la::Vector<3, double> offset;
      if (isVapor(globalCellIndex)) {
        offset = tarch::la::Vector<3,double>{0, /*MD size*/-115/*CFD size*/-15/*half a mamico cell*/+1.25, 0};
      } else if (isLiquid(globalCellIndex)) {
        offset = tarch::la::Vector<3, double>{0.0};
      } else {
        std::cout << "PreciceInterface::getMDToMacroscopicSolverMeshOffset: wrong cell index" << std::endl;
        exit(EXIT_FAILURE);
      }
      return offset;
    }

    void readVectorData(std::string meshName, std::string dataName, coupling::datastructures::MacroscopicCell<3>* const cell, const double vx, const double vy, const double vz) override {
      tarch::la::Vector<3, double> momentum{vx, vy, vz};
      momentum=momentum*_massCell;
      cell->setMicroscopicMass(_massCell);
      cell->setMicroscopicMomentum(momentum);
    }

    void writeVectorData(std::string meshName, std::string dataName, const coupling::datastructures::MacroscopicCell<3>* const cell, double& vx, double& vy, double& vz) {
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
      } else {
        std::cout << "PreciceInterface::readData: incorrect mesh name " << meshName << " or data name " << dataName << std::endl;
        exit(EXIT_FAILURE);
      }
    }
  };

  struct ScenarioConfig : public tarch::configuration::Configuration {
    ~ScenarioConfig() {}

    void parseSubtag(tinyxml2::XMLElement* node) override {
      tinyxml2::XMLElement* subtag = node->FirstChildElement("coupling");
      tarch::configuration::ParseConfiguration::readBoolMandatory(twoWayCoupling, subtag, "two-way-coupling");
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
    bool couetteAnalytical = false;
    double channelHeight;
    double wallVelocity;
    double kinematicViscosity;
  };

  coupling::preciceadapter::PreciceAdapter<3>* _preciceAdapter;
  PreciceInterface* _preciceInterface;
  coupling::InstanceHandling<MY_LINKEDCELL, 3>* _instanceHandling;
  tarch::utils::MultiMDService<3>* _multiMDService;
  coupling::services::MultiMDCellService<MY_LINKEDCELL, 3>* _multiMDCellService;
};
