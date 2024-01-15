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
  CouetteScenario() : Scenario("couette") {}
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
    simplemd::configurations::MolecularDynamicsConfiguration mdConfig;
    tarch::configuration::ParseConfiguration::parseConfiguration<simplemd::configurations::MolecularDynamicsConfiguration>(xmlConfigurationFilename,
                                                                                                                           "molecular-dynamics", mdConfig);
    coupling::configurations::MaMiCoConfiguration<3> mamicoConfig;
    tarch::configuration::ParseConfiguration::parseConfiguration<coupling::configurations::MaMiCoConfiguration<dim>>(xmlConfigurationFilename, "mamico",
                                                                                                                     mamicoConfig);
    ScenarioConfig scenarioConfig;
    tarch::configuration::ParseConfiguration::parseConfiguration<ScenarioConfig>(xmlConfigurationFilename, "scenario", scenarioConfig);
#if defined(LS1_MARDYN)
    assert((mamicoConfig.getMacroscopicCellConfiguration().getNumberLinkedCellsPerMacroscopicCell() == tarch::la::Vector<3, unsigned int>(1)));
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
    const double densityCell = mdConfig.getDomainConfiguration().getMoleculesPerDirection()[0] *
                               mdConfig.getDomainConfiguration().getMoleculesPerDirection()[1] *
                               mdConfig.getDomainConfiguration().getMoleculesPerDirection()[2] / (domainSize[0] * domainSize[1] * domainSize[2]);
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
    double mamico_dt = mdConfig.getSimulationConfiguration().getNumberOfTimesteps() * mdConfig.getSimulationConfiguration().getDt();
    int cycle = 0;
    while (_preciceAdapter->isCouplingOngoing()) {
      if (_preciceAdapter->requiresWritingCheckpoint()) {
        std::cout << "Writing checkpoint" << std::endl;
        _instanceHandling->writeCheckpoint("checkpoint_restart_" + rank, 0);
      }
      _preciceAdapter->readData(_preciceInterface);
      _multiMDCellService->sendFromMacro2MD(_preciceAdapter->getM2mCells(), _preciceAdapter->getM2mCellIndices());
      if (!scenarioConfig.couetteAnalytical) {
        _instanceHandling->simulateTimesteps(mdConfig.getSimulationConfiguration().getNumberOfTimesteps(), mdStepCounter, *_multiMDCellService);
        mdStepCounter += mdConfig.getSimulationConfiguration().getNumberOfTimesteps();
        _multiMDCellService->sendFromMD2Macro(_preciceAdapter->getm2MCells(), _preciceAdapter->getm2MCellIndices());
        _multiMDCellService->plotEveryMacroscopicTimestep(cycle);
      } else {
        coupling::solvers::CouetteSolver<3>* couetteSolver =
            new coupling::solvers::CouetteSolver<3>(scenarioConfig.channelHeight, scenarioConfig.wallVelocity, scenarioConfig.kinematicViscosity);
        couetteSolver->advance(mamico_dt * (cycle + 1));
        for (unsigned int i = 0; i < _preciceAdapter->getm2MCells().size(); i++) {
          tarch::la::Vector<3, double> cellMidPoint =
              getCellMidPoint(coupling::indexing::convertToVector<dim>({_preciceAdapter->getm2MCellIndices()[i]}), domainOffset, cellSize);
          const tarch::la::Vector<3, double> momentum{massCell * (couetteSolver->getVelocity(cellMidPoint))};
          _preciceAdapter->getm2MCells()[i]->setMacroscopicMass(massCell);
          _preciceAdapter->getm2MCells()[i]->setMacroscopicMomentum(momentum);
        }
      }
      _preciceAdapter->writeData(_preciceInterface);
      _preciceAdapter->advance(mamico_dt);
      if (_preciceAdapter->requiresReadingCheckpoint()) {
        std::cout << "Reading checkpoint" << std::endl;
      } else {
        cycle++;
        if (_preciceAdapter->getm2MCells().size() != 0 && scenarioConfig.csvEveryTimestep >= 1 && cycle % scenarioConfig.csvEveryTimestep == 0)
          write2CSV(_preciceAdapter->getm2MCells(), _preciceAdapter->getm2MCellIndices(), cycle, numberCells, domainOffset, cellSize, overLap, rank);
      } 
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

  void write2CSV(const std::vector<coupling::datastructures::MacroscopicCell<3>*>& m2MBuffer, const unsigned int* const m2MCellIndices, const int couplingCycle,
                 const tarch::la::Vector<3, int> numberCells, const tarch::la::Vector<3, double> domainOffset, const tarch::la::Vector<3, double> cellSize,
                 const int overlap, const unsigned int rank) const {
    std::stringstream ss;
    ss << "results_" << rank << "_" << couplingCycle << ".csv";
    std::ofstream file(ss.str().c_str());
    if (!file.is_open()) {
      exit(EXIT_FAILURE);
    }
    file << "i;j;k;x;y;z;v_x;v_y;v_z;T;m" << std::endl;
    tarch::la::Vector<3, int> cellIndex;
    for (unsigned int i = 0; i < m2MBuffer.size(); i++) {
      cellIndex = coupling::indexing::convertToVector<3>({m2MCellIndices[i]});
      bool isOuter = false;
      for (unsigned int currentDim = 0; currentDim < 3; currentDim++) {
        isOuter |= (int)cellIndex[currentDim] == overlap;
        isOuter |= (int)cellIndex[currentDim] == numberCells[currentDim] + 1 - overlap;
      }
      if (!isOuter) {
        tarch::la::Vector<3, double> cellMidPoint = getCellMidPoint(cellIndex, domainOffset, cellSize);
        tarch::la::Vector<3, double> vel(m2MBuffer[i]->getMacroscopicMomentum());
        if (m2MBuffer[i]->getMacroscopicMass() != 0.0) {
          vel = (1.0 / m2MBuffer[i]->getMacroscopicMass()) * vel;
        }
        file << cellIndex[0] << ";" << cellIndex[1] << ";" << cellIndex[2] << ";" 
             << cellMidPoint[0] << ";" << cellMidPoint[1] << ";" << cellMidPoint[2] << ";" 
             << vel[0] << ";" << vel[1] << ";" << vel[2] << ";" 
             << m2MBuffer[i]->getTemperature() << ";" << m2MBuffer[i]->getMacroscopicMass();
        file << std::endl;
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
