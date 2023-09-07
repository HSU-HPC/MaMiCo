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

namespace precice_scenario{
  class CouetteScenario;
}


class precice_scenario::CouetteScenario : public Scenario {
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
    if (_macroscopicSolverInterface != NULL) {
      delete _macroscopicSolverInterface;
      _macroscopicSolverInterface = NULL;
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

  void run() {
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

    _macroscopicSolverInterface = new precice_scenario::PreciceInterface<dim>(numberCells, overLap);

    coupling::indexing::IndexingService<dim>::getInstance().init(mdConfig, mamicoConfig, _macroscopicSolverInterface, rank);

    _multiMDCellService = new coupling::services::MultiMDCellService<MY_LINKEDCELL, dim>(
        _instanceHandling->getMDSolverInterface(), _macroscopicSolverInterface, mdConfig, mamicoConfig, xmlConfigurationFilename.c_str(), *_multiMDService);
    
    coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL, dim>::getInstance().setMacroscopicSolverInterface(_macroscopicSolverInterface);
    
    _instanceHandling->setMacroscopicCellServices(*_multiMDCellService);
    
    _multiMDCellService->computeAndStoreTemperature(scenarioConfig.temp);
    
    _preciceAdapter = new precice_scenario::PreciceAdapter<dim>();
    _preciceAdapter->setMeshes(_macroscopicSolverInterface, domainOffset, cellSize);
    double precice_dt = _preciceAdapter->initialize();
    double mamico_dt = mdConfig.getSimulationConfiguration().getNumberOfTimesteps() * mdConfig.getSimulationConfiguration().getDt();
    int cycle = 0;
    const double densityCell = mdConfig.getDomainConfiguration().getMoleculesPerDirection()[0] *
                               mdConfig.getDomainConfiguration().getMoleculesPerDirection()[1] *
                               mdConfig.getDomainConfiguration().getMoleculesPerDirection()[2] / (domainSize[0] * domainSize[1] * domainSize[2]);
    const double massCell = densityCell * cellSize[0] * cellSize[1] * cellSize[2];
    while (_preciceAdapter->isCouplingOngoing()) {
      _preciceAdapter->readData(massCell);
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
      _preciceAdapter->writeData();
      precice_dt = _preciceAdapter->advance(mamico_dt);
      if (precice_dt < mamico_dt)
        std::cout << "warning: rank " << rank << " maximum timestep from preCICE (" << precice_dt << ") is lower than MaMiCo timestep (" << mamico_dt << ")"
                  << std::endl;
      cycle++;
      if (_preciceAdapter->getm2MCells().size() != 0 && scenarioConfig.csvEveryTimestep >= 1 && cycle % scenarioConfig.csvEveryTimestep == 0)
        write2CSV(_preciceAdapter->getm2MCells(), _preciceAdapter->getm2MCellIndices(), cycle, numberCells, domainOffset, cellSize, overLap, rank);
    }
  }

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

  void deleteBuffer(std::vector<coupling::datastructures::MacroscopicCell<3>*>& buffer) const {
    for (unsigned int i = 0; i < buffer.size(); i++) {
      if (buffer[i] != NULL) {
        delete buffer[i];
        buffer[i] = NULL;
      }
    }
    buffer.clear();
  }

  struct ScenarioConfig : public tarch::configuration::Configuration {
    ~ScenarioConfig() {}

    void parseSubtag(tinyxml2::XMLElement* node) override {
      tinyxml2::XMLElement* subtag = node->FirstChildElement("coupling");
      tarch::configuration::ParseConfiguration::readIntMandatory(csvEveryTimestep, subtag, "write-csv-every-timestep");
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
  };

  precice_scenario::PreciceAdapter<3>* _preciceAdapter;
  precice_scenario::PreciceInterface<3>* _macroscopicSolverInterface;
  coupling::InstanceHandling<MY_LINKEDCELL, 3>* _instanceHandling;
  tarch::utils::MultiMDService<3>* _multiMDService;
  coupling::services::MultiMDCellService<MY_LINKEDCELL, 3>* _multiMDCellService;
};
