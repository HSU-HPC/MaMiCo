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
#include "coupling/solvers/PreciceAdapter.h"
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

class PreciceScenario : public Scenario {
public:
  PreciceScenario() : Scenario("PreciceScenario") {}
  ~PreciceScenario() {
    deleteBuffer(_buf._M2mBuffer);
    if (_buf._M2mCellIndices != NULL) {
      delete[] _buf._M2mCellIndices;
      _buf._M2mCellIndices = NULL;
    }
    deleteBuffer(_buf._m2MBuffer);
    if (_buf._m2MCellIndices != NULL) {
      delete[] _buf._m2MCellIndices;
      _buf._m2MCellIndices = NULL;
    }
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
    global_log = new Log::Logger(Log::Error); // Info
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

    coupling::IndexConversion<3>* indexConversion = initIndexConversion(mamicoConfig.getMacroscopicCellConfiguration().getMacroscopicCellSize(),
                                             _multiMDService->getNumberProcessesPerMDSimulation(), _multiMDService->getGlobalRank(),
                                             mdConfig.getDomainConfiguration().getGlobalDomainSize(),
                                             mdConfig.getDomainConfiguration().getGlobalDomainOffset(),
                                             mamicoConfig.getParallelTopologyConfiguration().getParallelTopologyType(), computeTopologyOffset());
    _macroscopicSolverInterface = new coupling::solvers::PreciceInterface<dim>(numberCells, overLap, rank, indexConversion);

    coupling::indexing::IndexingService<dim>::getInstance().init(mdConfig, mamicoConfig, _macroscopicSolverInterface, rank);

    _multiMDCellService = new coupling::services::MultiMDCellService<MY_LINKEDCELL, dim>(
        _instanceHandling->getMDSolverInterface(), _macroscopicSolverInterface, mdConfig, mamicoConfig, xmlConfigurationFilename.c_str(), *_multiMDService);
    
    coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL, dim>::getInstance().setMacroscopicSolverInterface(_macroscopicSolverInterface);
    
    _instanceHandling->setMacroscopicCellServices(*_multiMDCellService);
    
    _multiMDCellService->computeAndStoreTemperature(scenarioConfig.temp);
    allocateM2mBuffer(rank);
    allocatem2MBuffer(rank);
    
    _preciceAdapter = new coupling::solvers::PreciceAdapter<dim>("mamico-M2m-mesh", "mamico-m2M-mesh", "VelocityMacro", "VelocityMicro");
    _preciceAdapter->setMeshes(_macroscopicSolverInterface, domainOffset, cellSize);
    double precice_dt = _preciceAdapter->initialize();
    double mamico_dt = mdConfig.getSimulationConfiguration().getNumberOfTimesteps() * mdConfig.getSimulationConfiguration().getDt();
    int cycle = 0;
    const double densityCell = mdConfig.getDomainConfiguration().getMoleculesPerDirection()[0] *
                               mdConfig.getDomainConfiguration().getMoleculesPerDirection()[1] *
                               mdConfig.getDomainConfiguration().getMoleculesPerDirection()[2] / (domainSize[0] * domainSize[1] * domainSize[2]);
    const double massCell = densityCell * cellSize[0] * cellSize[1] * cellSize[2];
    std::cout << "1" << std::endl;
    while (_preciceAdapter->isCouplingOngoing()) {
      _preciceAdapter->readData();
      std::cout << "2" << std::endl;
      for (unsigned int i = 0; i < _buf._M2mBuffer.size(); i++) {
        tarch::la::Vector<3, double> cellMidPoint =
            getCellMidPoint(coupling::indexing::convertToVector<dim>({_buf._M2mCellIndices[i]}), domainOffset, cellSize);
        tarch::la::Vector<3, double> momentum(massCell * _preciceAdapter->getVelocity(cellMidPoint));
        _buf._M2mBuffer[i]->setMicroscopicMass(massCell);
        _buf._M2mBuffer[i]->setMicroscopicMomentum(momentum);
      }
      std::cout << "3" << std::endl;
      _multiMDCellService->sendFromMacro2MD(_buf._M2mBuffer, _buf._M2mCellIndices);
      if (!scenarioConfig.couetteAnalytical) {
        std::cout << "4" << std::endl;
        _instanceHandling->simulateTimesteps(mdConfig.getSimulationConfiguration().getNumberOfTimesteps(), mdStepCounter, *_multiMDCellService);
        mdStepCounter += mdConfig.getSimulationConfiguration().getNumberOfTimesteps();
        _multiMDCellService->sendFromMD2Macro(_buf._m2MBuffer, _buf._m2MCellIndices);
        _multiMDCellService->plotEveryMacroscopicTimestep(cycle);
        std::cout << "5" << std::endl;
      } else {
        coupling::solvers::CouetteSolver<3>* couetteSolver =
            new coupling::solvers::CouetteSolver<3>(scenarioConfig.channelHeight, scenarioConfig.wallVelocity, scenarioConfig.dynamicViscosity);
        couetteSolver->advance(mamico_dt * (cycle + 1));
        for (unsigned int i = 0; i < _buf._m2MBuffer.size(); i++) {
          tarch::la::Vector<3, double> cellMidPoint =
              getCellMidPoint(coupling::indexing::convertToVector<dim>({_buf._m2MCellIndices[i]}), domainOffset, cellSize);
          const tarch::la::Vector<3, double> momentum{massCell * (couetteSolver->getVelocity(cellMidPoint))};
          _buf._m2MBuffer[i]->setMacroscopicMass(massCell);
          _buf._m2MBuffer[i]->setMacroscopicMomentum(momentum);
        }
      }
      _preciceAdapter->setMDBoundaryValues(_buf._m2MBuffer, _buf._m2MCellIndices);
      _preciceAdapter->writeData();
      std::cout << "6" << std::endl;
      precice_dt = _preciceAdapter->advance(mamico_dt);
      std::cout << "7" << std::endl;
      if (precice_dt < mamico_dt)
        std::cout << "warning: rank " << rank << " maximum timestep from preCICE (" << precice_dt << ") is lower than MaMiCo timestep (" << mamico_dt << ")"
                  << std::endl;
      cycle++;
      if (_buf._m2MBuffer.size() != 0 && scenarioConfig.csvEveryTimestep >= 1 && cycle % scenarioConfig.csvEveryTimestep == 0)
        write2CSV(_buf._m2MBuffer, _buf._m2MCellIndices, cycle, numberCells, domainOffset, cellSize, overLap, rank);
    }
  }

private:
  unsigned int computeTopologyOffset() const {
    // determine topology offset of this rank
    const unsigned int intNumberProcesses = computeScalarNumberProcesses();
    const unsigned int topologyOffset = (_multiMDService->getGlobalRank() / intNumberProcesses) * intNumberProcesses;
    return topologyOffset;
  }

  unsigned int computeScalarNumberProcesses() const {
    unsigned int np = _multiMDService->getNumberProcessesPerMDSimulation()[0];
    for (unsigned int d = 1; d < 3; d++) {
      np = np * _multiMDService->getNumberProcessesPerMDSimulation()[d];
    }
    return np;
  }

  coupling::IndexConversion<3>* initIndexConversion(tarch::la::Vector<3, double> macroscopicCellSize, tarch::la::Vector<3, unsigned int> numberProcesses,
                                                      unsigned int rank, tarch::la::Vector<3, double> globalMDDomainSize,
                                                      tarch::la::Vector<3, double> globalMDDomainOffset,
                                                      coupling::paralleltopology::ParallelTopologyType parallelTopologyType,
                                                      unsigned int topologyOffset) const {

    tarch::la::Vector<3, unsigned int> globalNumberMacroscopicCells(0);
    for (unsigned int d = 0; d < 3; d++) {
      globalNumberMacroscopicCells[d] = (unsigned int)floor(globalMDDomainSize[d] / macroscopicCellSize[d] + 0.5);
      if (fabs(globalNumberMacroscopicCells[d] * macroscopicCellSize[d] - globalMDDomainSize[d]) > 1e-13) {
        std::cout << "coupling::services::MultiMDCellService::initIndexConversi"
                     "on(): Deviation of domain size > 1e-13!"
                  << std::endl;
      }
    }
    coupling::IndexConversion<3>* ic = new coupling::IndexConversion<3>(globalNumberMacroscopicCells, numberProcesses, rank, globalMDDomainSize,
                                                                            globalMDDomainOffset, parallelTopologyType, topologyOffset);
    if (ic == NULL) {
      std::cout << "coupling::services::MultiMDCellService::initIndexConversion"
                   "(): ic==NULL!"
                << std::endl;
      exit(EXIT_FAILURE);
    }
    return ic;
  }

  tarch::la::Vector<3, double> getCellMidPoint(const tarch::la::Vector<3, int> cellIndex, const tarch::la::Vector<3, double> domainOffset,
                                               const tarch::la::Vector<3, double> cellSize) const {
    tarch::la::Vector<3, double> cellMidPoint(domainOffset - 0.5 * cellSize);
    for (unsigned int d = 0; d < 3; d++) {
      cellMidPoint[d] = cellMidPoint[d] + cellIndex[d] * cellSize[d];
    }
    return cellMidPoint;
  }

  void allocateM2mBuffer(const unsigned int rank) {
    using namespace coupling::indexing;
    const unsigned int dim = 3;
    std::vector<unsigned int> M2mBufferCellIndices;
    for (CellIndex<dim, IndexTrait::vector> cellIndex : CellIndex<dim, IndexTrait::vector>()) {
      tarch::la::Vector<3, unsigned int> cellVectorIndex = static_cast<tarch::la::Vector<dim, unsigned int>>(cellIndex.get());
      if (_macroscopicSolverInterface->sendMacroscopicQuantityToMDSolver(cellVectorIndex)) {
        std::vector<unsigned int> ranks = _macroscopicSolverInterface->getSourceRanks(cellVectorIndex);
        bool containsThisRank = false;
        for (unsigned int k = 0; k < ranks.size(); k++) {
          containsThisRank = containsThisRank || (ranks[k] == rank);
        }
        if (containsThisRank) {
          _buf._M2mBuffer.push_back(new coupling::datastructures::MacroscopicCell<3>());
          M2mBufferCellIndices.push_back(convertToScalar<dim>(cellIndex));
        }
      }
    }
    std::cout << "@@@@ _buf._M2mBuffer.size():" << _buf._M2mBuffer.size() << std::endl;
    _buf._M2mCellIndices = new unsigned int[_buf._M2mBuffer.size()];
    std::copy(M2mBufferCellIndices.begin(), M2mBufferCellIndices.end(), _buf._M2mCellIndices);
  }

  void allocatem2MBuffer(const unsigned int rank) {
    using namespace coupling::indexing;
    const unsigned int dim = 3;
    std::vector<unsigned int> m2MBufferCellIndices;
    for (CellIndex<dim, IndexTrait::vector> cellIndex : CellIndex<dim, IndexTrait::vector>()) {
      tarch::la::Vector<3, unsigned int> cellVectorIndex = static_cast<tarch::la::Vector<dim, unsigned int>>(cellIndex.get());
      if (_macroscopicSolverInterface->receiveMacroscopicQuantityFromMDSolver(cellVectorIndex)) {
        std::vector<unsigned int> ranks = _macroscopicSolverInterface->getTargetRanks(cellVectorIndex);
        bool containsThisRank = false;
        for (unsigned int k = 0; k < ranks.size(); k++) {
          containsThisRank = containsThisRank || (ranks[k] == (unsigned int)rank);
        }
        if (containsThisRank) {
          _buf._m2MBuffer.push_back(new coupling::datastructures::MacroscopicCell<3>());
          m2MBufferCellIndices.push_back(convertToScalar<dim>(cellIndex));
        }
      }
    }
    std::cout << "@@@@ _buf._m2MBuffer.size():" << _buf._m2MBuffer.size() << std::endl;
    _buf._m2MCellIndices = new unsigned int[_buf._m2MBuffer.size()];
    std::copy(m2MBufferCellIndices.begin(), m2MBufferCellIndices.end(), _buf._m2MCellIndices);
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
      tarch::configuration::ParseConfiguration::readDoubleOptional(dynamicViscosity, subtag, "dynamic-viscosity");
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
    double dynamicViscosity;
  };

  struct CouplingBuffer {
    std::vector<coupling::datastructures::MacroscopicCell<3>*> _M2mBuffer;
    unsigned int* _M2mCellIndices;
    std::vector<coupling::datastructures::MacroscopicCell<3>*> _m2MBuffer;
    unsigned int* _m2MCellIndices;
  };

  coupling::solvers::PreciceAdapter<3>* _preciceAdapter;
  coupling::solvers::PreciceInterface<3>* _macroscopicSolverInterface;
  coupling::InstanceHandling<MY_LINKEDCELL, 3>* _instanceHandling;
  tarch::utils::MultiMDService<3>* _multiMDService;
  coupling::services::MultiMDCellService<MY_LINKEDCELL, 3>* _multiMDCellService;
  CouplingBuffer _buf;
};
