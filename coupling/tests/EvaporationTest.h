#ifndef _COUPLING_TESTS_EVAPORATIONTEST_H_
#define _COUPLING_TESTS_EVAPORATIONTEST_H_

#include "simplemd/configurations/MolecularDynamicsConfiguration.h"
#include "tarch/configuration/ParseConfiguration.h"
#include "tarch/utils/MultiMDService.h"
#include "coupling/configurations/EvaporationConfiguration.h"
#include "coupling/configurations/MaMiCoConfiguration.h"
#include "coupling/indexing/IndexingService.h"
#include "coupling/interface/MDSimulationFactory.h"
#include "coupling/interface/impl/SimpleMD/SimpleMDSolverInterface.h"
#include "coupling/services/MultiMDCellService.h"
#include "coupling/solvers/PreciceAdapter2.h"
#include "coupling/solvers/CouetteSolverInterface.h"
#include "coupling/CouplingMDDefinitions.h"
#include "coupling/tests/Test.h"
#include <mpi.h>
#include <random>
#include <sys/time.h>

class EvaporationTest : public Test {
public:
  EvaporationTest() : Test("EvaporationTest") {}
  virtual ~EvaporationTest() {}

  virtual void run() {
    init();
    for (int cycle = 0; cycle < _evaporationConfig.couplingCycles; cycle++) {
      advanceMacro(cycle);
      advanceMicro(cycle);
    }
    shutdown();
  }

private:
  void init() {
    const unsigned int dim = 3;
    // init MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
    // parse configurations
    std::string xmlConfigurationFilename("evaporation.xml");
    tarch::configuration::ParseConfiguration::parseConfiguration<simplemd::configurations::MolecularDynamicsConfiguration>(xmlConfigurationFilename, "molecular-dynamics", _mdConfig);
    if (!_mdConfig.isValid()) {
      std::cout << "ERROR EvaporationTest: Invalid molecular dynamics config." << std::endl;
      exit(EXIT_FAILURE);
    }
    tarch::configuration::ParseConfiguration::parseConfiguration<coupling::configurations::MaMiCoConfiguration<dim>>(xmlConfigurationFilename, "mamico", _mamicoConfig);
    if (!_mamicoConfig.isValid()) {
      std::cout << "ERROR EvaporationTest: Invalid MaMiCo config." << std::endl;
      exit(EXIT_FAILURE);
    }
    _evaporationConfig = coupling::configurations::EvaporationConfig::parseConfiguration(xmlConfigurationFilename);
    // init solvers
    _macroSolver = new coupling::solvers::PreciceAdapter<dim>();
    if (_macroSolver != NULL) std::cout << "Macro solver not null on rank: " << _rank << std::endl;
    _multiMDService = new tarch::utils::MultiMDService<dim>(_mdConfig.getMPIConfiguration().getNumberOfProcesses(), _evaporationConfig.totalNumberMDSimulations);
    _localMDInstances = _multiMDService->getLocalNumberOfMDSimulations();
    for (unsigned int i = 0; i < _localMDInstances; i++) {
      _mdSolvers.push_back(coupling::interface::SimulationAndInterfaceFactory::getInstance().getMDSimulation(_mdConfig, _mamicoConfig, _multiMDService->getLocalCommunicator()));
      if (_mdSolvers[i] == NULL) {
        std::cout << "ERROR EvaporationTest: _mdSolvers[" << i << "]==NULL!" << std::endl;
        exit(EXIT_FAILURE);
      }
    }
    if (_rank == 0) {
      gettimeofday(&_tv.start, NULL);
    }
    for (unsigned int i = 0; i < _localMDInstances; i++) {
      _mdSolvers[i]->init(*_multiMDService, _multiMDService->getGlobalNumberOfLocalMDSimulation(i));
      _mdSolversInterface.push_back(coupling::interface::SimulationAndInterfaceFactory::getInstance().getMDSolverInterface(_mdConfig, _mamicoConfig, _mdSolvers[i]));
      if (_mdSolversInterface[i] == NULL) {
        std::cout << "ERROR EvaporationTest: mdSolverInterface[" << i << "] == NULL!" << std::endl;
        exit(EXIT_FAILURE);
      }
    }
    const tarch::la::Vector<dim, double> mdOffset{_mdConfig.getDomainConfiguration().getGlobalDomainOffset()};
    const tarch::la::Vector<dim, double> mamicoMeshsize{_mamicoConfig.getMacroscopicCellConfiguration().getMacroscopicCellSize()};
    const tarch::la::Vector<dim, unsigned int> globalNumberMacroscopicCells{getGlobalNumberMacroscopicCells()};
    const unsigned int outerRegion{_mamicoConfig.getMomentumInsertionConfiguration().getInnerOverlap()};
    tarch::la::Vector<dim, unsigned int> offsetMDDomain(0);
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
    coupling::interface::MacroscopicSolverInterface<dim>* macroSolverInterface = new coupling::solvers::CouetteSolverInterface<3>(globalNumberMacroscopicCells, outerRegion);
    if (macroSolverInterface == NULL) {
      std::cout << "ERROR EvaporationTest: rank=" << _rank << ": macroSolverInterface==NULL!" << std::endl;
      exit(EXIT_FAILURE);
    }
    _multiMDCellService = new coupling::services::MultiMDCellService<MY_LINKEDCELL, dim>(_mdSolversInterface, macroSolverInterface, _mdConfig, _mamicoConfig, xmlConfigurationFilename.c_str(), *_multiMDService);
    coupling::indexing::IndexingService<dim>::getInstance().init(_mdConfig, _mamicoConfig, macroSolverInterface,  _rank);

    /*
    // TO DELETE
    const tarch::la::Vector<dim, double> macroscopicCellSize(_mamicoConfig.getMacroscopicCellConfiguration().getMacroscopicCellSize());
    using namespace coupling;
    using namespace indexing;
    for (CellIndex<dim, IndexTrait::vector> cellIndex : coupling::indexing::CellIndex<dim, IndexTrait::vector>()) {
      tarch::la::Vector<dim, int> cellVectorIndex = cellIndex.get();
      unsigned int cellScalarIndex = convertToScalar<dim, IndexTrait::vector>(cellIndex);
      if (cellScalarIndex < 70) {
        std::cout << "cellScalarIndex:" << cellScalarIndex << std::endl;
        std::cout << "cellVectorIndex:" << cellVectorIndex << std::endl;
        for (unsigned int currentDim = 0; currentDim < dim; currentDim++) {
          std::cout << "dim:" << currentDim << "," << cellVectorIndex[currentDim]*macroscopicCellSize[currentDim] + 0.5*macroscopicCellSize[currentDim] << std::endl;
        }
      }
    }
    // END TO DELETE
    */

    coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL, dim>::getInstance().setMacroscopicSolverInterface(macroSolverInterface);
    for (unsigned int i = 0; i < _localMDInstances; i++) {
      _mdSolvers[i]->setMacroscopicCellService(&(_multiMDCellService->getMacroscopicCellService(i)));
      _multiMDCellService->getMacroscopicCellService(i).computeAndStoreTemperature(_evaporationConfig.temp);
    }
    allocateCFDToMDBuffer(_multiMDCellService->getMacroscopicCellService(0).getIndexConversion(),*macroSolverInterface);
    allocateMDToCFDBuffer(_multiMDCellService->getMacroscopicCellService(0).getIndexConversion(),*macroSolverInterface);
    _macroSolver->setCouplingMesh();
    if (_rank == 0) {
      gettimeofday(&_tv.end, NULL);
      double runtime = (_tv.end.tv_sec - _tv.start.tv_sec) * 1000000 + (_tv.end.tv_usec - _tv.start.tv_usec);
      std::cout << "Initialization: " << (int)(runtime / 1000) << "ms" << std::endl;
      gettimeofday(&_tv.start_total, NULL);
    }
    std::cout << "Finish EvaporationTest::initSolvers() " << std::endl;
  }

  void advanceMacro(int cycle) {
    if (_macroSolver != NULL) {
      if (_rank == 0) {
        gettimeofday(&_tv.start, NULL);
      }
      _macroSolver->advance(_mdConfig.getSimulationConfiguration().getDt() * _mdConfig.getSimulationConfiguration().getNumberOfTimesteps());
      if (_rank == 0) {
        gettimeofday(&_tv.end, NULL);
        _tv.macro += (_tv.end.tv_sec - _tv.start.tv_sec) * 1000000 + (_tv.end.tv_usec - _tv.start.tv_usec);
      }
      fillCFDToMDBuffer();
    }
    _multiMDCellService->sendFromMacro2MD(_buf.CFDToMDBuffer, _buf.globalCellIndices4CFDToMDBuffer);
  }

  void advanceMicro(int cycle) {
    if (_rank == 0) {
      gettimeofday(&_tv.start, NULL);
    }
    for (unsigned int i = 0; i < _localMDInstances; i++) {
      coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL, 3>::getInstance().setMacroscopicCellService(&(_multiMDCellService->getMacroscopicCellService(i)));
      coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL, 3>::getInstance().setMDSolverInterface(_mdSolversInterface[i]);
      _mdSolvers[i]->simulateTimesteps(_mdConfig.getSimulationConfiguration().getNumberOfTimesteps(), _mdStepCounter);
      _multiMDCellService->getMacroscopicCellService(i).plotEveryMacroscopicTimestep(cycle);
    }
    _mdStepCounter += _mdConfig.getSimulationConfiguration().getNumberOfTimesteps();
    if (_rank == 0) {
      gettimeofday(&_tv.end, NULL);
      _tv.micro += (_tv.end.tv_sec - _tv.start.tv_sec) * 1000000 + (_tv.end.tv_usec - _tv.start.tv_usec);
    }
  }

  void shutdown() {
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
      // the shutdown operation may also delete the md solver interface; therefore, we update the MD solver interface in the vector _mdSolversInteface after the
      // shutdown is completed
      coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL, 3>::getInstance().setMDSolverInterface(_mdSolversInterface[i]);
      _mdSolvers[i]->shutdown();
      delete _mdSolvers[i];
      _mdSolvers[i] = NULL;
      _mdSolversInterface[i] = coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL, 3>::getInstance().getMDSolverInterface();
    }
    _mdSolvers.clear();
    for (unsigned int i = 0; i < _localMDInstances; i++) {
      if (_mdSolversInterface[i] != NULL) {
        delete _mdSolversInterface[i];
        _mdSolversInterface[i] = NULL;
      }
    }
    _mdSolversInterface.clear();

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

  tarch::la::Vector<3, unsigned int> getGlobalNumberMacroscopicCells() const {
    return tarch::la::Vector<3, unsigned int>{1,1,5};
  }

  void allocateCFDToMDBuffer(const coupling::IndexConversion<3>& indexConversion, coupling::interface::MacroscopicSolverInterface<3>& macroSolverInterface) {
    const unsigned int num = 3; // just 3 cells shall be sent in oneD
    deleteBuffer(_buf.CFDToMDBuffer);
    // TODO: Adapt for parallel setup, check where this is necessary
    // allocate array for cell indices
    unsigned int* indices = new unsigned int[num]{31, 40, 49}; // refers to third, fourth, fifth cell
    if (indices == NULL) {
      std::cout << "ERROR EvaporationTest::allocateCFDToMDBuffer(): indices==NULL!" << std::endl;
      exit(EXIT_FAILURE);
    }
    for (unsigned int i = 0; i < num; i++) {
      _buf.CFDToMDBuffer.push_back(new coupling::datastructures::MacroscopicCell<3>());
      if (_buf.CFDToMDBuffer[_buf.CFDToMDBuffer.size() - 1] == NULL) {
        std::cout << "ERROR EvaporationTest::allocateCFDToMDBuffer: CFDToMDBuffer[" << _buf.CFDToMDBuffer.size() - 1 << "]==NULL!" << std::endl;
        exit(EXIT_FAILURE);
      }
    }
    _buf.globalCellIndices4CFDToMDBuffer = indices;
  }


  void allocateMDToCFDBuffer(const coupling::IndexConversion<3>& indexConversion, coupling::interface::MacroscopicSolverInterface<3>& macroSolverInterface) {
    const unsigned int num = 2;
    deleteBuffer(_buf.MDToCFDBuffer);
    unsigned int* indices = new unsigned int[num]{58,67};
    if (indices == NULL) {
      std::cout << "ERROR EvaporationTest::allocateMDToCFDBuffer(): indices==NULL!" << std::endl;
      exit(EXIT_FAILURE);
    }
    for (unsigned int i = 0; i < num; i++) {
      _buf.MDToCFDBuffer.push_back(new coupling::datastructures::MacroscopicCell<3>());
      if (_buf.MDToCFDBuffer[_buf.MDToCFDBuffer.size() - 1] == NULL) {
        std::cout << "ERROR EvaporationTest::allocateMDToCFDBuffer: MDToCFDBuffer[" << _buf.MDToCFDBuffer.size() - 1 << "]==NULL!" << std::endl;
        exit(EXIT_FAILURE);
      }
    }
    _buf.globalCellIndices4MDToCFDBuffer = indices;
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

  void fillCFDToMDBuffer() const {
    std::cout << "fillCFDToMDBuffer" << std::endl;
    const coupling::IndexConversion<3>& indexConversion{_multiMDCellService->getIndexConversion()};
    const unsigned int size = _buf.CFDToMDBuffer.size();
    const tarch::la::Vector<3, double> domainOffset(indexConversion.getGlobalMDDomainOffset());
    const tarch::la::Vector<3, double> macroscopicCellSize(indexConversion.getMacroscopicCellSize());
    for (unsigned int i = 0; i < size; i++) {
      const tarch::la::Vector<3, unsigned int> globalIndex(indexConversion.getGlobalVectorCellIndex(_buf.globalCellIndices4CFDToMDBuffer[i]));
      tarch::la::Vector<3, double> cellMidPoint(domainOffset - 0.5 * macroscopicCellSize);
      for (unsigned int d = 0; d < 3; d++) {
        cellMidPoint[d] = cellMidPoint[d] + ((double)globalIndex[d]) * macroscopicCellSize[d];
      }
      std::cout << "cellMidPoint:" << cellMidPoint << std::endl;
      double mass = macroscopicCellSize[0] * macroscopicCellSize[1] * macroscopicCellSize[2];
      mass *= _macroSolver->getDensity(cellMidPoint);
      tarch::la::Vector<3, double> momentum(mass * _macroSolver->getVelocity(cellMidPoint));
      _buf.CFDToMDBuffer[i]->setMicroscopicMass(mass);
      _buf.CFDToMDBuffer[i]->setMicroscopicMomentum(momentum);
    }
  }

  struct CouplingBuffer {
    std::vector<coupling::datastructures::MacroscopicCell<3>*> CFDToMDBuffer;
    unsigned int* globalCellIndices4CFDToMDBuffer;
    std::vector<coupling::datastructures::MacroscopicCell<3>*> MDToCFDBuffer;
    unsigned int* globalCellIndices4MDToCFDBuffer;
  };

  struct TimingValues {
    timeval start_total;
    timeval start;
    timeval end;
    double micro{0};
    double macro{0};
  };


  int _rank{0};
  simplemd::configurations::MolecularDynamicsConfiguration _mdConfig;
  coupling::configurations::MaMiCoConfiguration<3> _mamicoConfig;
  coupling::configurations::EvaporationConfig _evaporationConfig;
  unsigned int _mdStepCounter{0};
  coupling::solvers::PreciceAdapter<3>* _macroSolver;
  tarch::utils::MultiMDService<3>* _multiMDService;
  coupling::services::MultiMDCellService<MY_LINKEDCELL, 3>* _multiMDCellService;
  CouplingBuffer _buf;
  unsigned int _localMDInstances;
  std::vector<coupling::interface::MDSolverInterface<MY_LINKEDCELL, 3>*> _mdSolversInterface;
  std::vector<coupling::interface::MDSimulation*> _mdSolvers;
  TimingValues _tv;
};

#endif // _COUPLING_TESTS_EVAPORATIONTEST_H_
