// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _COUPLING_TESTS_EVAPORATIONTEST_H_
#define _COUPLING_TESTS_EVAPORATIONTEST_H_

#include "tarch/utils/MultiMDService.h"
#include "simplemd/configurations/MolecularDynamicsConfiguration.h"
#include "coupling/tests/Test.h"
#include "coupling/CouplingMDDefinitions.h"
#include "tarch/configuration/ParseConfiguration.h"
#include "coupling/solvers/CouetteSolver.h"
#include "coupling/solvers/LBCouetteSolver.h"
#if(BUILD_WITH_OPENFOAM)
#include "coupling/solvers/FoamClass.h"
#include "coupling/solvers/FoamSolverInterface.h"
#endif
#include "coupling/solvers/FDCouetteSolver.h"
#include "coupling/solvers/CouetteSolverInterface.h"
#include "coupling/solvers/LBCouetteSolverInterface.h"
#include "coupling/interface/MDSimulationFactory.h"
#include "coupling/interface/impl/SimpleMD/SimpleMDSolverInterface.h"
#include "coupling/configurations/MaMiCoConfiguration.h"
#include "coupling/services/MultiMDCellService.h"
#if (COUPLING_MD_PARALLEL==COUPLING_MD_YES)
#include <mpi.h>
#endif
#include <sys/time.h>
#include <random>

/**
 * @brief Runs a evaporation scenario simulation
 * @author Helene Wittenberg
 */
class EvaporationTest: public Test {
public:
  /** @brief simple constructor */
  EvaporationTest(): Test("EvaporationTest"), _generator(0){}
  /** @brief a dummy destructor */
  virtual ~EvaporationTest(){}

  /** triggers void init(), runOneCouplingCycle() and shutdown()
   *  @brief runs the simulation */
  virtual void run(){
    init();
    for (int cycle = 0; cycle < _cfg.couplingCycles; cycle++) {
      runOneCouplingCycle(cycle);
	 }
    shutdown();
  }

private:
  /** Defines the type of continuum solver for the coupled simulation  */
  enum MacroSolverType{COUETTE_LB=0, ///< the Lattice-Bolzmann solver is used (coupling::solvers::LBCouetteSolver)
  };

  /** triggers initMPI(), parseConfiguration(), and initSolvers()
   *  @brief initialises everthing necessary for the test */
  void init(){
    initMPI();
    parseConfigurations();
    initSolvers();
  }

  /** it advances the macro (advanceMacro()) and micro solver (advanceMicro),
   *  computes the signal to noise ratio (computeSNR()) and sends the data from the macro to
   *  the micro solver (twoWayCoupling())
   *  @brief combines the functioniality necessary for a cycle of the coupled simulation  */
  void runOneCouplingCycle(int cycle){
    advanceMacro(cycle);
    advanceMicro(cycle);
  }

  /** @brief initialises all MPI variables  */
  void initMPI(){
    _rank = 0;
    #if (COUPLING_MD_PARALLEL==COUPLING_MD_YES)
    MPI_Comm_rank(MPI_COMM_WORLD,&_rank);
    #endif
  }

  /** @brief reads the configuration from the xml file and calls parseCouetteTestConfiguration() */
  void parseConfigurations(){
    tarch::configuration::ParseConfiguration::parseConfiguration<simplemd::configurations::MolecularDynamicsConfiguration>("couette.xml","molecular-dynamics",_simpleMDConfig);
    if (!_simpleMDConfig.isValid()){std::cout << "ERROR CouetteTest: Invalid SimpleMD config!" << std::endl; exit(EXIT_FAILURE);}
    tarch::configuration::ParseConfiguration::parseConfiguration<coupling::configurations::MaMiCoConfiguration<3> >("couette.xml","mamico",_mamicoConfig);
    if (!_mamicoConfig.isValid()){ std::cout << "ERROR CouetteTest: Invalid MaMiCo config!" << std::endl; exit(EXIT_FAILURE); }

    parseEvaporationTestConfiguration();
  }

  /** @brief reads the configuartion, checks that mandatory data is provided and stores the data in variables */
  void parseEvaporationTestConfiguration(){
    tinyxml2::XMLDocument conffile;
    tinyxml2::XMLElement *node = NULL;
    conffile.LoadFile("couette.xml");
    node = conffile.FirstChildElement("couette-test");
    if (node == NULL){
      std::cout << "Could not read input file couette.xml: missing element <couette-test>" << std::endl;
      exit(EXIT_FAILURE);
    }
    tinyxml2::XMLElement *n_mamico = node->NextSiblingElement();
    if(n_mamico == NULL){
      std::cout << "Could not read input file couette.xml: missing element <mamico>" << std::endl;
      exit(EXIT_FAILURE);
    }
	tinyxml2::XMLElement *n_md = n_mamico->NextSiblingElement();
    if(n_md == NULL){
      std::cout << "Could not read input file couette.xml: missing element <molecular-dynamics>" << std::endl;
      exit(EXIT_FAILURE);
    }
	tinyxml2::XMLElement *n_fp = n_md->NextSiblingElement();
    if(n_fp == NULL){
      std::cout << "Could not read input file couette.xml: missing element <filter-pipeline>" << std::endl;
      exit(EXIT_FAILURE);
    }
	tinyxml2::XMLElement *n_unexpected = n_fp->NextSiblingElement();
    if(n_unexpected != NULL){
      std::cout << "Could not read input file couette.xml: unknown element " << n_unexpected->Name() << std::endl;
      exit(EXIT_FAILURE);
    }
    tinyxml2::XMLElement* subtag = node->FirstChildElement("domain");
    if (subtag == NULL){
      std::cout << "Could not read input file couette.xml: Missing subtag: domain" << std::endl;
      exit(EXIT_FAILURE);
    }
    tarch::configuration::ParseConfiguration::readDoubleMandatory(_cfg.channelheight,subtag,"channelheight");
    tarch::configuration::ParseConfiguration::readVector<3,double>(_cfg.wallVelocity,subtag,"wall-velocity");

    subtag = node->FirstChildElement("coupling");
    if (subtag == NULL){
      std::cout << "Could not read input file couette.xml: Missing subtag: coupling" << std::endl;
      exit(EXIT_FAILURE);
    }
    tarch::configuration::ParseConfiguration::readIntMandatory(_cfg.couplingCycles,subtag,"coupling-cycles");
    tarch::configuration::ParseConfiguration::readIntMandatory(_cfg.csvEveryTimestep,subtag,"write-csv-every-timestep");
    subtag = node->FirstChildElement("microscopic-solver");
    if (subtag == NULL){
      std::cout << "Could not read input file couette.xml: Missing subtag: microscopic-solver" << std::endl;
      exit(EXIT_FAILURE);
    }
    std::string type;
    tarch::configuration::ParseConfiguration::readStringMandatory(type,subtag,"type");
    if(type == "md"){
      tarch::configuration::ParseConfiguration::readDoubleMandatory(_cfg.temp,subtag,"temperature");
      tarch::configuration::ParseConfiguration::readIntMandatory(_cfg.equSteps,subtag,"equilibration-steps");
      tarch::configuration::ParseConfiguration::readIntMandatory(_cfg.totalNumberMDSimulations,subtag,"number-md-simulations");
      if(_cfg.totalNumberMDSimulations < 1){
        std::cout << "Could not read input file couette.xml: number-md-simulations < 1" << std::endl;
        exit(EXIT_FAILURE);
      }
    }
    else{
      std::cout << "Could not read input file couette.xml: Unknown microscopic solver type!" << std::endl;
      exit(EXIT_FAILURE);
    }
    tarch::configuration::ParseConfiguration::readDoubleMandatory(_cfg.density,subtag,"density");
    subtag = node->FirstChildElement("macroscopic-solver");
    if (subtag == NULL){
      std::cout << "Could not read input file couette.xml: Missing subtag: macroscopic-solver" << std::endl;
      exit(EXIT_FAILURE);
    }
    double vis;
    tarch::configuration::ParseConfiguration::readDoubleMandatory(vis, subtag, "viscosity");
    _cfg.kinVisc = vis / _cfg.density;
  }

  /** @brief initialises the macro and micro solver according to the setup from the xml file and pre-proccses them */
  void initSolvers(){
    // for timing measurements
    _tv.micro = 0;
    _tv.macro = 0;

    // allocate solvers
    _couetteSolver = NULL;
    _couetteSolver = getCouetteSolver( _mamicoConfig.getMacroscopicCellConfiguration().getMacroscopicCellSize()[0],
      _simpleMDConfig.getSimulationConfiguration().getDt()*_simpleMDConfig.getSimulationConfiguration().getNumberOfTimesteps()
      );
	 if (_couetteSolver != NULL) std::cout << "Couette solver not null on rank: " << _rank << std::endl; //TODO: remove debug

    _multiMDService = new tarch::utils::MultiMDService<3>(_simpleMDConfig.getMPIConfiguration().getNumberOfProcesses(),_cfg.totalNumberMDSimulations);
    _localMDInstances = _multiMDService->getLocalNumberOfMDSimulations();

    for (unsigned int i = 0; i < _localMDInstances; i++){
      _simpleMD.push_back(coupling::interface::SimulationAndInterfaceFactory::getInstance().getMDSimulation(
        _simpleMDConfig,_mamicoConfig
        #if (COUPLING_MD_PARALLEL==COUPLING_MD_YES)
        , _multiMDService->getLocalCommunicator()
        #endif
      ));
      if (_simpleMD[i]==NULL){
        std::cout << "ERROR CouetteTest: _simpleMD[" << i << "]==NULL!" << std::endl;
        exit(EXIT_FAILURE);
      }
    }
    _mdStepCounter = 0;
    if (_rank == 0){ gettimeofday(&_tv.start,NULL); }
    for (unsigned int i = 0; i < _localMDInstances; i++){
      _simpleMD[i]->init(*_multiMDService,_multiMDService->getGlobalNumberOfLocalMDSimulation(i));
    }
    // allocate coupling interfaces
    for (unsigned int i = 0; i < _localMDInstances; i++){
      _mdSolverInterface.push_back(coupling::interface::SimulationAndInterfaceFactory::getInstance().
        getMDSolverInterface(_simpleMDConfig, _mamicoConfig, _simpleMD[i]));
      if (_mdSolverInterface[i] == NULL){
        std::cout << "ERROR CouetteTest: mdSolverInterface[" << i << "] == NULL!" << std::endl;
        exit(EXIT_FAILURE);
      }
    }
    coupling::interface::MacroscopicSolverInterface<3>* couetteSolverInterface = getCouetteSolverInterface(
      _couetteSolver, _simpleMDConfig.getDomainConfiguration().getGlobalDomainOffset(),
      _mamicoConfig.getMacroscopicCellConfiguration().getMacroscopicCellSize(),
      getGlobalNumberMacroscopicCells(_simpleMDConfig,_mamicoConfig),_mamicoConfig.getMomentumInsertionConfiguration().getInnerOverlap()
    );
    _multiMDCellService = new coupling::services::MultiMDCellService<MY_LINKEDCELL,3>(
      _mdSolverInterface,couetteSolverInterface, _simpleMDConfig.getMPIConfiguration().getNumberOfProcesses(), (unsigned int) _rank, _cfg.totalNumberMDSimulations,
      _mamicoConfig.getParticleInsertionConfiguration(), _mamicoConfig.getMomentumInsertionConfiguration(), _mamicoConfig.getBoundaryForceConfiguration(),
      _mamicoConfig.getTransferStrategyConfiguration(), _mamicoConfig.getParallelTopologyConfiguration(),  _mamicoConfig.getThermostatConfiguration(), _simpleMDConfig.getSimulationConfiguration().getNumberOfTimesteps(),
      _mamicoConfig.getMacroscopicCellConfiguration(), "couette.xml", *_multiMDService
    );
    // set couette solver interface in MamicoInterfaceProvider
    coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL,3>::getInstance().setMacroscopicSolverInterface(couetteSolverInterface);
    for (unsigned int i = 0; i < _localMDInstances; i++){
      _simpleMD[i]->setMacroscopicCellService(&(_multiMDCellService->getMacroscopicCellService(i)));
      // compute and store temperature in macroscopic cells (temp=1.1 everywhere)
      _multiMDCellService->getMacroscopicCellService(i).computeAndStoreTemperature(_cfg.temp);
    }
    // allocate buffers for send/recv operations
    allocateSendBuffer(_multiMDCellService->getMacroscopicCellService(0).getIndexConversion(),*couetteSolverInterface);
    allocateRecvBuffer(_multiMDCellService->getMacroscopicCellService(0).getIndexConversion(),*couetteSolverInterface);
    // finish time measurement for initialisation
    if(_rank == 0){
      gettimeofday(&_tv.end,NULL);
      double runtime = (_tv.end.tv_sec - _tv.start.tv_sec)*1000000 + (_tv.end.tv_usec - _tv.start.tv_usec);
      std::cout << "Initialization: " << (int)(runtime/1000) << "ms" << std::endl;
    }
    if (_rank == 0){ gettimeofday(&_tv.start_total,NULL); }
    std::cout << "Finish CouetteTest::initSolvers() " << std::endl;
  }

  /** @brief advances the continuum solver and collects data to send to md (fillSendBuffer())
   *  @param cycle the number of the current coupling time step */
  void advanceMacro(int cycle){
    if(_couetteSolver != NULL){
      if (_rank==0){ gettimeofday(&_tv.start,NULL); }
      _couetteSolver->advance(_simpleMDConfig.getSimulationConfiguration().getDt()*_simpleMDConfig.getSimulationConfiguration().getNumberOfTimesteps());
      if (_rank==0){
        gettimeofday(&_tv.end,NULL);
        _tv.macro += (_tv.end.tv_sec - _tv.start.tv_sec)*1000000 + (_tv.end.tv_usec - _tv.start.tv_usec);
        //std::cout << "Finish _couetteSolver->advance " << std::endl;
      }
    }
  }

  /** @brief advances the md solver for one coupling time step and collect the data for the coupling
   *  @param cycle the number of the current coupling time step  */
  void advanceMicro(int cycle){
    if (_rank==0){ gettimeofday(&_tv.start,NULL); }
    // run MD instances
    for (unsigned int i = 0; i < _localMDInstances; i++){
      // set macroscopic cell service and interfaces in MamicoInterfaceProvider
      coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL,3>::getInstance().setMacroscopicCellService(&(_multiMDCellService->getMacroscopicCellService(i)));
      coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL,3>::getInstance().setMDSolverInterface(_mdSolverInterface[i]);
      _simpleMD[i]->simulateTimesteps(_simpleMDConfig.getSimulationConfiguration().getNumberOfTimesteps(),_mdStepCounter);
      // plot macroscopic time step info in multi md service
      _multiMDCellService->getMacroscopicCellService(i).plotEveryMacroscopicTimestep(cycle);
    }
    _mdStepCounter += _simpleMDConfig.getSimulationConfiguration().getNumberOfTimesteps();

    if (_rank==0){
      gettimeofday(&_tv.end,NULL);
      _tv.micro += (_tv.end.tv_sec - _tv.start.tv_sec)*1000000 + (_tv.end.tv_usec - _tv.start.tv_usec);
    }
    // send back data from MD instances and merge it
    _multiMDCellService->sendFromMD2Macro(_buf.recvBuffer,_buf.globalCellIndices4RecvBuffer);
    write2CSV(_buf.recvBuffer,_buf.globalCellIndices4RecvBuffer,_multiMDCellService->getMacroscopicCellService(0).getIndexConversion(),cycle+1);
    std::cout << "Finished " << cycle << std::endl;
  }

  /** @brief finalize the time measurement, and cleans up at the end of the simulation */
  void shutdown(){
    // finish time measurement for coupled simulation
    if (_rank==0){
      gettimeofday(&_tv.end,NULL);
      double time_total = (_tv.end.tv_sec - _tv.start_total.tv_sec)*1000000 + (_tv.end.tv_usec - _tv.start_total.tv_usec);
      std::cout << "Finished all coupling cycles after " << time_total/1000000 << " s" << std::endl;
      std::cout << "Time percentages Micro, Macro, Other: " << std::endl;
      std::cout << _tv.micro/time_total*100 << ", " << _tv.macro/time_total*100 << ",  "
                << (1-(_tv.micro+_tv.macro)/time_total)*100 << std::endl;
    }

    // shutdown MD simulation
    for (unsigned int i = 0; i < _localMDInstances; i++) {
      // the shutdown operation may also delete the md solver interface; therefore, we update the MD solver interface in the vector _mdSolverInteface after the shutdown is completed
      coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL,3>::getInstance().setMDSolverInterface(_mdSolverInterface[i]);
      _simpleMD[i]->shutdown(); delete _simpleMD[i]; _simpleMD[i] = NULL;
      _mdSolverInterface[i] = coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL,3>::getInstance().getMDSolverInterface();
    }
    _simpleMD.clear();
    for (unsigned int i = 0; i < _localMDInstances; i++) {
      if (_mdSolverInterface[i] != NULL){
        delete _mdSolverInterface[i];
        _mdSolverInterface[i] = NULL;
      }
    }
    _mdSolverInterface.clear();

    coupling::interface::MacroscopicSolverInterface<3>* couetteSolverInterface =
      coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL,3>::getInstance().getMacroscopicSolverInterface();

    if( _multiMDService != NULL){delete _multiMDService; _multiMDService = NULL;}
    if (couetteSolverInterface != NULL){delete couetteSolverInterface; couetteSolverInterface = NULL;}
    if (_couetteSolver != NULL){delete _couetteSolver; _couetteSolver=NULL;}
    if(_multiMDCellService != NULL){delete _multiMDCellService; _multiMDCellService=NULL;}

    std::cout << "Finish CouetteTest::shutdown() " << std::endl;
  }

  /** computes global number of macroscopic cells from configs. Required by couette solver interface before MacroscopicCellService is initialised! */
  tarch::la::Vector<3,unsigned int> getGlobalNumberMacroscopicCells(
  const simplemd::configurations::MolecularDynamicsConfiguration &simpleMDConfig, const coupling::configurations::MaMiCoConfiguration<3> &mamicoConfig) const {
    tarch::la::Vector<3,double> domainSize(simpleMDConfig.getDomainConfiguration().getGlobalDomainSize());
    tarch::la::Vector<3,double> dx(mamicoConfig.getMacroscopicCellConfiguration().getMacroscopicCellSize());
    tarch::la::Vector<3,unsigned int> globalNumberMacroscopicCells(0);
    for (unsigned int d = 0;  d < 3; d++){
     int buf = floor(domainSize[d]/dx[d]+0.5);
     globalNumberMacroscopicCells[d] = (unsigned int) buf;
    }
    return globalNumberMacroscopicCells;
  }

  /** This is only done on rank 0.
   *  @brief allocates the send buffer (with values for all macroscopic cells).
   *  @param indexConversion instance of the indexConversion
   *  @param couetteSolverInterface interface for the continuum solver */
  void allocateSendBuffer(const coupling::IndexConversion<3>& indexConversion,
  coupling::interface::MacroscopicSolverInterface<3> &couetteSolverInterface) {
    // determine global number of cells
    const tarch::la::Vector<3,unsigned int> cells(indexConversion.getGlobalNumberMacroscopicCells()+tarch::la::Vector<3,unsigned int>(2));
    const unsigned int num = cells[0]*cells[1]*cells[2];
    // delete all potential entries of sendBuffer
    deleteBuffer(_buf.sendBuffer);
    // count number of cells to be sent from this process; therefore, loop over all global macroscopic cells...
    unsigned int numCellsSent=0;
    for (unsigned int i =0; i < num; i++){
       // ... and find out, if the current cell should be send to MD from this couette solver process
       if(couetteSolverInterface.sendMacroscopicQuantityToMDSolver(indexConversion.getGlobalVectorCellIndex(i))){
         std::vector<unsigned int> ranks = couetteSolverInterface.getSourceRanks(indexConversion.getGlobalVectorCellIndex(i));
         bool containsThisRank=false;
         for (unsigned int k = 0; k < ranks.size(); k++){
           containsThisRank = containsThisRank || (ranks[k]==(unsigned int)_rank);
         }
         if (containsThisRank){ numCellsSent++; }
      }
    }
    // allocate array for cell indices
    unsigned int* indices = new unsigned int [numCellsSent];
    if (indices==NULL){std::cout << "ERROR CouetteTest::allocateSendBuffer(): indices==NULL!" << std::endl; exit(EXIT_FAILURE); }
    // allocate sendBuffer and initialise all entries, incl. indices
    for (unsigned int i = 0; i < num; i++){
      if (couetteSolverInterface.sendMacroscopicQuantityToMDSolver(indexConversion.getGlobalVectorCellIndex(i))){
         std::vector<unsigned int> ranks = couetteSolverInterface.getSourceRanks(indexConversion.getGlobalVectorCellIndex(i));
         bool containsThisRank=false;
         for (unsigned int k = 0; k < ranks.size(); k++){
           containsThisRank = containsThisRank || (ranks[k]==(unsigned int)_rank);
         }
         if (containsThisRank){
           _buf.sendBuffer.push_back(new coupling::datastructures::MacroscopicCell<3>());
           if (_buf.sendBuffer[_buf.sendBuffer.size()-1]==NULL){std::cout << "ERROR CouetteTest::allocateSendBuffer: sendBuffer[" << _buf.sendBuffer.size()-1 << "]==NULL!" << std::endl; exit(EXIT_FAILURE);}
           indices[_buf.sendBuffer.size()-1] = i;
         }
      }
    }
    #if (COUPLING_MD_DEBUG==COUPLING_MD_YES)
    for (unsigned int i = 0; i < numCellsSent; i++){
      std::vector<unsigned int> ranks = couetteSolverInterface.getSourceRanks(indexConversion.getGlobalVectorCellIndex(indices[i]));
      std::cout << "Current rank= " << _rank << ", Send cell " << indexConversion.getGlobalVectorCellIndex(indices[i]) << ", ranks=";
      for (unsigned int j = 0; j < ranks.size(); j++){ std::cout << " " << ranks[j];}
      std::cout << std::endl;
    }
    #endif
    _buf.globalCellIndices4SendBuffer = indices;
  }

  /** allocates the recv-buffer. This buffer contains all global inner macroscopic cells, but only on rank 0. On all other ranks, no cells are stored and a NULL ptr is returned */
  void allocateRecvBuffer(const coupling::IndexConversion<3>& indexConversion,
  coupling::interface::MacroscopicSolverInterface<3> &couetteSolverInterface) {
    // determine global number of cells
    const tarch::la::Vector<3,unsigned int> cells(indexConversion.getGlobalNumberMacroscopicCells()+tarch::la::Vector<3,unsigned int>(2));
    const unsigned int num = cells[0]*cells[1]*cells[2];
    // delete all potential entries of sendBuffer
    deleteBuffer(_buf.recvBuffer);
    // determine number of cells that should be received
    unsigned int numCellsRecv = 0;
    for (unsigned int i = 0; i < num; i++){
      if(couetteSolverInterface.receiveMacroscopicQuantityFromMDSolver(indexConversion.getGlobalVectorCellIndex(i))){
        std::vector<unsigned int> ranks = couetteSolverInterface.getTargetRanks(indexConversion.getGlobalVectorCellIndex(i));
        bool containsThisRank=false;
        for (unsigned int k = 0; k < ranks.size(); k++){
          containsThisRank = containsThisRank || (ranks[k]==(unsigned int)_rank);
        }
        if (containsThisRank){ numCellsRecv++; }
      }
    }
    // allocate array for cell indices
    unsigned int* indices = new unsigned int [numCellsRecv];
    if (indices==NULL){std::cout << "ERROR CouetteTest::allocateRecvBuffer(): indices==NULL!" << std::endl; exit(EXIT_FAILURE); }
    // allocate recvBuffer and initialise all entries, incl. indices
    for (unsigned int i = 0; i < num; i++){
      if (couetteSolverInterface.receiveMacroscopicQuantityFromMDSolver(indexConversion.getGlobalVectorCellIndex(i))){
        std::vector<unsigned int> ranks = couetteSolverInterface.getTargetRanks(indexConversion.getGlobalVectorCellIndex(i));
        bool containsThisRank=false;
        for (unsigned int k = 0; k < ranks.size(); k++){
          containsThisRank = containsThisRank || (ranks[k]==(unsigned int)_rank);
        }
        if (containsThisRank){
          _buf.recvBuffer.push_back(new coupling::datastructures::MacroscopicCell<3>());
          if (_buf.recvBuffer[_buf.recvBuffer.size()-1]==NULL){std::cout << "ERROR CouetteTest::allocateRecvBuffer: recvBuffer[" << _buf.recvBuffer.size()-1 << "]==NULL!" << std::endl; exit(EXIT_FAILURE);}
          // set linearized index
          indices[_buf.recvBuffer.size()-1] = i;
        }
      }
    }
    #if (COUPLING_MD_DEBUG==COUPLING_MD_YES)
    for (unsigned int i = 0; i < numCellsRecv; i++){
      std::vector<unsigned int> ranks = couetteSolverInterface.getTargetRanks(indexConversion.getGlobalVectorCellIndex(indices[i]));
      std::cout << "Current rank= " << _rank << ", Recv cell " << indexConversion.getGlobalVectorCellIndex(indices[i]) << ", ranks=";
      for (unsigned int j = 0; j < ranks.size(); j++){ std::cout << " " << ranks[j];}
      std::cout << std::endl;
    }
    #endif
    _buf.globalCellIndices4RecvBuffer = indices;
  }

  /** @brief write macroscopic cells that have been received from MD to csv file
   *  @param recvBuffer the buffer for the data, which comes from md
   *  @param recvIndices the indices for the macr cells in the buffer
   *  @param indexConversion an instance of the indexConversion
   *  @param couplingCycle the current number of coupling cycle */
  void write2CSV(
    std::vector<coupling::datastructures::MacroscopicCell<3>* >& recvBuffer,const unsigned int * const recvIndices,
    const coupling::IndexConversion<3>& indexConversion, int couplingCycle
  ) const {
    if(recvBuffer.size() == 0) return;
    if(_cfg.csvEveryTimestep < 1 || couplingCycle % _cfg.csvEveryTimestep > 0) return;
    // form file name and open file
    std::stringstream ss;
    ss << "CouetteAvgMultiMDCells_" << _rank << "_" << couplingCycle << ".csv";
    std::ofstream file(ss.str().c_str());
    if (!file.is_open()){std::cout << "ERROR CouetteTest::write2CSV(): Could not open file " << ss.str() << "!" << std::endl; exit(EXIT_FAILURE);}

    // loop over received cells; read macroscopic mass+momentum buffers and write cell index, mass and velocity to one line in the csv-file
    const unsigned int numCellsRecv = recvBuffer.size();
    for (unsigned int i = 0; i < numCellsRecv; i++){
      tarch::la::Vector<3,double> vel(recvBuffer[i]->getMacroscopicMomentum());
      if (recvBuffer[i]->getMacroscopicMass()!=0.0){ vel = (1.0/recvBuffer[i]->getMacroscopicMass())*vel; }
      const tarch::la::Vector<3,unsigned int> counter(indexConversion.getGlobalVectorCellIndex(recvIndices[i]));
      file << counter[0] << " ; " << counter[1] << " ; " << counter[2] << " ; " << vel[0] << " ; " << vel[1] << " ; " << vel[2] << " ; " << recvBuffer[i]->getMacroscopicMass();
      file << std::endl;
    }

    // close file
    file.close();
  }


  /** @brief deletes the data in the buffer for the macro to md transfer
   *  @param sendBuffer the buffer to be cleaned */
  void deleteBuffer(std::vector<coupling::datastructures::MacroscopicCell<3>* >& sendBuffer) const {
    // delete all potential entries of sendBuffer
    for (unsigned int i = 0; i < sendBuffer.size(); i++){if (sendBuffer[i]!=NULL){ delete sendBuffer[i]; sendBuffer[i]=NULL;}}
    sendBuffer.clear();
  }

  /** @brief fills send buffer with data from macro/continuum solver
   *  @param density the general density of the fluid
   *  @param couetteSolver the continuum solver
   *  @param indexConversion an instance of the indexConversion
   *  @param sendBuffer the bufffer to send data from macro to micro
   *  @param globalCellIndices4SendBuffer the global linearized indices of the macroscopic cells in the buffer  */
  void fillSendBuffer(
    const double density, const coupling::solvers::AbstractCouetteSolver<3>& couetteSolver, const coupling::IndexConversion<3>& indexConversion,
    std::vector<coupling::datastructures::MacroscopicCell<3>* >& sendBuffer, const unsigned int * const globalCellIndices4SendBuffer
  ) const {
    const unsigned int size = sendBuffer.size();
    const tarch::la::Vector<3,double> domainOffset(indexConversion.getGlobalMDDomainOffset());
    const tarch::la::Vector<3,double> macroscopicCellSize(indexConversion.getMacroscopicCellSize());

    for (unsigned int i = 0; i < size; i++){
      // get global cell index vector
      const tarch::la::Vector<3,unsigned int> globalIndex(indexConversion.getGlobalVectorCellIndex(globalCellIndices4SendBuffer[i]));
      // determine cell midpoint
      tarch::la::Vector<3,double> cellMidPoint(domainOffset-0.5*macroscopicCellSize);
      for (unsigned int d = 0; d < 3; d++){ cellMidPoint[d] = cellMidPoint[d] + ((double)globalIndex[d])*macroscopicCellSize[d]; }

      double mass = density * macroscopicCellSize[0]*macroscopicCellSize[1]*macroscopicCellSize[2];
      mass *= static_cast<const coupling::solvers::LBCouetteSolver*>(&couetteSolver)->getDensity(cellMidPoint);

      // compute momentum
      tarch::la::Vector<3,double> momentum(mass*couetteSolver.getVelocity(cellMidPoint));
      sendBuffer[i]->setMicroscopicMass(mass);
      sendBuffer[i]->setMicroscopicMomentum(momentum);
    }
  }

  /** @returns the correct marco/continuum solver for the given setup
   *  @param dx the grid size (equidistant mesh)
   *  @param dt the time step  */
  coupling::solvers::AbstractCouetteSolver<3>* getCouetteSolver(const double dx, const double dt){
    coupling::solvers::AbstractCouetteSolver<3>* solver = NULL;
    // LB solver: active on lbNumberProcesses
    solver = new coupling::solvers::LBCouetteSolver(_cfg.channelheight,_cfg.wallVelocity,_cfg.kinVisc,dx,dt,-1,"LBCouette",tarch::la::Vector<3,unsigned int>(1));
    if (solver==NULL){
      std::cout << "ERROR CouetteTest::getCouetteSolver(): LB solver==NULL!" << std::endl;
      exit(EXIT_FAILURE);
    }
    return solver;
  }

  /** @returns the interface for the macro/continuum solver
   *  @param couetteSolver the macro/continuum solver
   *  @param mdOffset the offset of the md domain from (0.0.0)
   *  @param mamicoMeshsize
   *  @param globalNumberMacroscopicCells the total number macroscopic cells for the whole domain
   *  @param outerRegion
   *  @todo piet, what is the mamicoMeshsize & the outer layer  */
  coupling::interface::MacroscopicSolverInterface<3>* getCouetteSolverInterface(
    coupling::solvers::AbstractCouetteSolver<3>* couetteSolver,
    tarch::la::Vector<3,double> mdOffset,
    tarch::la::Vector<3,double> mamicoMeshsize,
    tarch::la::Vector<3,unsigned int> globalNumberMacroscopicCells,unsigned int outerRegion
  ){
    coupling::interface::MacroscopicSolverInterface<3>* interface = NULL;
    coupling::solvers::LBCouetteSolver *lbSolver = static_cast<coupling::solvers::LBCouetteSolver*>(couetteSolver);
    if (lbSolver==NULL){std::cout << "ERROR CouetteTest::getCouetteSolverInterface(...), rank=" << _rank << ": Could not convert abstract to LB solver!" << std::endl; exit(EXIT_FAILURE);}
    // compute number of cells of MD offset; detect any mismatches!
    tarch::la::Vector<3,unsigned int> offsetMDDomain(0);
    for (unsigned int d = 0; d < 3; d++){
      if (mdOffset[d] < 0.0){std::cout << "ERROR CouetteTest::getCouetteSolverInterface(...): mdOffset[" << d << "]<0.0!" << std::endl; exit(EXIT_FAILURE);}
      offsetMDDomain[d] = floor(mdOffset[d]/mamicoMeshsize[d] + 0.5);
      if (fabs( (offsetMDDomain[d]*mamicoMeshsize[d] - mdOffset[d])/mamicoMeshsize[d])>1.0e-8){
        std::cout << "ERROR CouetteTest::getCouetteSolverInterface: MD offset and mesh size mismatch!" << std::endl;
        exit(EXIT_FAILURE);
      }
    }
    interface = new coupling::solvers::CouetteSolverInterface<3>(globalNumberMacroscopicCells,outerRegion);
    if (interface==NULL){std::cout << "ERROR CouetteTest::getCouetteSolverInterface(...), rank=" << _rank << ": interface==NULL!" << std::endl; exit(EXIT_FAILURE);}
    return interface;
  }

  /** @brief all the variables necessary to define the couette scenario are stored in here */
  struct EvaporationConfig{
    /** @brief channel is always expected to have origin at (0.0,0.0,0.0) and to have the same width as height */
    double channelheight;
    /** @brief the total length of the channel */
    double channellength;
    /** @brief velocity of moving wall (lower boundary moves) */
    tarch::la::Vector<3,double> wallVelocity;
    /** @brief number of coupling cycles, that is continuum time steps; MD/DPD: 1000 */
    int couplingCycles;
    /** @brief the time step interval for writing the md data to csv files  */
    int csvEveryTimestep;
    /** @brief the general density of the fluid under consideration */
    double density;
    /** @brief the kinematic viscosity of the fluid under consideration */
    double kinVisc;
    /** @brief the start temperature for the fluid under consideration */
    double temp;
    /** @brief number of equilibartion time steps = number of time steps that the md will run before the coupling is enabled */
    int equSteps;
    /** @brief the number of md simulation instances in a multi-instance coupling  */
    int totalNumberMDSimulations;
  };

  /** the buffers store macroscopic cells, so momentum and density will be transferred
   *  @brief holds the buffers for the data transfer */
  struct CouplingBuffer{
    /** @brief the buffer for data transfer from macro to md */
    std::vector<coupling::datastructures::MacroscopicCell<3>* > sendBuffer;
    /** @brief the global indices of the macroscopic cells in the sendBuffer */
    unsigned int *globalCellIndices4SendBuffer;
    /** @brief the buffer for data transfer from md to macro */
    std::vector<coupling::datastructures::MacroscopicCell<3>* > recvBuffer;
    /** @brief the global indices of the macroscopic cells in the recvBuffer*/
    unsigned int *globalCellIndices4RecvBuffer;
  };

  /** @brief holds all the variables for the time measurement of a simulation
   *  @todo Piet*/
  struct TimingValues{
    /** @brief */
    timeval start_total;
    /** @brief */
    timeval start;
    /** @brief */
    timeval end;
    /** @brief */
    double micro;
    /** @brief */
    double macro;
  };

  /** @brief the rank of the current MPI process */
  int _rank;
  /** @brief the config data and information for SimpleMD */
  simplemd::configurations::MolecularDynamicsConfiguration _simpleMDConfig;
  /** @brief the config data and information for MaMiCo*/
  coupling::configurations::MaMiCoConfiguration<3> _mamicoConfig;
  /** @brief the CouetteConfig for the current setup */
  EvaporationConfig _cfg;
  /** @brief the counter for the time steps, which are done within the md */
  unsigned int _mdStepCounter;
  /** @brief the current macro/continuum solver */
  coupling::solvers::AbstractCouetteSolver<3> *_couetteSolver;
  /** @todo piet */
  tarch::utils::MultiMDService<3>* _multiMDService;
  /** @todo piet */
  coupling::services::MultiMDCellService<MY_LINKEDCELL,3> *_multiMDCellService;
  /** @brief the current buffer for the coupling */
  CouplingBuffer _buf;
  /** @todo piet */
  unsigned int _localMDInstances;
  /** @brief the interface to the md solver */
  std::vector<coupling::interface::MDSolverInterface<MY_LINKEDCELL,3>* > _mdSolverInterface;
  /** @brief the md solver */
  std::vector<coupling::interface::MDSimulation*> _simpleMD;
  /** @todo piet */
  std::default_random_engine _generator;
  TimingValues _tv;
};

#endif // _COUPLING_TESTS_COUETTETEST_H_
