// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _COUPLING_TESTS_COUETTETEST_H_
#define _COUPLING_TESTS_COUETTETEST_H_

#include "tarch/utils/MultiMDService.h"
#include "simplemd/configurations/MolecularDynamicsConfiguration.h"
#include "coupling/tests/Test.h"
#include "coupling/CouplingMDDefinitions.h"
#include "tarch/configuration/ParseConfiguration.h"
#include "coupling/solvers/CouetteSolver.h"
#include "coupling/solvers/LBCouetteSolver.h"
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

/** 
 * Versatile configurable Couette flow test for noise-filtered multi-instance Nie coupling. 
 *
 * Features, depending on couette.xml config file:
 * -> one-way or two-way coupling
 * -> statup or oscillating flow scenario
 * -> using analytical or Lattice Bolzmann Couette flow solver
 * -> using (multi-instance) SimpleMD or synthetic MD data (with Gaussian noise)
 * -> CSV output option for filtered data
 * -> ability to delay scenario startup or two-way coupling for filter initialisation time
 * -> signal-to-noise ratio (SNR) computation (between filter output and macroscopic solver)
 * -> runtime measurements and logging separately for coupled simulation components
 * -> 'tws-loop' feature for testing many POD time-window-size parameters in a single run
 *
 *  @author Piet Jarmatz
 */
class CouetteTest: public Test {
public:
  CouetteTest(): Test("CouetteTest"){}
  virtual ~CouetteTest(){}

  virtual void run(){
    init();
    for (int cycle = 0; cycle < _cfg.couplingCycles; cycle++)
      runOneCouplingCycle(cycle);
    shutdown();
  }

private:
  enum MacroSolverType{COUETTE_ANALYTICAL=0,COUETTE_LB=1};
  enum MicroSolverType{SIMPLEMD=0,SYNTHETIC=1};
  
  void init(){
    initMPI();
    parseConfigurations();

    // for time measurements
    //timeval start;
    //timeval end;

    // allocate solvers
    coupling::solvers::AbstractCouetteSolver<3> *couetteSolver = NULL;
    couetteSolver = getCouetteSolver( _mamicoConfig.getMacroscopicCellConfiguration().getMacroscopicCellSize()[0],
      _simpleMDConfig.getSimulationConfiguration().getDt()*_simpleMDConfig.getSimulationConfiguration().getNumberOfTimesteps()
      );
    tarch::utils::MultiMDService<3> multiMDService(_simpleMDConfig.getMPIConfiguration().getNumberOfProcesses(),_cfg.totalNumberMDSimulations);
    const unsigned int mdInstances = multiMDService.getLocalNumberOfMDSimulations(); // local number of MD simulations

    std::vector<coupling::interface::MDSimulation*> simpleMD;
    for (unsigned int i = 0; i < mdInstances; i++){
      simpleMD.push_back(coupling::interface::SimulationAndInterfaceFactory::getInstance().getMDSimulation(
        _simpleMDConfig,_mamicoConfig
        #if (COUPLING_MD_PARALLEL==COUPLING_MD_YES)
        , multiMDService.getLocalCommunicator()
        #endif
      ));
      if (simpleMD[i]==NULL){
        std::cout << "ERROR CouetteTest: simpleMD[" << i << "]==NULL!" << std::endl; 
        exit(EXIT_FAILURE);
      }
    }

    // equilibrate MD
    _mdStepCounter = 0;
    //if (_rank == 0){ gettimeofday(&start,NULL); }
    for (unsigned int i = 0; i < mdInstances; i++){
      simpleMD[i]->init(multiMDService,multiMDService.getGlobalNumberOfLocalMDSimulation(i));
      simpleMD[i]->switchOffCoupling();
      simpleMD[i]->simulateTimesteps(_cfg.equSteps,_mdStepCounter);
    }
    _mdStepCounter += _cfg.equSteps;
    // finish time measurement for initialisation(incl equilibration) and start time measurement for coupling
    //if(rank==0){
    //  gettimeofday(&end,NULL);
    //  double runtime = (end.tv_sec-start.tv_sec)*1000000 + (end.tv_usec-start.tv_usec);
    //  std::cout << "NieTest-Equilibration: " << (int)(runtime/1000) << "ms" << std::endl;
    //}

    // allocate coupling interfaces
    std::vector<coupling::interface::MDSolverInterface<MY_LINKEDCELL,3>* > mdSolverInterface;
    for (unsigned int i = 0; i < mdInstances; i++){
      simpleMD[i]->switchOnCoupling();
      mdSolverInterface.push_back(coupling::interface::SimulationAndInterfaceFactory::getInstance().
        getMDSolverInterface(_simpleMDConfig, _mamicoConfig, simpleMD[i]));
      if (mdSolverInterface[i] == NULL){
        std::cout << "ERROR CouetteTest: mdSolverInterface[" << i << "] == NULL!" << std::endl; 
        exit(EXIT_FAILURE);
      }
    }
    coupling::interface::MacroscopicSolverInterface<3>* couetteSolverInterface = getCouetteSolverInterface(
      couetteSolver, _simpleMDConfig.getDomainConfiguration().getGlobalDomainOffset(),
      _mamicoConfig.getMacroscopicCellConfiguration().getMacroscopicCellSize(),
      getGlobalNumberMacroscopicCells(_simpleMDConfig,_mamicoConfig),_mamicoConfig.getMomentumInsertionConfiguration().getInnerOverlap()
    );

    // initialise macroscopic cell service for multi-MD case and set single cell services in each MD simulation
    coupling::services::MultiMDCellService<MY_LINKEDCELL,3> multiMDCellService(
      mdSolverInterface,couetteSolverInterface, _simpleMDConfig.getMPIConfiguration().getNumberOfProcesses(), (unsigned int) _rank, _cfg.totalNumberMDSimulations,
      _mamicoConfig.getParticleInsertionConfiguration(), _mamicoConfig.getMomentumInsertionConfiguration(), _mamicoConfig.getBoundaryForceConfiguration(),
      _mamicoConfig.getTransferStrategyConfiguration(), _mamicoConfig.getNoiseReductionConfiguration(), _mamicoConfig.getParallelTopologyConfiguration(), _simpleMDConfig.getSimulationConfiguration().getNumberOfTimesteps(),
      _mamicoConfig.getMacroscopicCellConfiguration(), multiMDService
    );
    // set couette solver interface in MamicoInterfaceProvider
    coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL,3>::getInstance().setMacroscopicSolverInterface(couetteSolverInterface);

    for (unsigned int i = 0; i < mdInstances; i++){
      simpleMD[i]->setMacroscopicCellService(&(multiMDCellService.getMacroscopicCellService(i)));
      // compute and store temperature in macroscopic cells (temp=1.1 everywhere)
      multiMDCellService.getMacroscopicCellService(i).computeAndStoreTemperature(_cfg.temp);
    }

    // allocate buffers for send/recv operations
    std::vector<coupling::datastructures::MacroscopicCell<3>* > sendBuffer;
    unsigned int *globalCellIndices4SendBuffer = allocateSendBuffer(sendBuffer,multiMDCellService.getMacroscopicCellService(0).getIndexConversion(),_rank,*couetteSolverInterface);
    std::vector<coupling::datastructures::MacroscopicCell<3>* > recvBuffer;
    unsigned int *globalCellIndices4RecvBuffer = allocateRecvBuffer(recvBuffer,multiMDCellService.getMacroscopicCellService(0).getIndexConversion(),_rank,*couetteSolverInterface);

    //std::default_random_engine generator(0);

    //if (rank==0){ gettimeofday(&start,NULL); }
  }

  void initMPI(){
    _rank = 0;
    #if (COUPLING_MD_PARALLEL==COUPLING_MD_YES)
    MPI_Comm_rank(MPI_COMM_WORLD,&_rank);
    #endif
  }

  void parseConfigurations(){
    tarch::configuration::ParseConfiguration::parseConfiguration<simplemd::configurations::MolecularDynamicsConfiguration>("couette_simplemd.xml","molecular-dynamics",_simpleMDConfig);
    if (!_simpleMDConfig.isValid()){std::cout << "ERROR NieTest: Invalid SimpleMD config!" << std::endl; exit(EXIT_FAILURE);}
    
    tarch::configuration::ParseConfiguration::parseConfiguration<coupling::configurations::MaMiCoConfiguration<3> >("couette_mamico.xml","mamico",_mamicoConfig);
    if (!_mamicoConfig.isValid()){ std::cout << "ERROR NieTest: Invalid MaMiCo config!" << std::endl; exit(EXIT_FAILURE); }

    parseCouetteTestConfiguration();
  }

  void parseCouetteTestConfiguration(){
    tinyxml2::XMLDocument conffile;
    tinyxml2::XMLElement *node = NULL;
    conffile.LoadFile("couette.xml");
    node = conffile.FirstChildElement("couette-test");
    if (node == NULL){
      std::cout << "Could not read input file couette.xml" << std::endl;
      exit(EXIT_FAILURE);
    }

    tinyxml2::XMLElement* subtag = node->FirstChildElement("domain");
    if (subtag == NULL){
      std::cout << "Could not read input file couette.xml: Missing subtag: domain" << std::endl;
      exit(EXIT_FAILURE);
    }
    tarch::configuration::ParseConfiguration::readDoubleMandatory(_cfg.channelheight,subtag,"channelheight");
    tarch::configuration::ParseConfiguration::readVector<3,double>(_cfg.wallVelocity,subtag,"wall-velocity");
    _cfg.wallInitCycles = 0;
    tarch::configuration::ParseConfiguration::readIntOptional(_cfg.wallInitCycles,subtag,"wall-init-cycles");
    if(_cfg.wallInitCycles > 0)
      tarch::configuration::ParseConfiguration::readVector<3,double>(_cfg.wallInitVelocity,subtag,"wall-init-velocity");
    _cfg.wallOscillations = 0;
    tarch::configuration::ParseConfiguration::readDoubleOptional(_cfg.wallOscillations,subtag,"wall-oscillations");

    subtag = node->FirstChildElement("coupling");
    if (subtag == NULL){
      std::cout << "Could not read input file couette.xml: Missing subtag: coupling" << std::endl;
      exit(EXIT_FAILURE);
    }
    tarch::configuration::ParseConfiguration::readIntMandatory(_cfg.couplingCycles,subtag,"coupling-cycles");
    tarch::configuration::ParseConfiguration::readBoolMandatory(_cfg.md2Macro,subtag,"send-from-md-to-macro");
    tarch::configuration::ParseConfiguration::readBoolMandatory(_cfg.macro2Md,subtag,"send-from-macro-to-md");
    _cfg.filterInitCycles = 0;
    tarch::configuration::ParseConfiguration::readIntOptional(_cfg.filterInitCycles,subtag,"filter-init-cycles");
    _cfg.csvEveryTimestep = 0;
    tarch::configuration::ParseConfiguration::readIntOptional(_cfg.csvEveryTimestep,subtag,"write-csv-every-timestep");
    _cfg.computeSNR = false;
    tarch::configuration::ParseConfiguration::readBoolOptional(_cfg.computeSNR,subtag,"compute-snr");

    subtag = node->FirstChildElement("microscopic-solver");
    if (subtag == NULL){
      std::cout << "Could not read input file couette.xml: Missing subtag: microscopic-solver" << std::endl;
      exit(EXIT_FAILURE);
    }
    std::string type;
    tarch::configuration::ParseConfiguration::readStringMandatory(type,subtag,"type");
    if(type == "md"){ 
      _cfg.miSolverType = SIMPLEMD;
      tarch::configuration::ParseConfiguration::readDoubleMandatory(_cfg.temp,subtag,"temperature");
      _cfg.equSteps = 0;
      tarch::configuration::ParseConfiguration::readIntOptional(_cfg.equSteps,subtag,"equilibration-steps");
      _cfg.totalNumberMDSimulations = 1;
      tarch::configuration::ParseConfiguration::readIntOptional(_cfg.totalNumberMDSimulations,subtag,"number-md-simulations");
    }
    else if(type == "synthetic"){
      _cfg.miSolverType = SYNTHETIC;
      tarch::configuration::ParseConfiguration::readDoubleMandatory(_cfg.noiseSigma,subtag,"noise-sigma");
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
    tarch::configuration::ParseConfiguration::readStringMandatory(type,subtag,"type");
    if(type == "lb"){ 
      _cfg.maSolverType = COUETTE_LB;
      tarch::configuration::ParseConfiguration::readVector<3,unsigned int>(_cfg.lbNumberProcesses,subtag,"number-of-processes");
      _cfg.plotEveryTimestep = 0;
      tarch::configuration::ParseConfiguration::readIntOptional(_cfg.plotEveryTimestep,subtag,"plot-every-timestep");
    }
    else if(type == "analytical"){
      _cfg.maSolverType = COUETTE_ANALYTICAL;
      if(!(_cfg.wallVelocity[1] == 0.0 && _cfg.wallVelocity[2] == 0.0)){
        std::cout << "analytic solver only supports flow in x-direction" << std::endl;
        exit(EXIT_FAILURE);
      }
    }
    else{
      std::cout << "Could not read input file couette.xml: Unknown macroscopic solver type!" << std::endl;
      exit(EXIT_FAILURE);
    }
    double vis;
    tarch::configuration::ParseConfiguration::readDoubleMandatory(vis, subtag, "viscosity");
    _cfg.kinVisc = vis / _cfg.density;
    _cfg.initAdvanceCycles = 0;
    tarch::configuration::ParseConfiguration::readIntOptional(_cfg.initAdvanceCycles,subtag,"init-advance-cycles"); 

    subtag = node->FirstChildElement("tws-loop");
    if (subtag == NULL) _cfg.twsLoop = false;
    else{
      _cfg.twsLoop = true;
      tarch::configuration::ParseConfiguration::readIntMandatory(_cfg.twsLoopMin,subtag,"min");
      tarch::configuration::ParseConfiguration::readIntMandatory(_cfg.twsLoopMax,subtag,"max");
      _cfg.twsLoopStep = 1;
      tarch::configuration::ParseConfiguration::readIntOptional(_cfg.twsLoopStep,subtag,"step");
    }
  }

  void runOneCouplingCycle(unsigned int cycle){/*
    // run one time step for couette solver
      if (couetteSolver!=NULL){
        if ( (solverType==COUETTE_LB) && cycles == wallInitCycles){
          static_cast<coupling::solvers::LBCouetteSolver*>(couetteSolver)->setWallVelocity(wallVelocity2);
        }
        couetteSolver->advance(simpleMDConfig.getSimulationConfiguration().getDt()*simpleMDConfig.getSimulationConfiguration().getNumberOfTimesteps());
        std::cout << "Finish couetteSolver->advance " << std::endl;
        //couetteSolver2->advance(simpleMDConfig.getSimulationConfiguration().getDt()*simpleMDConfig.getSimulationConfiguration().getNumberOfTimesteps());
      }
      // extract data from couette solver and send them to MD (can take any index-conversion object)
      fillSendBuffer(density,*couetteSolver,multiMDCellService.getMacroscopicCellService(0).getIndexConversion(),sendBuffer,globalCellIndices4SendBuffer);
      multiMDCellService.sendFromMacro2MD(sendBuffer,globalCellIndices4SendBuffer);
      std::cout << "Finish multiMDCellService.sendFromMacro2MD " << std::endl;
      // run MD instances
      for (unsigned int i = 0; i < mdInstances; i++){
        // set macroscopic cell service and interfaces in MamicoInterfaceProvider
        coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL,3>::getInstance().setMacroscopicCellService(&(multiMDCellService.getMacroscopicCellService(i)));
        coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL,3>::getInstance().setMDSolverInterface(mdSolverInterface[i]);
        
        simpleMD[i]->simulateTimesteps(simpleMDConfig.getSimulationConfiguration().getNumberOfTimesteps(),mdStepCounter);
        std::cout << "Finish simpleMD[i]->simulateTimesteps " << std::endl;
        mdStepCounter+=simpleMDConfig.getSimulationConfiguration().getNumberOfTimesteps();
        // reset mdStepCounter unless this is the last MD instance
        if (i<mdInstances-1){mdStepCounter -= simpleMDConfig.getSimulationConfiguration().getNumberOfTimesteps();}

        // plot macroscopic time step info in multi md service
        multiMDCellService.getMacroscopicCellService(i).plotEveryMacroscopicTimestep(cycles);
      }
      // send back data from MD instances and merge it
      multiMDCellService.sendFromMD2Macro(recvBuffer,globalCellIndices4RecvBuffer);
      std::cout << "Finish multiMDCellService.sendFromMD2Macro " << std::endl;

      //fillRecvBuffer(density,*couetteSolver2,multiMDCellService.getMacroscopicCellService(0).getIndexConversion(),recvBuffer,globalCellIndices4RecvBuffer, generator);


      if ( (solverType==COUETTE_LB) && twoWayCoupling && cycles == twoWayCouplingInitCycles){
        static_cast<coupling::solvers::LBCouetteSolver*>(couetteSolver)->setMDBoundary(simpleMDConfig.getDomainConfiguration().getGlobalDomainOffset(),
          simpleMDConfig.getDomainConfiguration().getGlobalDomainSize(),
          mamicoConfig.getMomentumInsertionConfiguration().getInnerOverlap());
      }

      if ( (solverType==COUETTE_LB) && twoWayCoupling && cycles >= twoWayCouplingInitCycles){
        static_cast<coupling::solvers::LBCouetteSolver*>(couetteSolver)->setMDBoundaryValues(recvBuffer,globalCellIndices4RecvBuffer,multiMDCellService.getMacroscopicCellService(0).getIndexConversion());
      }
      // write data to csv-compatible file for evaluation
      write2CSV(recvBuffer,globalCellIndices4RecvBuffer,multiMDCellService.getMacroscopicCellService(0).getIndexConversion(),rank,cycles);
      if (rank==0){std::cout << "Finish coupling cycle " << cycles << std::endl;}*/
  }

  void shutdown(){/*
    // finish time measurement for coupled simulation; start time measurement for shut down of simulation
    if(rank==0){
      gettimeofday(&end,NULL);
      double runtime = (end.tv_sec-start.tv_sec)*1000000 + (end.tv_usec-start.tv_usec);
      std::cout << "NieTest-Coupling: " << (int)(runtime/1000) << "ms for " << couplingCycles << " coupling cycles and " << totalNumberMDSimulations << " MD simulations" << std::endl;
      gettimeofday(&start,NULL);
    }

    // shutdown simulation and free buffers/arrays
    deleteBuffer(sendBuffer);
    if (globalCellIndices4SendBuffer!=NULL){delete [] globalCellIndices4SendBuffer; globalCellIndices4SendBuffer = NULL; }
    deleteBuffer(recvBuffer);
    if (globalCellIndices4RecvBuffer!=NULL){delete [] globalCellIndices4RecvBuffer; globalCellIndices4RecvBuffer = NULL; }
    for (unsigned int i = 0; i < mdInstances; i++) {
      // the shutdown operation may also delete the md solver interface; therefore, we update the MD solver interface in the vector mdSolverInteface after the shutdown is completed
      coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL,3>::getInstance().setMDSolverInterface(mdSolverInterface[i]);
      simpleMD[i]->shutdown(); delete simpleMD[i]; simpleMD[i]=NULL;
      mdSolverInterface[i] = coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL,3>::getInstance().getMDSolverInterface();
    }
    simpleMD.clear();
    for (unsigned int i = 0; i < mdInstances; i++) {if (mdSolverInterface[i]!=NULL){delete mdSolverInterface[i]; mdSolverInterface[i]=NULL;}}
    mdSolverInterface.clear();
    if (couetteSolverInterface!=NULL){delete couetteSolverInterface; couetteSolverInterface = NULL;}
    if (couetteSolver!=NULL){delete couetteSolver; couetteSolver=NULL;}

    // end time measurement for shutdown
    if(rank==0){
      gettimeofday(&end,NULL);
      double runtime = (end.tv_sec-start.tv_sec)*1000000 + (end.tv_usec-start.tv_usec);
      std::cout << "NieTest-Shutdown: " << (int)(runtime/1000) << "ms" << std::endl;
    }*/
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

  /** allocates the send buffer (with values for all macroscopic cells) and returns indices. This is only done on rank 0. */
  unsigned int* allocateSendBuffer(std::vector<coupling::datastructures::MacroscopicCell<3>* >& sendBuffer, const coupling::IndexConversion<3>& indexConversion,int rank,
  coupling::interface::MacroscopicSolverInterface<3> &couetteSolverInterface) const {
    // determine global number of cells
    const tarch::la::Vector<3,unsigned int> cells(indexConversion.getGlobalNumberMacroscopicCells()+tarch::la::Vector<3,unsigned int>(2));
    const unsigned int num = cells[0]*cells[1]*cells[2];

    // delete all potential entries of sendBuffer
    deleteBuffer(sendBuffer);

    // count number of cells to be sent from this process; therefore, loop over all global macroscopic cells...
    unsigned int numCellsSent=0;
    for (unsigned int i =0; i < num; i++){
       // ... and find out, if the current cell should be send to MD from this couette solver process
       if(couetteSolverInterface.sendMacroscopicQuantityToMDSolver(indexConversion.getGlobalVectorCellIndex(i))){
         std::vector<unsigned int> ranks = couetteSolverInterface.getSourceRanks(indexConversion.getGlobalVectorCellIndex(i));
         bool containsThisRank=false;
         for (unsigned int k = 0; k < ranks.size(); k++){
           containsThisRank = containsThisRank || (ranks[k]==(unsigned int)rank);
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
           containsThisRank = containsThisRank || (ranks[k]==(unsigned int)rank);
         }
         if (containsThisRank){
           sendBuffer.push_back(new coupling::datastructures::MacroscopicCell<3>());
           if (sendBuffer[sendBuffer.size()-1]==NULL){std::cout << "ERROR CouetteTest::allocateSendBuffer: sendBuffer[" << sendBuffer.size()-1 << "]==NULL!" << std::endl; exit(EXIT_FAILURE);}
           indices[sendBuffer.size()-1] = i;
         }
      }
    }

    #if (COUPLING_MD_DEBUG==COUPLING_MD_YES)
    for (unsigned int i = 0; i < numCellsSent; i++){
      std::vector<unsigned int> ranks = couetteSolverInterface.getSourceRanks(indexConversion.getGlobalVectorCellIndex(indices[i]));
      std::cout << "Current rank= " << rank << ", Send cell " << indexConversion.getGlobalVectorCellIndex(indices[i]) << ", ranks=";
      for (unsigned int j = 0; j < ranks.size(); j++){ std::cout << " " << ranks[j];}
      std::cout << std::endl;
    }
    #endif
    return indices;
  }

  /** allocates the recv-buffer. This buffer contains all global inner macroscopic cells, but only on rank 0. On all other ranks, no cells are stored and a NULL ptr is returned */
  unsigned int* allocateRecvBuffer(
  std::vector<coupling::datastructures::MacroscopicCell<3>* >& recvBuffer,const coupling::IndexConversion<3>& indexConversion,int rank,
  coupling::interface::MacroscopicSolverInterface<3> &couetteSolverInterface) const {

    // determine global number of cells
    const tarch::la::Vector<3,unsigned int> cells(indexConversion.getGlobalNumberMacroscopicCells()+tarch::la::Vector<3,unsigned int>(2));
    const unsigned int num = cells[0]*cells[1]*cells[2];

    // delete all potential entries of sendBuffer
    deleteBuffer(recvBuffer);

    // determine number of cells that should be received
    unsigned int numCellsRecv = 0;
    for (unsigned int i = 0; i < num; i++){
      if(couetteSolverInterface.receiveMacroscopicQuantityFromMDSolver(indexConversion.getGlobalVectorCellIndex(i))){
        std::vector<unsigned int> ranks = couetteSolverInterface.getTargetRanks(indexConversion.getGlobalVectorCellIndex(i));
        bool containsThisRank=false;
        for (unsigned int k = 0; k < ranks.size(); k++){
          containsThisRank = containsThisRank || (ranks[k]==(unsigned int)rank);
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
          containsThisRank = containsThisRank || (ranks[k]==(unsigned int)rank);
        }
        if (containsThisRank){
          recvBuffer.push_back(new coupling::datastructures::MacroscopicCell<3>());
          if (recvBuffer[recvBuffer.size()-1]==NULL){std::cout << "ERROR CouetteTest::allocateRecvBuffer: recvBuffer[" << recvBuffer.size()-1 << "]==NULL!" << std::endl; exit(EXIT_FAILURE);}
          // set linearized index
          indices[recvBuffer.size()-1] = i;
        }
      }
    }

    #if (COUPLING_MD_DEBUG==COUPLING_MD_YES)
    for (unsigned int i = 0; i < numCellsRecv; i++){
      std::vector<unsigned int> ranks = couetteSolverInterface.getTargetRanks(indexConversion.getGlobalVectorCellIndex(indices[i]));
      std::cout << "Current rank= " << rank << ", Recv cell " << indexConversion.getGlobalVectorCellIndex(indices[i]) << ", ranks=";
      for (unsigned int j = 0; j < ranks.size(); j++){ std::cout << " " << ranks[j];}
      std::cout << std::endl;
    }
    #endif
    return indices;
  }


  /** write cells that have been received from MD to csv file */
  void write2CSV(
    std::vector<coupling::datastructures::MacroscopicCell<3>* >& recvBuffer,const unsigned int * const recvIndices,
    const coupling::IndexConversion<3>& indexConversion,int rank, int couplingCycle
  ) const {
    // form file name and open file
    std::stringstream ss;
    ss << "CouetteAvgMultiMDCells_" << rank << "_" << couplingCycle << ".csv";
    std::ofstream file(ss.str().c_str());
    if (!file.is_open()){std::cout << "ERROR CouetteTest::write2CSV(): Could not open file " << ss.str() << "!" << std::endl; exit(EXIT_FAILURE);}

    // loop over received cells; read macroscopic mass+momentum buffers and write cell index, mass and velocity to one line in the csv-file
    const unsigned int numCellsRecv = recvBuffer.size();
    for (unsigned int i = 0; i < numCellsRecv; i++){
      tarch::la::Vector<3,double> vel(recvBuffer[i]->getMacroscopicMomentum());
      if (recvBuffer[i]->getMacroscopicMass()!=0.0){ vel = (1.0/recvBuffer[i]->getMacroscopicMass())*vel; }
      const tarch::la::Vector<3,unsigned int> counter(indexConversion.getGlobalVectorCellIndex(recvIndices[i]));
      file << counter[0] << " ; " << counter[1] << " ; " << counter[2] << " ; " << vel[0] << " ; " << vel[1] << " ; " << vel[2] << " ; " << recvBuffer[i]->getMacroscopicMass() << ";";
      file << std::endl;
    }

    // close file
    file.close();
  }


  /** deletes the send buffer */
  void deleteBuffer(std::vector<coupling::datastructures::MacroscopicCell<3>* >& sendBuffer) const {
    // delete all potential entries of sendBuffer
    for (unsigned int i = 0; i < sendBuffer.size(); i++){if (sendBuffer[i]!=NULL){ delete sendBuffer[i]; sendBuffer[i]=NULL;}}
    sendBuffer.clear();
  }

  /** fills send buffer with data from couette solver */
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

      // TODO check solver type or add getDensity to AbstractCouetteSolver
      double mass = density * macroscopicCellSize[0]*macroscopicCellSize[1]*macroscopicCellSize[2]
        * static_cast<const coupling::solvers::LBCouetteSolver*>(&couetteSolver)->getDensity(cellMidPoint);

      // compute momentum
      const tarch::la::Vector<3,double> momentum(mass*couetteSolver.getVelocity(cellMidPoint));
      sendBuffer[i]->setMicroscopicMass(mass);
      sendBuffer[i]->setMicroscopicMomentum(momentum);
    }
  }

  /*void fillRecvBuffer(
    const double density, const coupling::solvers::AbstractCouetteSolver<3>& couetteSolver, const coupling::IndexConversion<3>& indexConversion,
    std::vector<coupling::datastructures::MacroscopicCell<3>* >& sendBuffer, const unsigned int * const globalCellIndices4SendBuffer,
    std::default_random_engine& generator
  ) const {
    const unsigned int size = sendBuffer.size();
    const tarch::la::Vector<3,double> domainOffset(indexConversion.getGlobalMDDomainOffset());
    const tarch::la::Vector<3,double> macroscopicCellSize(indexConversion.getMacroscopicCellSize());
    const double mass = density*macroscopicCellSize[0]*macroscopicCellSize[1]*macroscopicCellSize[2];

    std::normal_distribution<double> distribution (0.0,1.0/6.0);
    for (unsigned int i = 0; i < size; i++){
      // get global cell index vector
      const tarch::la::Vector<3,unsigned int> globalIndex(indexConversion.getGlobalVectorCellIndex(globalCellIndices4SendBuffer[i]));
      // determine cell midpoint
      tarch::la::Vector<3,double> cellMidPoint(domainOffset-0.5*macroscopicCellSize);
      for (unsigned int d = 0; d < 3; d++){ cellMidPoint[d] = cellMidPoint[d] + ((double)globalIndex[d])*macroscopicCellSize[d]; }
      // compute momentum
      const tarch::la::Vector<3,double> noise(distribution(generator),distribution(generator),distribution(generator));
      const tarch::la::Vector<3,double> momentum(mass*(couetteSolver.getVelocity(cellMidPoint)+noise));
      sendBuffer[i]->setMacroscopicMass(mass);
      sendBuffer[i]->setMacroscopicMomentum(momentum);
    }
  }*/

  coupling::solvers::AbstractCouetteSolver<3>* getCouetteSolver(const double dx, const double dt){
    coupling::solvers::AbstractCouetteSolver<3>* solver = NULL;
    tarch::la::Vector<3,double> vel = _cfg.wallInitCycles > 0 ? _cfg.wallInitVelocity : _cfg.wallVelocity;
    // analytical solver: is only active on rank 0
    if (_cfg.maSolverType == COUETTE_ANALYTICAL){
      if (_rank == 0){
        solver = new coupling::solvers::CouetteSolver<3>(_cfg.channelheight,vel[0],_cfg.kinVisc);
        if (solver==NULL){
          std::cout << "ERROR CouetteTest::getCouetteSolver(): Analytic solver==NULL!" << std::endl; 
          exit(EXIT_FAILURE);
        }
      }
    // LB solver: active on lbNumberProcesses
    } else if(_cfg.maSolverType == COUETTE_LB){
      solver = new coupling::solvers::LBCouetteSolver(_cfg.channelheight,vel,_cfg.kinVisc,dx,dt,_cfg.plotEveryTimestep,"LBCouette",_cfg.lbNumberProcesses,1);
      if (solver==NULL){
        std::cout << "ERROR CouetteTest::getCouetteSolver(): LB solver==NULL!" << std::endl;
        exit(EXIT_FAILURE);
      }
    }
    else{
      std::cout << "ERROR CouetteTest::getCouetteSolver(): Unknown solver type!" << std::endl; 
      exit(EXIT_FAILURE);
    }
    return solver;
  }

  coupling::interface::MacroscopicSolverInterface<3>* getCouetteSolverInterface(
    coupling::solvers::AbstractCouetteSolver<3>* couetteSolver,
    tarch::la::Vector<3,double> mdOffset,
    tarch::la::Vector<3,double> mamicoMeshsize,
    tarch::la::Vector<3,unsigned int> globalNumberMacroscopicCells,unsigned int outerRegion
  ){
    coupling::interface::MacroscopicSolverInterface<3>* interface = NULL;
    if (_cfg.maSolverType == COUETTE_ANALYTICAL){
      interface = new coupling::solvers::CouetteSolverInterface<3>(globalNumberMacroscopicCells,outerRegion);
    } else if (_cfg.maSolverType == COUETTE_LB){
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
      interface = new coupling::solvers::LBCouetteSolverInterface(lbSolver->getAvgNumberLBCells(),lbSolver->getNumberProcesses(),offsetMDDomain,globalNumberMacroscopicCells,outerRegion);
    }
    if (interface==NULL){std::cout << "ERROR CouetteTest::getCouetteSolverInterface(...), rank=" << _rank << ": interface==NULL!" << std::endl; exit(EXIT_FAILURE);}
    return interface;
  }

  struct CouetteConfig{
    // channel is always expected to have origin at (0.0,0.0,0.0) and to be cubic (MD 30: 50.0, MD 60: 100.0, MD 120: 200.0)
    double channelheight;   
    // velocity of moving wall (lower boundary moves)
    tarch::la::Vector<3,double> wallVelocity;
    int wallInitCycles;
    tarch::la::Vector<3,double> wallInitVelocity;
    double wallOscillations;
    // number of coupling cycles, that is continuum time steps; MD/DPD: 1000
    int couplingCycles;
    bool md2Macro, macro2Md, computeSNR;   
    int filterInitCycles;
    int csvEveryTimestep;    
    MacroSolverType maSolverType;
    MicroSolverType miSolverType;
    double density, kinVisc;
    // only for LB couette solver: number of processes 
    tarch::la::Vector<3,unsigned int> lbNumberProcesses; 
    // only for LB couette solver: VTK plotting per time step
    int plotEveryTimestep;
    int initAdvanceCycles;                           
    double temp;
    int equSteps;
    int totalNumberMDSimulations;
    double noiseSigma;
    bool twsLoop;
    int twsLoopMin,twsLoopMax,twsLoopStep;
  };

  int _rank;
  simplemd::configurations::MolecularDynamicsConfiguration _simpleMDConfig;
  coupling::configurations::MaMiCoConfiguration<3> _mamicoConfig;
  CouetteConfig _cfg;
  unsigned int _mdStepCounter;
};
#endif // _COUPLING_TESTS_COUETTETEST_H_
