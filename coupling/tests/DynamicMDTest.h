// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _COUPLING_TESTS_DYNAMICMDTEST_H_
#define _COUPLING_TESTS_DYNAMICMDTEST_H_

#include "tarch/utils/MultiMDService.h"
#include "simplemd/configurations/MolecularDynamicsConfiguration.h"
#include "coupling/tests/Test.h"
#include "coupling/CouplingMDDefinitions.h"
#include "tarch/configuration/ParseConfiguration.h"
#include "tarch/utils/RandomNumberService.h"
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
#include <math.h>
#include <sys/time.h>
#include <random>

/** 
 * This is a copy of the CouetteTest.h (2020-07-17)
 * Adaptions are made to the run phase, in order to test the
 * shutdown of md simulations.
 *
 *  @author Niklas Wittmer
 */
class DynamicMDTest: public Test {
public:
  DynamicMDTest(): Test("DynamicMDTest"), _generator(0){}
  virtual ~DynamicMDTest(){}

  virtual void run(){
  #if defined(COUPLING_DYNAMIC_MD_SUDDEN )
  std::cout << "Running with MPI Shutdown test sudden." << std::endl;
  #endif
  #if defined(COUPLING_DYNAMIC_MD_SUCCESSIVE)
  std::cout << "Running with MPI Shutdown test successive." << std::endl;
  #endif
    init();

    if(_cfg.twsLoop){twsLoop();return;}
    for (int cycle = 0; cycle < _cfg.couplingCycles; cycle++) {
      runOneCouplingCycle(cycle);

#if defined(COUPLING_DYNAMIC_MD_SUDDEN)
      // Drop 50 random md instances in cycle 249
      if(cycle == 0) {
        for(int c=1;c<2;c++) {
          int iMD = c; // Global MD index to be shut down
          if(_rank == 0) std::cout << "Delete global md simulation " << iMD << std::endl;
          
          int iSim = _multiMDService->getLocalNumberOfGlobalMDSimulation(iMD);
          _multiMDCellService->rmMDSimulation(iSim, iMD);
          if(iSim >= 0 && iSim < (int)_localMDInstances) {
            _simpleMD[iSim]->shutdown();
            _simpleMD.erase(_simpleMD.begin()+iSim);
            delete _mdSolverInterface[iSim];
            _mdSolverInterface.erase(_mdSolverInterface.begin()+iSim);
            _localMDInstances -= 1;
          }
        }
      }
#endif
#if defined(COUPLING_DYNAMIC_MD_SUCCESSIVE)
      // After cycle 120 delete one md instance per cycle
      if(cycle >= 120 && cycle < 220) {
        int iMD; // Global MD index to be shut down
	iMD = 30 + cycle - 120;
        
        int iSim = _multiMDService->getLocalNumberOfGlobalMDSimulation(iMD);
        _multiMDCellService->rmMDSimulation(iSim, iMD);
        if(iSim >= 0 && iSim < (int)_localMDInstances) {
          _simpleMD[iSim]->shutdown();
          _simpleMD.erase(_simpleMD.begin()+iSim);
          _mdSolverInterface.erase(_mdSolverInterface.begin()+iSim);
          _localMDInstances -= 1;
        }
      }
#endif
    }
    shutdown();
  }

private:
  enum MacroSolverType{COUETTE_ANALYTICAL=0,COUETTE_LB=1};
  enum MicroSolverType{SIMPLEMD=0,SYNTHETIC=1};
  
  void init(){
    initMPI();
    parseConfigurations();
    initSolvers();
  }

  void runOneCouplingCycle(int cycle){
    advanceMacro(cycle);
    if(_localMDInstances > 0) advanceMicro(cycle);
    computeSNR(cycle);
    twoWayCoupling(cycle);
    // write data to csv-compatible file for evaluation
    write2CSV(_buf.recvBuffer,_buf.globalCellIndices4RecvBuffer,_multiMDCellService->getIndexConversion(),cycle);
    if(_rank==0) {std::cout << "Finish coupling cycle " << cycle << std::endl;}
  }

  void initMPI(){
    _rank = 0;
    #if (COUPLING_MD_PARALLEL==COUPLING_MD_YES)
    MPI_Comm_rank(MPI_COMM_WORLD,&_rank);
    #endif
  }

  void twsLoop(){
    for(_tws = _cfg.twsLoopMin; _tws <= _cfg.twsLoopMax; _tws += _cfg.twsLoopStep){
      init();
      for (int cycle = 0; cycle < _cfg.couplingCycles; cycle++)
        runOneCouplingCycle(cycle);
      shutdown();
    }
  }

  void parseConfigurations(){
    tarch::configuration::ParseConfiguration::parseConfiguration<simplemd::configurations::MolecularDynamicsConfiguration>("couette_simplemd.xml","molecular-dynamics",_simpleMDConfig);
    if (!_simpleMDConfig.isValid()){std::cout << "ERROR CouetteTest: Invalid SimpleMD config!" << std::endl; exit(EXIT_FAILURE);}
    
    tarch::configuration::ParseConfiguration::parseConfiguration<coupling::configurations::MaMiCoConfiguration<3> >("couette_mamico.xml","mamico",_mamicoConfig);
    if (!_mamicoConfig.isValid()){ std::cout << "ERROR CouetteTest: Invalid MaMiCo config!" << std::endl; exit(EXIT_FAILURE); }

    parseCouetteTestConfiguration();
  }

  void parseCouetteTestConfiguration(){
    tinyxml2::XMLDocument conffile;
    tinyxml2::XMLElement *node = NULL;
    conffile.LoadFile("couette.xml");
    node = conffile.FirstChildElement("couette-test");
    if (node == NULL){
      std::cout << "Could not read input file couette.xml: missing element <couette-test>" << std::endl;
      exit(EXIT_FAILURE);
    }
    tinyxml2::XMLElement *n2 = node->NextSiblingElement();
    if(n2 != NULL){
      std::cout << "Could not read input file couette.xml: unknown element " << n2->Name() << std::endl;
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
    tarch::configuration::ParseConfiguration::readBoolMandatory(_cfg.twoWayCoupling,subtag,"two-way-coupling");
    tarch::configuration::ParseConfiguration::readBoolMandatory(_cfg.md2Macro,subtag,"send-from-md-to-macro");
    tarch::configuration::ParseConfiguration::readBoolMandatory(_cfg.macro2Md,subtag,"send-from-macro-to-md");
    tarch::configuration::ParseConfiguration::readIntMandatory(_cfg.filterInitCycles,subtag,"filter-init-cycles");
    tarch::configuration::ParseConfiguration::readIntMandatory(_cfg.csvEveryTimestep,subtag,"write-csv-every-timestep");
    tarch::configuration::ParseConfiguration::readBoolMandatory(_cfg.computeSNR,subtag,"compute-snr");

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
      tarch::configuration::ParseConfiguration::readIntMandatory(_cfg.equSteps,subtag,"equilibration-steps");
      tarch::configuration::ParseConfiguration::readIntMandatory(_cfg.totalNumberMDSimulations,subtag,"number-md-simulations");
      if(_cfg.totalNumberMDSimulations < 1){
        std::cout << "Could not read input file couette.xml: number-md-simulations < 1" << std::endl;
        exit(EXIT_FAILURE);
      }
    }
    else if(type == "synthetic"){
      _cfg.miSolverType = SYNTHETIC;
      tarch::configuration::ParseConfiguration::readDoubleMandatory(_cfg.noiseSigma,subtag,"noise-sigma");
      _cfg.totalNumberMDSimulations = 1;
      tarch::configuration::ParseConfiguration::readIntOptional(_cfg.totalNumberMDSimulations,subtag,"number-md-simulations");
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
    _cfg.lbNumberProcesses = tarch::la::Vector<3,unsigned int>(1);
    tarch::configuration::ParseConfiguration::readStringMandatory(type,subtag,"type");
    if(type == "lb"){ 
      _cfg.maSolverType = COUETTE_LB;
      tarch::configuration::ParseConfiguration::readVector<3,unsigned int>(_cfg.lbNumberProcesses,subtag,"number-of-processes");
      tarch::configuration::ParseConfiguration::readIntMandatory(_cfg.plotEveryTimestep,subtag,"plot-every-timestep");
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
    tarch::configuration::ParseConfiguration::readIntMandatory(_cfg.initAdvanceCycles,subtag,"init-advance-cycles"); 

    subtag = node->FirstChildElement("tws-loop");
    if (subtag == NULL) _cfg.twsLoop = false;
    else{
      _cfg.twsLoop = true;
      tarch::configuration::ParseConfiguration::readIntMandatory(_cfg.twsLoopMin,subtag,"min");
      tarch::configuration::ParseConfiguration::readIntMandatory(_cfg.twsLoopMax,subtag,"max");
      _cfg.twsLoopStep = 1;
      tarch::configuration::ParseConfiguration::readIntOptional(_cfg.twsLoopStep,subtag,"step");
    }

    if(_cfg.miSolverType == SYNTHETIC){
      if(_cfg.md2Macro || _cfg.macro2Md || _cfg.totalNumberMDSimulations > 1 || 
        _cfg.lbNumberProcesses[0] != 1 || _cfg.lbNumberProcesses[1] != 1 || _cfg.lbNumberProcesses[2] != 1){
        std::cout << "Invalid configuration: Synthetic MD runs sequentially on rank 0 only. "
          << "It does neither support parallel communication nor multi-instance sampling" << std::endl;
        exit(EXIT_FAILURE);
      }
    }
    if(_cfg.maSolverType == COUETTE_ANALYTICAL){
      if(_cfg.twoWayCoupling){
        std::cout << "Invalid configuration: COUETTE_ANALYTICAL does not support twoWayCoupling" << std::endl;
        exit(EXIT_FAILURE);
      }
    }
  }

  void initSolvers(){
    // for timing measurements
    _tv.micro = 0;
    _tv.macro = 0;
    _tv.filter = 0;

    // allocate solvers
    _couetteSolver = NULL;
    _couetteSolver = getCouetteSolver( _mamicoConfig.getMacroscopicCellConfiguration().getMacroscopicCellSize()[0],
      _simpleMDConfig.getSimulationConfiguration().getDt()*_simpleMDConfig.getSimulationConfiguration().getNumberOfTimesteps()
      );

    // even if _cfg.miSolverType == SYNTHETIC then 
    // multiMDService, _simpleMD, _mdSolverInterface etc need to be initialized anyway,
    // so that we can finally obtain getIndexConversion from MultiMDCellService, 
    // because SYNTHETIC should fill the same cells which SIMPLEMD would use
    // depending e.g. on domain size and offset from simpleMDConfig

    _multiMDService = new tarch::utils::MultiMDService<3>(_simpleMDConfig.getMPIConfiguration().getNumberOfProcesses(), 
                                                            _cfg.totalNumberMDSimulations);
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
    if(_cfg.miSolverType == SIMPLEMD){
      // equilibrate MD
      for (unsigned int i = 0; i < _localMDInstances; i++){
        _simpleMD[i]->switchOffCoupling();
        _simpleMD[i]->simulateTimesteps(_cfg.equSteps,_mdStepCounter);
      }
      _mdStepCounter += _cfg.equSteps;
    }

    // allocate coupling interfaces
    for (unsigned int i = 0; i < _localMDInstances; i++){
      _simpleMD[i]->switchOnCoupling();
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

    if(_cfg.twsLoop){
      // initialise macroscopic cell service for multi-MD case and set single cell services in each MD simulation
      _multiMDCellService = new coupling::services::MultiMDCellService<MY_LINKEDCELL,3>(
        _mdSolverInterface,couetteSolverInterface, _simpleMDConfig.getMPIConfiguration().getNumberOfProcesses(), (unsigned int) _rank, _cfg.totalNumberMDSimulations,
        _mamicoConfig.getParticleInsertionConfiguration(), _mamicoConfig.getMomentumInsertionConfiguration(), _mamicoConfig.getBoundaryForceConfiguration(),
        _mamicoConfig.getTransferStrategyConfiguration(), _mamicoConfig.getNoiseReductionConfiguration(), _mamicoConfig.getParallelTopologyConfiguration(), _simpleMDConfig.getSimulationConfiguration().getNumberOfTimesteps(),
        _mamicoConfig.getMacroscopicCellConfiguration(), *_multiMDService, _tws
      );
    }
    else{
      // initialise macroscopic cell service for multi-MD case and set single cell services in each MD simulation
      _multiMDCellService = new coupling::services::MultiMDCellService<MY_LINKEDCELL,3>(
        _mdSolverInterface,couetteSolverInterface, _simpleMDConfig.getMPIConfiguration().getNumberOfProcesses(), (unsigned int) _rank, _cfg.totalNumberMDSimulations,
        _mamicoConfig.getParticleInsertionConfiguration(), _mamicoConfig.getMomentumInsertionConfiguration(), _mamicoConfig.getBoundaryForceConfiguration(),
        _mamicoConfig.getTransferStrategyConfiguration(), _mamicoConfig.getNoiseReductionConfiguration(), _mamicoConfig.getParallelTopologyConfiguration(), _simpleMDConfig.getSimulationConfiguration().getNumberOfTimesteps(),
        _mamicoConfig.getMacroscopicCellConfiguration(), *_multiMDService
      );
    }

    if(_cfg.miSolverType == SIMPLEMD){
      // set couette solver interface in MamicoInterfaceProvider
      coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL,3>::getInstance().setMacroscopicSolverInterface(couetteSolverInterface);

      for (unsigned int i = 0; i < _localMDInstances; i++){
        _simpleMD[i]->setMacroscopicCellService(&(_multiMDCellService->getMacroscopicCellService(i)));
        // compute and store temperature in macroscopic cells (temp=1.1 everywhere)
        _multiMDCellService->getMacroscopicCellService(i).computeAndStoreTemperature(_cfg.temp);
      }
    }

    // allocate buffers for send/recv operations
    allocateSendBuffer(_multiMDCellService->getMacroscopicCellService(0).getIndexConversion(),*couetteSolverInterface);
    allocateRecvBuffer(_multiMDCellService->getMacroscopicCellService(0).getIndexConversion(),*couetteSolverInterface);

    // manually allocate noise reduction if necessary
    if(_cfg.miSolverType == SYNTHETIC){
      if(_cfg.twsLoop){
        _noiseReduction = _mamicoConfig.getNoiseReductionConfiguration().interpreteConfiguration<3>(
          _multiMDCellService->getMacroscopicCellService(0).getIndexConversion(), *_multiMDService, _tws);
      }
      else{
        _noiseReduction = _mamicoConfig.getNoiseReductionConfiguration().interpreteConfiguration<3>(
          _multiMDCellService->getMacroscopicCellService(0).getIndexConversion(), *_multiMDService); 
      }
    }

    if(_cfg.initAdvanceCycles > 0 && _couetteSolver != NULL)
      _couetteSolver->advance(_cfg.initAdvanceCycles * _simpleMDConfig.getSimulationConfiguration().getDt()
        * _simpleMDConfig.getSimulationConfiguration().getNumberOfTimesteps());

    // finish time measurement for initialisation
    if(_rank == 0){
      gettimeofday(&_tv.end,NULL);
      double runtime = (_tv.end.tv_sec - _tv.start.tv_sec)*1000000 + (_tv.end.tv_usec - _tv.start.tv_usec);
      std::cout << "Initialization: " << (int)(runtime/1000) << "ms" << std::endl;
    }

    if(_cfg.computeSNR){
      std::cout << "Output for every coupling cycle, for the cell 87 in recvBuffer:" << std::endl;
      std::cout << "cycle number (after filter-init-cycles), vel_x macroscopic solver, vel_x filter output" << std::endl;
      _sum_signal = 0;
      _sum_noise = 0;
    }

    if (_rank == 0){ gettimeofday(&_tv.start_total,NULL); }
    std::cout << "Finish CouetteTest::initSolvers() " << std::endl;
  }

  void advanceMacro(int cycle){
    if(_couetteSolver != NULL){
      if (_rank==0){ gettimeofday(&_tv.start,NULL); }

      // run one time step for macroscopic couette solver
      if( _cfg.wallInitCycles > 0 && cycle == _cfg.wallInitCycles){
        _couetteSolver->setWallVelocity(_cfg.wallVelocity);
      }
      if(_cfg.wallOscillations != 0){
        tarch::la::Vector<3,double> vel = cycle < _cfg.wallInitCycles ? _cfg.wallInitVelocity : _cfg.wallVelocity;
        vel = vel * cos(2 * M_PI * _cfg.wallOscillations * cycle / _cfg.couplingCycles);
        _couetteSolver->setWallVelocity(vel);
      }
      _couetteSolver->advance(_simpleMDConfig.getSimulationConfiguration().getDt()*_simpleMDConfig.getSimulationConfiguration().getNumberOfTimesteps());
      if (_rank==0){
        gettimeofday(&_tv.end,NULL);
        _tv.macro += (_tv.end.tv_sec - _tv.start.tv_sec)*1000000 + (_tv.end.tv_usec - _tv.start.tv_usec);
        //std::cout << "Finish _couetteSolver->advance " << std::endl;
      }
    }

    // extract data from couette solver and send them to MD (can take any index-conversion object)
    fillSendBuffer(_cfg.density,*_couetteSolver,_multiMDCellService->getIndexConversion(),_buf.sendBuffer,_buf.globalCellIndices4SendBuffer);
    if(_cfg.macro2Md){
      _multiMDCellService->sendFromMacro2MD(_buf.sendBuffer,_buf.globalCellIndices4SendBuffer);
      //std::cout << "Finish _multiMDCellService->sendFromMacro2MD " << std::endl;
    }
    
  }

  void advanceMicro(int cycle){
    if (_rank==0){ gettimeofday(&_tv.start,NULL); }
    if(_cfg.miSolverType == SIMPLEMD){
      // run MD instances
      for (unsigned int i = 0; i < _localMDInstances; i++){
        // set macroscopic cell service and interfaces in MamicoInterfaceProvider
        coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL,3>::getInstance().setMacroscopicCellService(&(_multiMDCellService->getMacroscopicCellService(i)));
        coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL,3>::getInstance().setMDSolverInterface(_mdSolverInterface[i]);

        _simpleMD[i]->simulateTimesteps(_simpleMDConfig.getSimulationConfiguration().getNumberOfTimesteps(),_mdStepCounter);
        //std::cout << "Finish _simpleMD[i]->simulateTimesteps " << std::endl;

        // plot macroscopic time step info in multi md service
        _multiMDCellService->getMacroscopicCellService(i).plotEveryMacroscopicTimestep(cycle);
      }
      _mdStepCounter += _simpleMDConfig.getSimulationConfiguration().getNumberOfTimesteps();

      if (_rank==0){
        gettimeofday(&_tv.end,NULL);
        _tv.micro += (_tv.end.tv_sec - _tv.start.tv_sec)*1000000 + (_tv.end.tv_usec - _tv.start.tv_usec);
      }

      // send back data from MD instances and merge it
      if(_cfg.md2Macro){
        _tv.filter += _multiMDCellService->sendFromMD2Macro(_buf.recvBuffer,_buf.globalCellIndices4RecvBuffer);
        //std::cout << "Finish _multiMDCellService->sendFromMD2Macro " << std::endl;
      }
    }
    
    if(_cfg.miSolverType == SYNTHETIC){
      fillRecvBuffer(_cfg.density,*_couetteSolver,_multiMDCellService->getMacroscopicCellService(0).getIndexConversion(),_buf.recvBuffer,_buf.globalCellIndices4RecvBuffer);

      if (_rank==0){
        gettimeofday(&_tv.end,NULL);
        _tv.micro += (_tv.end.tv_sec - _tv.start.tv_sec)*1000000 + (_tv.end.tv_usec - _tv.start.tv_usec);
        gettimeofday(&_tv.start,NULL);
      }

      //call noise filter on recvBuffer
      _noiseReduction->beginProcessInnerMacroscopicCells();
      for (unsigned int i = 0; i < _buf.recvBuffer.size(); i++){
        _noiseReduction->processInnerMacroscopicCell(*_buf.recvBuffer[i],_buf.globalCellIndices4RecvBuffer[i]);
      }
      _noiseReduction->endProcessInnerMacroscopicCells();
      if(_noiseReduction->_doubleTraversal){
        _noiseReduction->beginProcessInnerMacroscopicCells();
        for (unsigned int i = 0; i < _buf.recvBuffer.size(); i++){
          _noiseReduction->processInnerMacroscopicCell(*_buf.recvBuffer[i],_buf.globalCellIndices4RecvBuffer[i]);
        }
        _noiseReduction->endProcessInnerMacroscopicCells();
      }

      if (_rank==0){
        gettimeofday(&_tv.end,NULL);
        _tv.filter += (_tv.end.tv_sec - _tv.start.tv_sec)*1000000 + (_tv.end.tv_usec - _tv.start.tv_usec);
      }
    }
  }

  void computeSNR(int cycle){
    if(_cfg.computeSNR && cycle >= _cfg.filterInitCycles){
      std::cout << cycle - _cfg.filterInitCycles << ", ";
      const coupling::IndexConversion<3>& indexConversion = _multiMDCellService->getMacroscopicCellService(0).getIndexConversion();
      const tarch::la::Vector<3,double> domainOffset(indexConversion.getGlobalMDDomainOffset());
      const tarch::la::Vector<3,double> macroscopicCellSize(indexConversion.getMacroscopicCellSize());
      const double mass = _cfg.density*macroscopicCellSize[0]*macroscopicCellSize[1]*macroscopicCellSize[2];
      for (unsigned int i = 0; i < _buf.recvBuffer.size(); i++){

        // TODO use more cells
        if (i==87){

        // get global cell index vector
        const tarch::la::Vector<3,unsigned int> globalIndex(indexConversion.getGlobalVectorCellIndex(_buf.globalCellIndices4RecvBuffer[i]));
        // determine cell midpoint
        tarch::la::Vector<3,double> cellMidPoint(domainOffset-0.5*macroscopicCellSize);
        for (unsigned int d = 0; d < 3; d++){ cellMidPoint[d] = cellMidPoint[d] + ((double)globalIndex[d])*macroscopicCellSize[d]; }
        double vx_macro = _couetteSolver->getVelocity(cellMidPoint)[0];
        double vx_filter = (1/mass * _buf.recvBuffer[i]->getMacroscopicMomentum())[0];
        _sum_noise += (vx_macro - vx_filter) * (vx_macro - vx_filter);
        _sum_signal += vx_macro * vx_macro;
        
          std::cout << vx_macro << ", " << vx_filter << std::endl;
        }
      }
    }
  }

  void twoWayCoupling(int cycle){
    if ( (_cfg.maSolverType==COUETTE_LB) && _cfg.twoWayCoupling && cycle == _cfg.filterInitCycles){
      static_cast<coupling::solvers::LBCouetteSolver*>(_couetteSolver)->setMDBoundary(_simpleMDConfig.getDomainConfiguration().getGlobalDomainOffset(),
        _simpleMDConfig.getDomainConfiguration().getGlobalDomainSize(),
        _mamicoConfig.getMomentumInsertionConfiguration().getInnerOverlap());
    }
    if ( (_cfg.maSolverType==COUETTE_LB) && _cfg.twoWayCoupling && cycle >= _cfg.filterInitCycles){
      static_cast<coupling::solvers::LBCouetteSolver*>(_couetteSolver)->setMDBoundaryValues(_buf.recvBuffer,_buf.globalCellIndices4RecvBuffer,_multiMDCellService->getMacroscopicCellService(0).getIndexConversion());
    }
    
  }

  void shutdown(){
    if(_cfg.computeSNR){
      std::cout << "SNR = " << 10 * log10(_sum_signal / _sum_noise) << std::endl;
    }

    // finish time measurement for coupled simulation
    if (_rank==0){
      gettimeofday(&_tv.end,NULL);
      double time_total = (_tv.end.tv_sec - _tv.start_total.tv_sec)*1000000 + (_tv.end.tv_usec - _tv.start_total.tv_usec);
      std::cout << "Finished all coupling cycles after " << time_total/1000000 << " s" << std::endl;
      if(_cfg.twsLoop) std::cout << "TWS = " << _tws << std::endl;
      std::cout << "Time percentages Micro, Macro, Filter, Other: " << std::endl;
      std::cout << _tv.micro/time_total*100 << ", " << _tv.macro/time_total*100 << ",  "
                << _tv.filter/time_total*100 << ", " << (1-(_tv.micro+_tv.macro+_tv.filter)/time_total)*100 << std::endl;
    }

    // free buffers/arrays
    deleteBuffer(_buf.sendBuffer);
    if (_buf.globalCellIndices4SendBuffer != NULL){
      delete [] _buf.globalCellIndices4SendBuffer; 
      _buf.globalCellIndices4SendBuffer = NULL; 
    }
    deleteBuffer(_buf.recvBuffer);
    if (_buf.globalCellIndices4RecvBuffer != NULL){
      delete [] _buf.globalCellIndices4RecvBuffer; 
      _buf.globalCellIndices4RecvBuffer = NULL; 
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
    if(_noiseReduction != NULL){delete _noiseReduction; _noiseReduction=NULL;}

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

  /** allocates the send buffer (with values for all macroscopic cells) and returns indices. This is only done on rank 0. */
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


  /** write cells that have been received from MD to csv file */
  void write2CSV(
    std::vector<coupling::datastructures::MacroscopicCell<3>* >& recvBuffer,const unsigned int * const recvIndices,
    const coupling::IndexConversion<3>& indexConversion, int couplingCycle
  ) const {
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

      double mass = density * macroscopicCellSize[0]*macroscopicCellSize[1]*macroscopicCellSize[2];
      if(_cfg.maSolverType == COUETTE_LB)
        mass *= static_cast<const coupling::solvers::LBCouetteSolver*>(&couetteSolver)->getDensity(cellMidPoint);

      // compute momentum
      const tarch::la::Vector<3,double> momentum(mass*couetteSolver.getVelocity(cellMidPoint));
      sendBuffer[i]->setMicroscopicMass(mass);
      sendBuffer[i]->setMicroscopicMomentum(momentum);
    }
  }

  void fillRecvBuffer(
    const double density, const coupling::solvers::AbstractCouetteSolver<3>& couetteSolver, const coupling::IndexConversion<3>& indexConversion,
    std::vector<coupling::datastructures::MacroscopicCell<3>* >& sendBuffer, const unsigned int * const globalCellIndices4SendBuffer
  ) {
    const unsigned int size = sendBuffer.size();
    const tarch::la::Vector<3,double> domainOffset(indexConversion.getGlobalMDDomainOffset());
    const tarch::la::Vector<3,double> macroscopicCellSize(indexConversion.getMacroscopicCellSize());
    const double mass = density*macroscopicCellSize[0]*macroscopicCellSize[1]*macroscopicCellSize[2];

    std::normal_distribution<double> distribution (0.0,_cfg.noiseSigma);
    for (unsigned int i = 0; i < size; i++){
      // get global cell index vector
      const tarch::la::Vector<3,unsigned int> globalIndex(indexConversion.getGlobalVectorCellIndex(globalCellIndices4SendBuffer[i]));
      // determine cell midpoint
      tarch::la::Vector<3,double> cellMidPoint(domainOffset-0.5*macroscopicCellSize);
      for (unsigned int d = 0; d < 3; d++){ cellMidPoint[d] = cellMidPoint[d] + ((double)globalIndex[d])*macroscopicCellSize[d]; }
      // compute momentum
      const tarch::la::Vector<3,double> noise(distribution(_generator),distribution(_generator),distribution(_generator));
      const tarch::la::Vector<3,double> momentum(mass*(couetteSolver.getVelocity(cellMidPoint)+noise));
      sendBuffer[i]->setMacroscopicMass(mass);
      sendBuffer[i]->setMacroscopicMomentum(momentum);
    }
  }

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
    bool md2Macro, macro2Md, computeSNR, twoWayCoupling;   
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

  struct CouplingBuffer{
    std::vector<coupling::datastructures::MacroscopicCell<3>* > sendBuffer;
    unsigned int *globalCellIndices4SendBuffer;
    std::vector<coupling::datastructures::MacroscopicCell<3>* > recvBuffer;
    unsigned int *globalCellIndices4RecvBuffer;
  };

  struct TimingValues{
    timeval start_total;
    timeval start;
    timeval end;
    double micro;
    double macro;
    double filter;
  };

  int _rank, _tws;
  simplemd::configurations::MolecularDynamicsConfiguration _simpleMDConfig;
  coupling::configurations::MaMiCoConfiguration<3> _mamicoConfig;
  CouetteConfig _cfg;
  unsigned int _mdStepCounter;
  coupling::solvers::AbstractCouetteSolver<3> *_couetteSolver;
  tarch::utils::MultiMDService<3>* _multiMDService;
  coupling::services::MultiMDCellService<MY_LINKEDCELL,3> *_multiMDCellService;
  CouplingBuffer _buf;
  unsigned int _localMDInstances;
  std::vector<coupling::interface::MDSolverInterface<MY_LINKEDCELL,3>* > _mdSolverInterface;
  std::vector<coupling::interface::MDSimulation*> _simpleMD;
  std::default_random_engine _generator;
  coupling::noisereduction::NoiseReduction<3>* _noiseReduction;
  double _sum_signal, _sum_noise;
  TimingValues _tv;
};
#endif // _COUPLING_TESTS_COUETTETEST_H_