// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _COUPLING_TESTS_COUETTETEST_H_
#define _COUPLING_TESTS_COUETTETEST_H_

#include "coupling/indexing/IndexingService.h"

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

//This is ignored if you dont use synthetic MD. For further instructions cf. SYNTHETIC part of initSolvers().
#define SYNTHETICMD_SEQUENCE "SYNTHETIC-MD"

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
  CouetteTest(): Test("CouetteTest"), _generator(0){}
  virtual ~CouetteTest(){}

  virtual void run(){
    init();
    if(_cfg.twsLoop){twsLoop();return;}
    for (int cycle = 0; cycle < _cfg.couplingCycles; cycle++) {
      runOneCouplingCycle(cycle);
	 }
    shutdown();
  }

private:
  enum MacroSolverType{COUETTE_ANALYTICAL=0,COUETTE_LB=1,COUETTE_FD=2, COUETTE_FOAM=3};
  enum MicroSolverType{SIMPLEMD=0,SYNTHETIC=1};

  void init(){
	using idx_scalar = coupling::indexing::CellIndex<3>;
	using idx_vector = coupling::indexing::CellIndex<3, {.vector=true}>;
	idx_scalar::setDomainParameters(); 
	idx_vector::setDomainParameters(); 
	//using idx_m2m_scalar = coupling::indexing::CellIndex<3, {.md2macro=true}>;
	using idx_m2m_vector = coupling::indexing::CellIndex<3, {.vector=true, .md2macro=true}>;
	idx_m2m_vector::setDomainParameters();

	//TODO: error. incomplete type
	//idx_global::lowerBoundary = 10;

	idx_vector idx1 = coupling::indexing::CellIndex<3, {.vector=true}>( {0,1,0} );
	idx_scalar idx2 {14};

	std::cout << "idx1: " << idx1 << std::endl;
	std::cout << "idx1 to scalar: " << (idx_scalar) idx1 << std::endl;
	std::cout << "idx1 to m2m vector: " << (idx_m2m_vector) idx1 << std::endl;
	std::cout << "idx2: " << idx2 << std::endl;
	std::cout << "idx2 to vector: " << (idx_vector) idx2 << std::endl;
	//std::cout << (idx_m2m_vector) (idx_m2m_scalar) idx1 << std::endl;

	exit(0);

    parseConfigurations();
    initSolvers();
  }

  void runOneCouplingCycle(int cycle){
    advanceMacro(cycle);
    advanceMicro(cycle);
    computeSNR(cycle);
    twoWayCoupling(cycle);
    // if(_rank==0 && cycle%20==0) {
      // std::cout << "Finish coupling cycle " << cycle << std::endl;
    // }
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
    tarch::configuration::ParseConfiguration::parseConfiguration<simplemd::configurations::MolecularDynamicsConfiguration>("couette.xml","molecular-dynamics",_simpleMDConfig);
    if (!_simpleMDConfig.isValid()){std::cout << "ERROR CouetteTest: Invalid SimpleMD config!" << std::endl; exit(EXIT_FAILURE);}
    tarch::configuration::ParseConfiguration::parseConfiguration<coupling::configurations::MaMiCoConfiguration<3> >("couette.xml","mamico",_mamicoConfig);
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
    else if(type == "fd"){
      _cfg.maSolverType = COUETTE_FD;
      tarch::configuration::ParseConfiguration::readVector<3,unsigned int>(_cfg.lbNumberProcesses,subtag,"number-of-processes");
      tarch::configuration::ParseConfiguration::readIntMandatory(_cfg.plotEveryTimestep,subtag,"plot-every-timestep");
    }
    #if(BUILD_WITH_OPENFOAM)
    else if(type == "foam"){
      _cfg.maSolverType = COUETTE_FOAM;
      tarch::configuration::ParseConfiguration::readIntMandatory(_cfg.plotEveryTimestep,subtag,"plot-every-timestep");
      tarch::configuration::ParseConfiguration::readStringMandatory(_foam.directory,subtag,"foam-setup-directory");
      tarch::configuration::ParseConfiguration::readStringMandatory(_foam.folder,subtag,"foam-setup-folder");
      tarch::configuration::ParseConfiguration::readVector<12,unsigned int>(_foam.boundariesWithMD,subtag,"boundaries-with-MD");
    }
    #endif
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
      if(/*_cfg.md2Macro ||*/ _cfg.macro2Md || _cfg.totalNumberMDSimulations > 1 ||
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
	if (_couetteSolver != NULL) std::cout << "Couette solver not null on rank: " << _rank << std::endl; //TODO: remove debug

    // even if _cfg.miSolverType == SYNTHETIC then
    // multiMDService, _simpleMD, _mdSolverInterface etc need to be initialized anyway,
    // so that we can finally obtain getIndexConversion from MultiMDCellService,
    // because SYNTHETIC should fill the same cells which SIMPLEMD would use
    // depending e.g. on domain size and offset from simpleMDConfig

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
        _mamicoConfig.getTransferStrategyConfiguration(), _mamicoConfig.getParallelTopologyConfiguration(), _simpleMDConfig.getSimulationConfiguration().getNumberOfTimesteps(),
        _mamicoConfig.getMacroscopicCellConfiguration(), "couette.xml", *_multiMDService, _tws
      );
    }
    else{
      // initialise macroscopic cell service for multi-MD case and set single cell services in each MD simulation
      _multiMDCellService = new coupling::services::MultiMDCellService<MY_LINKEDCELL,3>(
        _mdSolverInterface,couetteSolverInterface, _simpleMDConfig.getMPIConfiguration().getNumberOfProcesses(), (unsigned int) _rank, _cfg.totalNumberMDSimulations,
        _mamicoConfig.getParticleInsertionConfiguration(), _mamicoConfig.getMomentumInsertionConfiguration(), _mamicoConfig.getBoundaryForceConfiguration(),
        _mamicoConfig.getTransferStrategyConfiguration(), _mamicoConfig.getParallelTopologyConfiguration(), _simpleMDConfig.getSimulationConfiguration().getNumberOfTimesteps(),
        _mamicoConfig.getMacroscopicCellConfiguration(), "couette.xml", *_multiMDService
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

	/*
	 * A synthethic solver is modeled using a dynamically linked filter,
	 * i.e. a lambda function producing artifical data in every filter step.
	 *
	 * This is how to properly instanciate and use a synthethic solver:
	 * - Create a sequence named like the macro SYNTHETICMD_SEQUENCE in xml. Use whatever input, but make sure the filter system's output is set to this sequence.
	 * - Set filtered-values = "macro-mass macro-momentum" for that sequence.
	 * - Use that sequence as input for all sequences that want (unfiltered) MD input.
	 *
	 * Note that this synthethic solver is designed to be used on sequential mode only and with only one MD instance.
	 *
	 * TODO
	 * - major bug when there is ONLY a FFF in a sequence
	 * - reduce capture: most variables in lambda can be defined beforehand as they are const (e.g. everything coming from cfg)
	 * - totalNumberMDSimulations > 1 is theoretically possible with this redesign. test it and remove the restriction for it to be 1
	 */

	else if(_cfg.miSolverType == SYNTHETIC)	{
		
		/*
		 * Synthetic MD runs sequentially only, as described above.
		 */
        if(/*_cfg.md2Macro ||*/ _cfg.macro2Md || _cfg.totalNumberMDSimulations > 1 ||
        _cfg.lbNumberProcesses[0] != 1 || _cfg.lbNumberProcesses[1] != 1 || _cfg.lbNumberProcesses[2] != 1) {
			throw std::runtime_error("ERROR: Syntethic MD is only available in sequential mode!");
		}

		/*
		 * Create new FilterFromFunction instance and insert it into Filtering System.
		 */
		try {
			for(unsigned int i = 0; i < _localMDInstances; i++) 
				_multiMDCellService->getMacroscopicCellService(i).getFilterPipeline().getSequence(SYNTHETICMD_SEQUENCE)->addFilter(
				new std::function<std::vector<double> (std::vector<double>, std::vector<std::array<unsigned int, 3>>)>{ //applyScalar
				[this] (
   	 				std::vector<double> inputScalars, //doesnt get used: matching MCS's addFilterToSequence(...) signature
					std::vector<std::array<unsigned int, 3>> cellIndices //only gets used to determine "size" (see below)
  				) {
					if (_rank==0){ gettimeofday(&_tv.start,NULL); }

					//sanity check
					if(inputScalars.size() != cellIndices.size())
						throw std::runtime_error("ERROR: Cell data and indexing of non-matching shapes!");

					//std::cout << "Entering synthetic MD scalar..." << std::endl;
					const coupling::IndexConversion<3>& indexConversion = _multiMDCellService->getMacroscopicCellService(0).getIndexConversion();
					const unsigned int size = cellIndices.size();
					const tarch::la::Vector<3,double> macroscopicCellSize(indexConversion.getMacroscopicCellSize());
					const double mass = (_cfg.density)*macroscopicCellSize[0]*macroscopicCellSize[1]*macroscopicCellSize[2];

					std::vector<double> syntheticMasses;
					for (unsigned int i = 0; i < size; i++){
						syntheticMasses.push_back(mass);
					}
					//std::cout << "Generated masses!" << std::endl;
					
					if (_rank==0){
						gettimeofday(&_tv.end,NULL);
						_tv.micro += (_tv.end.tv_sec - _tv.start.tv_sec)*1000000 + (_tv.end.tv_usec - _tv.start.tv_usec);
					}

					return syntheticMasses;

  				}},
				new std::function<std::vector<std::array<double, 3>> (std::vector<std::array<double,3>>, std::vector<std::array<unsigned int, 3>>)> { //applyVector
				[this] (	
					std::vector<std::array<double, 3>> inputVectors, //same for these 2
					std::vector<std::array<unsigned int, 3>> cellIndices
  				) {
					if (_rank==0){ gettimeofday(&_tv.start,NULL); }

					//sanity check
					if(inputVectors.size() != cellIndices.size())
						throw std::runtime_error("ERROR: Cell data and indexing of non-matching shapes!");

					//std::cout << "Entering synthetic MD vector." << std::endl;

					//unlike for the scalar case, we need the MD2Macro version of IC to calculate correct offsets
					const coupling::IndexConversionMD2Macro<3>* indexConversionMD2Macro =
						_multiMDCellService->getMacroscopicCellService(0).getFilterPipeline().getICM2M();

					const unsigned int size = cellIndices.size();
					const tarch::la::Vector<3,double> md2MacroDomainOffset = indexConversionMD2Macro->getGlobalMD2MacroDomainOffset();

					const tarch::la::Vector<3,double> macroscopicCellSize(indexConversionMD2Macro->getBaseIC()->getMacroscopicCellSize());
					const double mass = (_cfg.density)*macroscopicCellSize[0]*macroscopicCellSize[1]*macroscopicCellSize[2];

					std::normal_distribution<double> distribution (0.0,_cfg.noiseSigma);
					std::vector<std::array<double, 3>> syntheticMomenta;
					for (unsigned int i = 0; i < size; i++){
						// determine cell midpoint
						const tarch::la::Vector<3,unsigned int> globalIndex(indexConversionMD2Macro->getBaseIC()->getGlobalVectorCellIndex(_buf.globalCellIndices4RecvBuffer[i]));
						tarch::la::Vector<3,double> cellMidPoint(md2MacroDomainOffset-0.5*macroscopicCellSize);
						for (unsigned int d = 0; d < 3; d++){ cellMidPoint[d] = cellMidPoint[d] + ((double)globalIndex[d])*macroscopicCellSize[d]; }

						//compute momentum
						const tarch::la::Vector<3,double> noise(distribution(_generator),distribution(_generator),distribution(_generator));
						const tarch::la::Vector<3,double> momentum(mass*((*_couetteSolver).getVelocity(cellMidPoint)+noise));

						//conversion from tarch::la::Vector to std::array
						syntheticMomenta.push_back({momentum[0], momentum[1], momentum[2]});
					}
					//std::cout << "Generated momenta!" << std::endl;
					
					if (_rank==0){
						gettimeofday(&_tv.end,NULL);
						_tv.micro += (_tv.end.tv_sec - _tv.start.tv_sec)*1000000 + (_tv.end.tv_usec - _tv.start.tv_usec);
					}

					return syntheticMomenta;
  				}},
			   	0 //filterIndex
				);
		}
		catch(std::runtime_error& e) {
			auto expectedError = std::string("ERROR: Could not find Filter Sequence named ").append(SYNTHETICMD_SEQUENCE);
			if(expectedError.compare(e.what()) == 0)  {
				std::cout << "ERROR: Synthetic MD solver selected without providing filter sequence '" << SYNTHETICMD_SEQUENCE <<"' in config." << std::endl;
				exit(EXIT_FAILURE);
			}
			else throw e;
		}
	}

    // allocate buffers for send/recv operations
    allocateSendBuffer(_multiMDCellService->getMacroscopicCellService(0).getIndexConversion(),*couetteSolverInterface);
    allocateRecvBuffer(_multiMDCellService->getMacroscopicCellService(0).getIndexConversion(),*couetteSolverInterface);

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
		//When using Synthetic MD, 
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
    fillSendBuffer(_cfg.density,*_couetteSolver,_multiMDCellService->getMacroscopicCellService(0).getIndexConversion(),_buf.sendBuffer,_buf.globalCellIndices4SendBuffer);
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

	//Now that synthetic filters are modelled as filters, this looks very similar to the SimpleMD case above...
    if(_cfg.miSolverType == SYNTHETIC){
      //fillRecvBuffer(_cfg.density,*_couetteSolver,_multiMDCellService->getMacroscopicCellService(0).getIndexConversion(),_buf.recvBuffer,_buf.globalCellIndices4RecvBuffer);
	  for (unsigned int i = 0; i < _localMDInstances; i++){
        // set macroscopic cell service and interfaces in MamicoInterfaceProvider
        coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL,3>::getInstance().setMacroscopicCellService(&(_multiMDCellService->getMacroscopicCellService(i)));
        coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL,3>::getInstance().setMDSolverInterface(_mdSolverInterface[i]);
      }

      // send back data from MD instances and merge it
      if(_cfg.md2Macro){
		//_buf does not get used here: Instead, the synthetic MD in the SYNTHETICMD_SEQUENCE generates values. To prevent segfaults, it has to be nonempty, though.
		/*std::cout << "RECV BUF ADRESSES:" << std::endl;
		for (auto cell : _buf.recvBuffer) {
			std::cout << cell << ", ";
		}
		std::cout << std::endl;*/
        _tv.filter += _multiMDCellService->sendFromMD2Macro(_buf.recvBuffer,_buf.globalCellIndices4RecvBuffer);
        //std::cout << "Finish _multiMDCellService->sendFromMD2Macro " << std::endl;
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
    if (( _cfg.maSolverType==COUETTE_LB || _cfg.maSolverType==COUETTE_FD) && _cfg.twoWayCoupling && cycle == _cfg.filterInitCycles){
      static_cast<coupling::solvers::LBCouetteSolver*>(_couetteSolver)->setMDBoundary(_simpleMDConfig.getDomainConfiguration().getGlobalDomainOffset(),
        _simpleMDConfig.getDomainConfiguration().getGlobalDomainSize(),
        _mamicoConfig.getMomentumInsertionConfiguration().getInnerOverlap(),
        _multiMDCellService->getMacroscopicCellService(0).getIndexConversion(),  _buf.globalCellIndices4RecvBuffer, _buf.recvBuffer.size());
    }
    #if(BUILD_WITH_OPENFOAM)
    else if ( (_cfg.maSolverType==COUETTE_FOAM) && _cfg.twoWayCoupling && cycle == _cfg.filterInitCycles){
      static_cast<coupling::solvers::IcoFoam*>(_couetteSolver)->setMDBoundary(_simpleMDConfig.getDomainConfiguration().getGlobalDomainOffset(),
      _simpleMDConfig.getDomainConfiguration().getGlobalDomainSize(), _mamicoConfig.getMomentumInsertionConfiguration().getInnerOverlap(),
      _multiMDCellService->getMacroscopicCellService(0).getIndexConversion(),  _buf.globalCellIndices4RecvBuffer, _buf.recvBuffer.size());
    }
    #endif
    if (( _cfg.maSolverType==COUETTE_LB || _cfg.maSolverType==COUETTE_FD) && _cfg.twoWayCoupling && cycle >= _cfg.filterInitCycles){
      static_cast<coupling::solvers::LBCouetteSolver*>(_couetteSolver)->setMDBoundaryValues(_buf.recvBuffer,_buf.globalCellIndices4RecvBuffer,_multiMDCellService->getMacroscopicCellService(0).getIndexConversion());
    }
    #if(BUILD_WITH_OPENFOAM)
    else if (_cfg.maSolverType==COUETTE_FOAM && _cfg.twoWayCoupling && cycle >= _cfg.filterInitCycles){
      static_cast<coupling::solvers::IcoFoam*>(_couetteSolver)->setMDBoundaryValues(_buf.recvBuffer,_buf.globalCellIndices4RecvBuffer,_multiMDCellService->getMacroscopicCellService(0).getIndexConversion());
    }
    #endif
    // write data to csv-compatible file for evaluation
    write2CSV(_buf.recvBuffer,_buf.globalCellIndices4RecvBuffer,_multiMDCellService->getMacroscopicCellService(0).getIndexConversion(),cycle+1);
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
      if(_cfg.maSolverType == COUETTE_LB || _cfg.maSolverType == COUETTE_FD)
        mass *= static_cast<const coupling::solvers::LBCouetteSolver*>(&couetteSolver)->getDensity(cellMidPoint);

      // compute momentum
      tarch::la::Vector<3,double> momentum(mass*couetteSolver.getVelocity(cellMidPoint));
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
	//In case of synthetic MD, each rank needs access to the analytical solution. In that case we thus initialize analytical solvers on each rank other than 0
	else if( _cfg.miSolverType == SYNTHETIC){ //rank != 0
      if (_rank == 0){
        solver = new coupling::solvers::CouetteSolver<3>(_cfg.channelheight,vel[0],_cfg.kinVisc);
        if (solver==NULL){ //How is this even reachable? Copied it from above...
          std::cout << "ERROR CouetteTest::getCouetteSolver(): Analytic solver==NULL!" << std::endl;
          exit(EXIT_FAILURE);
        }
      }


	  }
    }
    #if(BUILD_WITH_OPENFOAM)
    else if(_cfg.maSolverType == COUETTE_FOAM){
      solver = new coupling::solvers::IcoFoam(_rank, _cfg.plotEveryTimestep, _cfg.channelheight, _foam.directory, _foam.folder, _foam.boundariesWithMD);
      if (solver==NULL){
        std::cout << "ERROR CouetteTest::getCouetteSolver(): IcoFoam solver==NULL!" << std::endl;
        exit(EXIT_FAILURE);
      }
    }
    #endif
    // LB solver: active on lbNumberProcesses
    else if(_cfg.maSolverType == COUETTE_LB){
      solver = new coupling::solvers::LBCouetteSolver(_cfg.channelheight,vel,_cfg.kinVisc,dx,dt,_cfg.plotEveryTimestep,"LBCouette",_cfg.lbNumberProcesses,1);
      if (solver==NULL){
        std::cout << "ERROR CouetteTest::getCouetteSolver(): LB solver==NULL!" << std::endl;
        exit(EXIT_FAILURE);
      }
    }
    else if(_cfg.maSolverType == COUETTE_FD){
      solver = new coupling::solvers::FiniteDifferenceSolver(_cfg.channelheight,vel,_cfg.kinVisc,dx,dt,_cfg.plotEveryTimestep,"FDCouette",_cfg.lbNumberProcesses,1);
      if (solver==NULL){
        std::cout << "ERROR CouetteTest::getCouetteSolver(): FD solver==NULL!" << std::endl;
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
    #if(BUILD_WITH_OPENFOAM)
    else if (_cfg.maSolverType == COUETTE_FOAM){
      interface = new coupling::solvers::FoamSolverInterface<3>(globalNumberMacroscopicCells,outerRegion);
    }
    #endif
    else if (_cfg.maSolverType == COUETTE_FD){
      coupling::solvers::FiniteDifferenceSolver *fdSolver = static_cast<coupling::solvers::FiniteDifferenceSolver*>(couetteSolver);
      if (fdSolver==NULL){std::cout << "ERROR CouetteTest::getCouetteSolverInterface(...), rank=" << _rank << ": Could not convert abstract to LB solver!" << std::endl; exit(EXIT_FAILURE);}
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
      interface = new coupling::solvers::LBCouetteSolverInterface(fdSolver->getAvgNumberLBCells(),fdSolver->getNumberProcesses(),offsetMDDomain,globalNumberMacroscopicCells,outerRegion);
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
  #if(BUILD_WITH_OPENFOAM)
  struct FoamConfig{
    std::string directory;
    std::string folder;
    tarch::la::Vector<12,unsigned int> boundariesWithMD;
  };
  #endif

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
  double _sum_signal, _sum_noise;
  TimingValues _tv;
  #if(BUILD_WITH_OPENFOAM)
  FoamConfig _foam;
  #endif
};

#endif // _COUPLING_TESTS_COUETTETEST_H_
