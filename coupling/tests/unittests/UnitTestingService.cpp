// This file is part of the MaMiCo project. For conditions of distribution
// and use, please see the copyright notice in MaMiCo's main folder

//constructors of testing::ut::UnitTestingService
testing::ut::UnitTestingService::UnitTestingService(MPI_Comm comm): _comm(comm) {

	MPI_Comm_size(comm, &_comm_size);
	MPI_Comm_rank(comm, &_rank);

	//Initialize MS instances for non-Mamico (e.g. built-in or STL) types
	//TODO: other primitives, STL types like vector
	
	//INT
	std::vector<int *> intMocks;
	for(int i = -1000; i < 1000; i++)
		intMocks.push_back( new int(i));
	addMockService<int>(intMocks);

	//UNSIGNED INT
	std::vector<unsigned int *> uintMocks;
	for(int i = 0; i < 1000; i++)
		uintMocks.push_back( new unsigned int(i));
	addMockService<unsigned int>(uintMocks);


	//BOOL
	std::vector<bool *> boolMocks;
	boolMocks.push_back(new bool(false));
	boolMocks.push_back(new bool(true));
	addMockService<bool>(boolMocks);

	//DOUBLE
	std::vector<double *> doubleMocks;
	for(double d = -500.0; d < 500.0; d += 0.1)
		doubleMocks.push_back( new double(d));
	addMockService<double>(doubleMocks);
	

	//STD::STRING
	std::vector<std::string *> stdstringMocks;
	stdstringMocks.push_back(new std::string(""));
	stdstringMocks.push_back(new std::string("Testing"));
	stdstringMocks.push_back(new std::string(" "));
	stdstringMocks.push_back(new std::string("\u03BC"));
	stdstringMocks.push_back(new std::string("äöü"));
	addMockService<std::string>(stdstringMocks);

	
		
	//Initialize Instances of (Simple)MD (& MD interface & IndexConversion & Mamico-Config)
	simplemd::configurations::MolecularDynamicsConfiguration simpleMDConfig;
	//Add all .xmls to this initializer list
	for(auto xml : { "test.xml" }) {

		//Init SimpleMD config. We won't use this later, so we don't allocate it on heap.
		tarch::configuration::ParseConfiguration::parseConfiguration<simplemd::configurations::MolecularDynamicsConfiguration>("test.xml","molecular-dynamics",simpleMDConfig);
    	if (!simpleMDConfig.isValid()) {
			std::cout << "ERROR UnitTesting: Invalid SimpleMD config in: " << xml << std::endl;
			exit(EXIT_FAILURE);
		}

		//Init new MamicoConfiguration 
		auto newMamicoConfig = new coupling::configurations::MaMiCoConfiguration<3>(); //TODO: other dims
		tarch::configuration::ParseConfiguration::parseConfiguration<coupling::configurations::MaMiCoConfiguration<3> >("test.xml","mamico",*newMamicoConfig);
    	if (!newMamicoConfig->isValid()) {
			std::cout << "ERROR UnitTesting: Invalid Mamico config (for SimpleMD) in: " << xml << std::endl; 
			delete newMamicoConfig;
			exit(EXIT_FAILURE);
		}

		_mamicoConfigs.push_back(newMamicoConfig);


		//Init new SimpleMD
		auto newSimpleMDInstance = coupling::interface::SimulationAndInterfaceFactory::getInstance().getMDSimulation(
        	simpleMDConfig,*newMamicoConfig, MPI_COMM_WORLD
      	);
		newSimpleMDInstance->init();

		//Init new MDSolverInterface
		_simpleMDs.push_back(newSimpleMDInstance);

		auto newMDInterface = coupling::interface::SimulationAndInterfaceFactory::getInstance().
        getMDSolverInterface(simpleMDConfig, *newMamicoConfig, newSimpleMDInstance);

		_mdSolverInterfaces.push_back(newMDInterface);

      	if (_simpleMDs.back() == nullptr or _mdSolverInterfaces.back() == nullptr){
			//TODO: More verbose error message
        	std::cout << "ERROR UnitTesting: SimpleMD or MDSolverInterface factory returned nullptr." << std::endl;
        	exit(EXIT_FAILURE);
      	}

		//Init new IndexConversion. This part is largely inspired by MacroscopicCellService::initIndexConversion
		auto globalMDDomainSize = newMDInterface->getGlobalMDDomainSize(); 
		auto globalMDDomainOffset = newMDInterface->getGlobalMDDomainOffset();
		auto macroscopicCellSize = newMamicoConfig->getMacroscopicCellConfiguration().getMacroscopicCellSize();
		auto parallelTopologyType = newMamicoConfig->getParallelTopologyConfiguration().getParallelTopologyType(); 

		tarch::la::Vector<3,unsigned int> globalNumberMacroscopicCells(0);
  		for (unsigned int d = 0; d < 3; d++) {
    		globalNumberMacroscopicCells[d] = (unsigned int) floor( globalMDDomainSize[d]/macroscopicCellSize[d] + 0.5 );
    		if ( fabs(globalNumberMacroscopicCells[d]*macroscopicCellSize[d] - globalMDDomainSize[d]) > 1e-13 ){
      			std::cout << "WARNING UnitTesting: Deviation of domain size MD vs macro cells > 1e-13!" << std::endl;
    		}
  		}

  		_indexConversions.push_back(new coupling::IndexConversion<3>(globalNumberMacroscopicCells,1,1,globalMDDomainSize,globalMDDomainOffset,parallelTopologyType));
	}
	

	//Initialize instances of CS
	//TODO: Currently not required by any Unit Test already implemented. Will need this in the future.

	//Initialize UT instances of MaMiCo classes that have Unit Tests
	_uts.push_back(new tarch::la::VectorUT<2, int>(this));
	_uts.push_back(new tarch::la::VectorUT<2, unsigned int>(this));
	_uts.push_back(new tarch::la::VectorUT<3, unsigned int>(this));
	_uts.push_back(new tarch::la::VectorUT<2, std::string>(this));
	_uts.push_back(new tarch::la::VectorUT<5, bool>(this));
	_uts.push_back(new tarch::la::VectorUT<3, tarch::la::Vector<2, int> >(this));
	_uts.push_back(new tarch::la::VectorUT<2, double>(this));
	_uts.push_back(new tarch::la::VectorUT<3, double>(this));

	_uts.push_back(new coupling::datastructures::MacroscopicCellUT<2>(this));
	_uts.push_back(new coupling::datastructures::MacroscopicCellUT<3>(this));

	_uts.push_back(new coupling::datastructures::MacroscopicCellWithLinkedCellsUT<MY_LINKEDCELL,2>(this));
	_uts.push_back(new coupling::datastructures::MacroscopicCellWithLinkedCellsUT<MY_LINKEDCELL,3>(this));

	_uts.push_back(new coupling::paralleltopology::XYZTopologyUT<2>(this));
	_uts.push_back(new coupling::paralleltopology::XYZTopologyUT<3>(this));

	_uts.push_back(new coupling::cellmappings::ComputeMassMappingUT<MY_LINKEDCELL,3>(this));

	//Note: Creating Usher UTs with other template parameters is currently not supported.
	_uts.push_back(new coupling::UsherParticleInsertionUT<MY_LINKEDCELL,3>(this));

	//...
	//
	#ifdef DEBUG_UTS
		std::cout << "\x1B[1mUnitTestingService:\x1B[0m Constructed." << std::endl;
	#endif
}


//member functions of testing::ut::UnitTestingService
void testing::ut::UnitTestingService::runAllUnitTests() {
	#ifdef DEBUG_UTS
		std::cout << "\x1B[1mUnitTestingService:\x1B[0m Now running tests for all classes... " << std::endl;
	#endif
	for(unsigned int uts_index = 0; uts_index < _uts.size() ;uts_index++) {
		runUnitTest(uts_index);
	}
	#ifdef DEBUG_UTS
		std::cout << "Done testing all classes!" << std::endl;
	#endif

}

//Assumes uts_index to not be out-of-bounds.
void testing::ut::UnitTestingService::runUnitTest(unsigned int uts_index) {
	#ifdef DEBUG_UTS
		std::cout << "\x1B[1mUnitTestingService:\x1B[0m Now running tests for class "<< _uts[uts_index]->getClassIdentifier_pretty() << "... " << std::endl;
	#endif
	try {
		_uts[uts_index]->runAllTests();
	}
	catch(const std::exception& e){
		std::cout << "|Rank: " << _rank << "| \x1B[1mUnitTestingService:\x1B[0m Caught error in function " << e.what() << " while testing " << _uts[uts_index]->getClassIdentifier_pretty() << std::endl;
	}
	#ifdef DEBUG_UTS
		std::cout << "Done running all unit tests for class " << _uts[uts_index]->getClassIdentifier_pretty() << "!" << std::endl;
	#endif

}

template<class T>
std::optional<testing::ut::MockService *> testing::ut::UnitTestingService::getMockService() {
	//requires c++20
	if ( !_mockServices.contains( typeid(T).name() )) 
		return std::nullopt;
	else 
		//create optional monad from mockService pointer
		return std::optional<testing::ut::MockService *>(_mockServices[typeid(T).name()]);
}

template<class T>
testing::ut::MockService* testing::ut::UnitTestingService::addMockService(std::vector<T *> mockValues){
	auto msT = new testing::ut::MockService();
	for (auto mock : mockValues) msT->addMock(mock);
	_mockServices[typeid(T).name()] = msT;
	return msT;
}
