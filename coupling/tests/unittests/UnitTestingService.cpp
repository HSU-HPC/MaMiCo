// This file is part of the MaMiCo project. For conditions of distribution
// and use, please see the copyright notice in MaMiCo's main folder

//constructors of testing::ut::UnitTestingService
testing::ut::UnitTestingService::UnitTestingService()
{
	//init MS instances for non-Mamico (e.g. built-in or STL) types
	//TODO: other primitives
	//TODO: confused about const-correctness. should double check
	
	//INT
	std::vector<int *> intMocks;
	for(int i = -1000; i < 1000; i++)
		intMocks.push_back( new int(i));
	addMockService<int>(intMocks);

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

	//Initialize UT instances of MaMiCo classes that have Unit Tests
	_uts.push_back(new tarch::la::VectorUT<2, int>(this));
	_uts.push_back(new tarch::la::VectorUT<2, std::string>(this));
	_uts.push_back(new tarch::la::VectorUT<5, bool>(this));
	_uts.push_back(new tarch::la::VectorUT<3, tarch::la::Vector<2, int> >(this));
	_uts.push_back(new tarch::la::VectorUT<2, double>(this));
	_uts.push_back(new tarch::la::VectorUT<3, double>(this));

	_uts.push_back(new coupling::datastructures::MacroscopicCellUT<2>(this));
	_uts.push_back(new coupling::datastructures::MacroscopicCellUT<3>(this));
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
		std::cout << "\x1B[1mUnitTestingService:\x1B[0m Caught error in function " << e.what() << " while testing " << _uts[uts_index]->getClassIdentifier_pretty() << std::endl;
		//TODO: Handle this case somehow. Aborting? Git revert? 
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
