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
	for(int i = 0; i < 100; i++)
		intMocks.push_back( new int(i));
	addMockService<int>(intMocks);

	//BOOL
	std::vector<bool *> boolMocks;
	boolMocks.push_back(new bool(false));
	boolMocks.push_back(new bool(true));
	addMockService<bool>(boolMocks);

	//STD::STRING
	std::vector<std::string *> stdstringMocks;
	stdstringMocks.push_back(new std::string(""));
	stdstringMocks.push_back(new std::string("Testing"));
	stdstringMocks.push_back(new std::string(" "));
	stdstringMocks.push_back(new std::string("\n"));
	stdstringMocks.push_back(new std::string("äöü"));
	addMockService<std::string>(stdstringMocks);

	//Initialize UT instances of MaMiCo classes that have Unit Tests
	_uts.push_back(new tarch::la::VectorUT<2, int>(this));
	_uts.push_back(new tarch::la::VectorUT<2, std::string>(this));
	_uts.push_back(new tarch::la::VectorUT<30, bool>(this));

	_uts.push_back(new tarch::la::VectorUT<3, tarch::la::Vector<2, int> >(this));
	//...
}


//member functions of testing::ut::UnitTestingService
void testing::ut::UnitTestingService::runAllUnitTests() {
	#ifdef DEBUG_UTS
		std::cout << "UnitTestingService: Now running tests for all classes... ";
	#endif
	for(auto ut : _uts) {
		try {
			ut->runAllTests();
		}
		catch(const std::exception& e) {
			std::cout << std::endl << "UnitTestingService: Caught " << e.what() << " while testing " << ut->getClassIdentifier() << std::endl;
			//TODO: Handle this case somehow. Aborting? Git revert? 
		}
	}
	#ifdef DEBUG_UTS
		std::cout << "Done!" << std::endl;
	#endif

}

//Assumes uts_index to not be out-of-bounds.
void testing::ut::UnitTestingService::runUnitTest(unsigned int uts_index) {
	#ifdef DEBUG_UTS
		std::cout << "UnitTestingService: Now running tests for class "<< _uts[uts_index]->getClassIdentifier() << "... ";
	#endif
	try {
		_uts[uts_index]->runAllTests();
	}
	catch(const std::exception& e){
		std::cout << std::endl << "UnitTestingService: Caught " << e.what() << " while testing " << _uts[uts_index]->getClassIdentifier() << std::endl;
		//TODO: Handle this case somehow. Aborting? Git revert? 
	}
	#ifdef DEBUG_UTS
		std::cout << "Done!" << std::endl;
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
