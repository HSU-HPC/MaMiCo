// This file is part of the MaMiCo project. For conditions of distribution
// and use, please see the copyright notice in MaMiCo's main folder

#pragma once

#include<vector>
#include<tuple>
#include<iostream> 
#include<functional>
#include<exception>
#include<stdexcept>
#include<cmath>

//#define DEBUG_UT

/*
 *	Expansion of the UnitTestInterface interface.
 *	Contains all a unit test class testing T should have that depends on the type T.
 *	Split from the other non-template-dependant half to make sure all UTs share a (non-generic) type.
 *
 * @Author Leonard Hannen, Felix Maurer
 */

namespace testing { namespace ut {
		template<class T>
		class UnitTestImpl;
}}

template<class T>
class testing::ut::UnitTestImpl : public testing::ut::UnitTestInterface {
	public:
		UnitTestImpl(int rank = 0, int comm_size = 0):
		UnitTestInterface(typeid(T).name(), rank, comm_size)
		{	
			#ifdef DEBUG_UT
				std::cout << this->PRINT_PREFIX() << "Constructed." << std::endl;
			#endif
		}


		virtual ~UnitTestImpl() {
			//delete all functions pointers stored in testFuncs
			for(auto func : _testFuncs) delete std::get<0>(func);
			#ifdef DEBUG_UT
				std::cout << this->PRINT_PREFIX() << "Deconstructed." << std::endl;
			#endif
		}

		//Executes all unit test functions of this UT
		virtual void runAllTests() {

			#ifdef DEBUG_UT
			std::cout << this->PRINT_PREFIX() << ": Now running all test functions..." << std::endl;
			#endif

			for(unsigned int tf_i = 0; tf_i < _testFuncs.size(); tf_i++) 
				try {
					runTest(tf_i);
				}
				catch (...) {
					std::rethrow_exception(std::current_exception());
				}

			#ifdef DEBUG_UT
			std::cout << this->PRINT_PREFIX() << ": ...done running all test functions!" << std::endl;
			#endif
		}

		//Executes a single test.
		virtual void runTest(unsigned int testFunc_index) {
			try {
				//Gets the first part of the testFunc tuple, i.e. a function pointer, and then dereferences it
				auto testFunctionPtr = std::get<0>(_testFuncs[testFunc_index]);

				//for all possible mocks of tested type T...
				auto mocksT = this->_mockService->getMocksCopy<T>();

				//...then select subset of those mocks belonging to this rank...
				unsigned int mocksT_perRank = floor(mocksT.size()/_comm_size);
				unsigned int mocksT_startIndex = mocksT_perRank * _rank;
				unsigned int mocksT_endIndex = (_rank == _comm_size - 1) ? mocksT.size() : mocksT_startIndex+(mocksT_perRank);

				unsigned int mockCount = 0; //Counting how many mocks belong to this rank
				for(unsigned int i = mocksT_startIndex; i < mocksT_endIndex; i++) {
					//...execute the function
					(*testFunctionPtr)(mocksT[i]);
					mockCount++;
				}

				std::cout 
				<< this->PRINT_PREFIX() //usual prefix
				<< std::get<1>(_testFuncs[testFunc_index]) //function name
				<< "[Mock Count: " << mockCount //mock indices of this rank
				<< "]: \x1B[032m Success! \x1B[0m" << std::endl; //actual message

			}
			catch (std::exception& e) {
				//Enrich exception info with idenfier of failed function (i.e second part of testFunc tuple)
				throw std::runtime_error(std::string("\x1B[31m").append(std::get<1>(_testFuncs[testFunc_index])).append(": ").append(e.what()).append("\x1B[0m"));
			}
		}

	protected:	
		/*
		 * List of tuples: 
		 * (0) 	function pointers matching the following interface:
		 * 		(T&) -> void
		 *
		 * 		Expected to throw exception if assertions fail.
		 * 		That exception then gets passed to the UTS, cf. this->run(All)Tests()
		 * (1)	Identifer string for that function. Usually its source code name.
		 */
		std::vector<std::tuple<std::function<void (T&)>*, std::string>> _testFuncs;

};
