// This file is part of the MaMiCo project. For conditions of distribution
// and use, please see the copyright notice in MaMiCo's main folder

#pragma once

#include<vector>
#include<iostream> 
#include<functional>
#include<exception>
#include<stdexcept>

#define DEBUG_UT

/*
 *TODO: explanatory interface comment
 *
 * @Author Leonard Hannen, Felix Maurer
 */

namespace testing { namespace ut {
		class UnitTestInterface;
}}

class testing::ut::UnitTestInterface {
	public:
		UnitTestInterface(std::string classIdentifier):
		_classIdentifier(classIdentifier)
		{	
			#ifdef DEBUG_UT
				std::cout << PRINT_PREFIX() << "Constructed." << std::endl;
			#endif
		}


		virtual ~UnitTestInterface() {
			for(auto func : _testFuncs) delete func;
			#ifdef DEBUG_UT
				std::cout << PRINT_PREFIX() << "Deconstructed." << std::endl;
			#endif
		}

		//Executes all unit test functions of this UT
		virtual void runAllTests() {
			for(auto tf : _testFuncs) 
				try {
					(*tf)();
				}
				catch (...) {
					std::rethrow_exception(std::current_exception());
				}
		}

		//Executes a single test.
		virtual void runTest(unsigned int testFunc_index) {
			try {
				(*_testFuncs[testFunc_index])();
			}
			catch (...) {
				std::rethrow_exception(std::current_exception());
			}
		}

		std::string getClassIdentifier() const { return _classIdentifier; }

	protected:

		virtual std::string PRINT_PREFIX() const {
			return std::string("UnitTest( ").std::string::append(_classIdentifier).std::string::append(" ): ");
		}

		/*
		 * List of function pointers matching the following interface:
		 * (void) -> void
		 *
		 * Expected to throw exception if assertions fail.
		 * That exception then gets passed to the UTS, cf. this->run(All)Tests()
		 */
		std::vector<std::function<void ()>*> _testFuncs;

		/*
		 * Each UnitTest must provide a MockService of the class it is testing.
		 * This should be done during instanciation of the UnitTest.
		 */
		testing::ut::MockService* _thisMs;

		/*
		 * This should be typeid(T).name() where T is the tested class
		 */
		std::string _classIdentifier;
};

