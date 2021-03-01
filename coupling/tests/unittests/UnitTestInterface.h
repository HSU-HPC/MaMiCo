// This file is part of the MaMiCo project. For conditions of distribution
// and use, please see the copyright notice in MaMiCo's main folder

#pragma once

#include<vector>
#include<tuple>
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
		template<class T>
		class UnitTestImpl;
}}

//TODO: A LOT OF COMMENTS

//EVERYTHING NOT TEMPLATE DEPENDANT
class testing::ut::UnitTestInterface {
	public:
		UnitTestInterface(std::string classIdentifier): _classIdentifier(classIdentifier) {}
		virtual void runAllTests() {}
		virtual void runTest() {}

		//used when handling strings
		std::string getClassIdentifier() const { return _classIdentifier; }

		//used when outputting identifer
		std::string getClassIdentifier_pretty() const { return std::string("\x1B[36m").append(_classIdentifier).append("\x1B[0m"); }

	protected:
		//prefix for all output of classes deriving this
		virtual std::string PRINT_PREFIX() const {
			return std::string("UnitTest(").append(getClassIdentifier_pretty()).append("): ");
		}

		/*
		 * Each UnitTest must provide a MockService of the class it is testing.
		 * This should be done during instanciation of the UnitTest.
		 */
		testing::ut::MockService* _mockService;

		/*
		 * This is typeid(T).name() where T is the tested class
		 */
		std::string _classIdentifier;

};

//EVERYTHING TEMPLATE DEPENDANT
template<class T>
class testing::ut::UnitTestImpl : public testing::ut::UnitTestInterface {
	public:
		UnitTestImpl():
		UnitTestInterface(typeid(T).name())
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

		//TODO: move to UTInterface?
		//Executes all unit test functions of this UT
		virtual void runAllTests() {
			for(unsigned int tf_i = 0; tf_i < _testFuncs.size(); tf_i++) 
				try {
					runTest(tf_i);
				}
				catch (...) {
					std::rethrow_exception(std::current_exception());
				}
		}

		//Executes a single test.
		virtual void runTest(unsigned int testFunc_index) {
			try {
				//TODO: time measurement
				//Gets the first part of the testFunc tuple, i.e. a function pointer, and then dereferences it
				auto testFunctionPtr = std::get<0>(_testFuncs[testFunc_index]);

				//for all possible mocks of tested type T...
				for(auto mock : this->_mockService->getMocksCopy<T>()) {
					//...execute the function
					(*testFunctionPtr)(mock);
				}

				#ifdef DEBUG_UT
					std::cout << this->PRINT_PREFIX() << std::get<1>(_testFuncs[testFunc_index]) << ": \x1B[032m Success! \x1B[0m" << std::endl;
				#endif
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

