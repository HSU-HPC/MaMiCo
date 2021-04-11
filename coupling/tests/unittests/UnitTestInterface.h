// This file is part of the MaMiCo project. For conditions of distribution
// and use, please see the copyright notice in MaMiCo's main folder

#pragma once

#include<vector>
#include<iostream> 

//#define DEBUG_UT

/*
 * An abstraction of the interface of a Unit Test should to make UTs managable in data structures.
 * All Unit Tests should implement not only this, but also testing::ut::UnitTestImpl.
 *
 * @Author Leonard Hannen, Felix Maurer
 */

namespace testing { namespace ut {
		class UnitTestInterface;
}}

class testing::ut::UnitTestInterface {
	public:
		UnitTestInterface(std::string classIdentifier, int rank = 0, int comm_size = 0): 
			_classIdentifier(classIdentifier),
			_rank(rank),
			_comm_size(comm_size)
		{}

		/*
		 * These two functions are implemented by the expanded interface UnitTestImpl.
		 * All unit tests should implement not only UnitTestInterface, but also UnitTestImpl.
		 */
		virtual void runAllTests(){ throw std::runtime_error("UnitTestInterface: Tried to execute unimplemented version of runAllTests!"); }
		virtual void runTest(){ throw std::runtime_error("UnitTestInterface: Tried to execute unimplemented version of runTest!"); }

		//used when handling strings
		std::string getClassIdentifier() const { return _classIdentifier; }

		//used when outputting identifer
		std::string getClassIdentifier_pretty() const { return std::string("\x1B[36m").append(_classIdentifier).append("\x1B[0m"); }

	protected:
		//prefix for all output of classes deriving this
		virtual std::string PRINT_PREFIX() const {
			return std::string("|Rank: ").append(std::to_string(_rank)).append("| ").append("UnitTest(").append(getClassIdentifier_pretty()).append("): ");
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

		/*
		 * Used for MPI communication
		 */
		int _rank;
		int _comm_size;

};

/*
 * The class UnitTestImpl inherits from this and introduces template dependant runTest functions.
 */
#include"UnitTestImpl.h"
