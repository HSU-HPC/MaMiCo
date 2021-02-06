// This file is part of the MaMiCo project. For conditions of distribution
// and use, please see the copyright notice in MaMiCo's main folder

#pragma once

#include<vector>
#include<map>
#include<any>
#include<optional>
#include<iostream> 
#include"MockService.cpph"
#include"UnitTestInterface.h"

#define DEBUG_UTS

/*
 *TODO: explanatory interface comment
 *
 * @Author Leonard Hannen, Felix Maurer
 */

namespace testing { namespace ut {
		class UnitTestingService;
}}


//Make this static?
class testing::ut::UnitTestingService {
	public:
		UnitTestingService();

		~UnitTestingService() {
			for (auto ut : _uts) delete ut;

			for (auto ms_tuple : _mockServices)delete ms_tuple.second;

			#ifdef DEBUG_UTS
				std::cout << "UnitTestingService: Deconstructed." << std::endl;
			#endif
		}

		//Executes all unit test functions of all UTs in _uts
		void runAllUnitTests();

		//Executes all unit test functions of a single UT in _uts.
		//Specified via index in _uts.
		void runUnitTest(unsigned int uts_index); 

		//Returns std::nullopt if type T is not found in _mockServices
		template<class T>
		std::optional<testing::ut::MockService*> getMockService();

		//TODO: Takes ownership -> std::move?
		/*
		 * Takes a set of mock values and creates a corresponding MockService.
		 * This service is then added to _mockServices and returned.
		 */
		template<class T>
		testing::ut::MockService* addMockService(std::vector<T *> mockValues);

	private:
		//list of UT pointers
		std::vector<testing::ut::UnitTestInterface *> _uts;

		//map of "Type"->MockService<Type>
		std::map<std::string, testing::ut::MockService *> _mockServices;
		
};

//INCLUDE ALL INDIVIDUAL UNIT TESTS HERE
#include"uts/tarchVectorUT.cpph"

//Include implementation
#include"UnitTestingService.cpp"
