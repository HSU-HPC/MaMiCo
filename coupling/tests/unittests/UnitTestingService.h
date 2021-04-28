// This file is part of the MaMiCo project. For conditions of distribution
// and use, please see the copyright notice in MaMiCo's main folder

#pragma once

//external headers
#include<vector>
#include<map>
#include<any>
#include<optional>
#include<iostream>
#include<mpi.h>

//General Mamico includes
#include "tarch/configuration/ParseConfiguration.h"
#include "simplemd/configurations/MolecularDynamicsConfiguration.h"
#include "coupling/interface/MDSimulationFactory.h"
#include "coupling/interface/impl/SimpleMD/SimpleMDSolverInterface.h"
#include "coupling/configurations/MaMiCoConfiguration.h"

//Unit Testing includes
#include"MockService.cpph"
#include"UnitTestInterface.h"


//#define DEBUG_UTS


/*
 *
 * Service class managing Unit Test objects and their corresponding MockServices.
 * This is the top-level user interface: To run tests, call runUnitTest(...).
 *
 * All Unit Test instances are created during instanciation of this (cf. .cpp file).
 * Mocks of some types that don't have a MockService are initialized there as well.
 * Currently, those are:
 * 	- primitive C++ types
 * 	- STL types
 * 	- SimpleMD
 *
 * Features currently missing:
 * - simpleMD dummy instances with dim != 3
 * - CS dummy instances
 * - some primitive + a lot of STL types
 *
 * @Author Leonard Hannen, Felix Maurer
 */

namespace testing { namespace ut {
		class UnitTestingService;
}}


//Make this static?
class testing::ut::UnitTestingService {
	public:
		UnitTestingService(MPI_Comm comm = MPI_COMM_WORLD);

		~UnitTestingService() {
			for (auto ut : _uts) delete ut;
			for (auto ms_tuple : _mockServices) delete ms_tuple.second;
			for (auto ic : _indexConversions) delete ic;
			//TODO: do i have to delete all MDs and their interfaces?

			#ifdef DEBUG_UTS
				std::cout << "\x1B[1mUnitTestingService:\x1B[0m Deconstructed." << std::endl;
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

		//Special case: SimpleMD "mocks"
		std::vector<coupling::interface::MDSimulation*> getSimpleMDs() { return _simpleMDs; }
		std::vector<coupling::interface::MDSolverInterface<MY_LINKEDCELL,3> *> getMDInterfaces() { return _mdSolverInterfaces; }
		std::vector<coupling::IndexConversion<3> *> getIndexConversions() { return _indexConversions; }
		std::vector<coupling::configurations::MaMiCoConfiguration<3> *> getMamicoConfigs() { return _mamicoConfigs; }

		/*
		 * Takes a set of mock values and creates a corresponding MockService.
		 * This service is then added to _mockServices and returned.
		 *
		 * The MockService gets freed during this's deconstruction.
		 */
		template<class T>
		testing::ut::MockService* addMockService(std::vector<T *> mockValues);

		//Referring to MPI ranks
		int getRank() { return _rank; }
		int getCommSize() { return _comm_size; }

	private:

		//used for MPI communication
		int _rank;
		MPI_Comm _comm;
		int _comm_size;

		//list of UT pointers
		std::vector<testing::ut::UnitTestInterface *> _uts;

		//map of "Type"->MockService<Type>
		std::map<std::string, testing::ut::MockService *> _mockServices;

		//dummy instances of simpleMD. the following two map 1:1:1, i.e the nth index conversion belongs to the nth MD
  		std::vector<coupling::interface::MDSimulation*> _simpleMDs;
  		std::vector<coupling::interface::MDSolverInterface<MY_LINKEDCELL,3> *> _mdSolverInterfaces;
  		std::vector<coupling::IndexConversion<3> *> _indexConversions;
  		std::vector<coupling::configurations::MaMiCoConfiguration<3> *> _mamicoConfigs;

		//dummy instances of CS	
		//TODO
};

//INCLUDE ALL INDIVIDUAL UNIT TESTS HERE
#include"uts/tarch/la/VectorUT.cpph"
#include"uts/coupling/datastructures/MacroscopicCellUT.cpph"
#include"uts/coupling/datastructures/MacroscopicCellWithLinkedCellsUT.cpph"
#include"uts/coupling/paralleltopology/XYZTopologyUT.cpph"
#include"uts/coupling/cell-mappings/ComputeMassMappingUT.cpph"
#include"uts/coupling/UsherParticleInsertionUT.cpph"

//Include implementation
#include"UnitTestingService.cpp"
