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
#include"tarch/tinyxml2/tinyxml2.h"

//General Mamico includes
#include "tarch/configuration/ParseConfiguration.h"
#include "simplemd/configurations/MolecularDynamicsConfiguration.h"
#include "coupling/interface/MDSimulationFactory.h"
#include "coupling/interface/impl/SimpleMD/SimpleMDSolverInterface.h"
#include "coupling/configurations/MaMiCoConfiguration.h"

/*
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
*/

//Unit Testing includes
#include"MockService.cpph"
#include"UnitTestInterface.h"


//#define DEBUG_UTS


/*
 *TODO: explanatory interface comment
 *
 * TODO:
 * - simpleMD dummy instances with dim != 3
 * - CS dummy instances
 * @Author Leonard Hannen, Felix Maurer
 */

namespace testing { namespace ut {
		class UnitTestingService;
}}


//Make this static?
class testing::ut::UnitTestingService {
	public:
		UnitTestingService(
			std::vector<std::pair<std::string, std::string>> simplemd_xmls,
			MPI_Comm comm = MPI_COMM_WORLD
		);

		~UnitTestingService() {
			for (auto ut : _uts) delete ut;

			for (auto ms_tuple : _mockServices)delete ms_tuple.second;

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
		int _comm_size;

		//list of UT pointers
		std::vector<testing::ut::UnitTestInterface *> _uts;

		//map of "Type"->MockService<Type>
		std::map<std::string, testing::ut::MockService *> _mockServices;

		//dummy instances of simpleMD
  		std::vector<coupling::interface::MDSimulation*> _simpleMDs;
  		std::vector<coupling::interface::MDSolverInterface<MY_LINKEDCELL,3>* > _mdSolverInterfaces;

		//dummy instances of CS	
};

//INCLUDE ALL INDIVIDUAL UNIT TESTS HERE
#include"uts/tarch/la/VectorUT.cpph"
#include"uts/coupling/datastructures/MacroscopicCellUT.cpph"
#include"uts/coupling/paralleltopology/XYZTopologyUT.cpph"

//Include implementation
#include"UnitTestingService.cpp"
