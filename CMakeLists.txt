cmake_minimum_required(VERSION 3.0)
project(mamico VERSION 2.0 LANGUAGES CXX)

# Make our cmake modules visible to CMake
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

# CMake Policies
# TODO: check policies
if(POLICY CMP0057)
  cmake_policy(SET CMP0057 NEW)
endif()

# CTest
option(BUILD_TESTING "Build testing" OFF)
message(STATUS ">> Build with testing: " ${BUILD_TESTING})
if (BUILD_TESTING)
  include(CTest)
endif()

#
# Setup Options
#
# Build type
if(NOT CMAKE_BUILD_TYPE)
  message(WARNING "CMAKE_BUILD_TYPE was not specified. DebugOptimized build is assumed.")
  set(CMAKE_BUILD_TYPE DebugOptimized CACHE STRING "Build type." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug DebugOptimized RelWithDebInfo Release)
endif()
message(STATUS ">> Build type: " ${CMAKE_BUILD_TYPE})
# Automatically add the current source and build directories to the include path
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# Enable the MPI library for parallel build, default is OFF
option(BUILD_WITH_MPI "Enable the MPI library for parallel build, default is OFF." OFF)
message(STATUS ">> Build with MPI: " ${BUILD_WITH_MPI})
# Dimension for MD build, default is MDDim3
set(MDDim MDDim3 CACHE STRING "Set to MDDim3 or MDDim2, default is MDDim3.")
message(STATUS ">> MD dimension : " ${MDDim})
set_property(CACHE MDDim PROPERTY STRINGS MDDim3 MDDim2)
# Set MD Simulation to use
set(MD_SIM SIMPLE_MD CACHE STRING "Choose SIMPLE_MD, LAMMPS, or LS1_MARDYN")
set_property(CACHE MD_SIM PROPERTY STRINGS SIMPLE_MD LAMMPS_MD LAMMPS_DPD LS1_MARDYN)
message(STATUS ">> MD simulation used : " ${MD_SIM})
# Enable the Eigen library used (so far) only for POD filtering, default is off
option(BUILD_WITH_EIGEN "Enable the Eigen library used (so far) only for POD filtering, default is off." OFF)
message(STATUS ">> Build with Eigen : " ${BUILD_WITH_EIGEN})
# Enable MaMiCo python bindings
option(BUILD_WITH_PYBIND11 "Enable the MaMiCo python binding, default is off." OFF)
message(STATUS ">> Build with pybind11 : " ${BUILD_WITH_PYBIND11})
if (BUILD_WITH_PYBIND11 AND NOT BUILD_WITH_MPI)
	message(STATUS ">> Build with pybind11 was enabled but build with MPI was not. Enabling build with MPI.")
	set(BUILD_WITH_MPI ON CACHE BOOL "Enable the MPI library for parallel build, default is OFF." FORCE)
	message(STATUS ">> Build with MPI: " ${BUILD_WITH_MPI})
endif()
# Enable coupling the OpenFOAM API, default is off
option(BUILD_WITH_OPENFOAM "Enable coupling the OpenFOAM API, default is off." OFF)
message(STATUS ">> Build with OpenFOAM : " ${BUILD_WITH_OPENFOAM})
# Enable coupling with the preCICE API, default is off
option(BUILD_WITH_PRECICE "Enable coupling with the preCICE API, default is off." OFF)
message(STATUS ">> Build with preCICE : " ${BUILD_WITH_PRECICE})
# Enable coupling with the LAMMPS, default is off
option(BUILD_WITH_LAMMPS "Enable coupling with LAMMPS MD, default is off." OFF)
message(STATUS ">> Build with LAMMPS : " ${BUILD_WITH_LAMMPS})
if(("${MD_SIM}" STREQUAL "LAMMPS_MD" OR "${MD_SIM}" STREQUAL "LAMMPS_DPD") AND NOT BUILD_WITH_LAMMPS)
  message(STATUS ">> MD_SIM is LAMMPS but BUILD_WITH_LAMMPS was disabled. Enabling build with LAMMPS.")
  set(BUILD_WITH_LAMMPS ON CACHE BOOL "Enable coupling with LAMMPS MD, default is off." FORCE)
  message(STATUS ">> Build with LAMMPS : " ${BUILD_WITH_LAMMPS})
endif()
if("${MD_SIM}" STREQUAL "SIMPLE_MD" AND BUILD_WITH_LAMMPS)
  message(WARNING ">> BUILD_WITH_LAMMPS is enabled but MD_SIM is SIMPLE_MD. Are you sure?")
endif()
# Compiler list to be used for compiletest
set(COMPILETEST_COMPILERS "g++ clang++" CACHE STRING "Compiler list to be used for compiletest")
message(STATUS ">> Compiler list to be used for compiletest : " ${COMPILETEST_COMPILERS})
# Enable coupling with the preCICE API, default is off
option(BUILD_WITH_TEST_INDEXING "Enable indexing tests, default is off." OFF)
message(STATUS ">> Build with indexing tests : " ${BUILD_WITH_TEST_INDEXING})
# Enable Adios2 output for megamol
option(BUILD_WITH_ADIOS2 "Enable output in the Adios2 format." OFF)
message(STATUS ">> Build with Adios2: " ${BUILD_WITH_ADIOS2})

#
# Find dependencies
#
if(BUILD_WITH_EIGEN)
	find_package(Eigen3 REQUIRED NO_MODULE)
endif()

if(BUILD_WITH_MPI)
	find_package(MPI REQUIRED)
endif()

if(BUILD_WITH_ADIOS2)
  enable_language(C)
  find_package(ADIOS2 REQUIRED)
endif()

if(BUILD_WITH_PYBIND11)
	find_package(Python COMPONENTS Interpreter Development REQUIRED)
	find_package(pybind11 REQUIRED)
endif()

if(BUILD_WITH_OPENFOAM)
	# Check valid OpenFOAM
	if(DEFINED ENV{WM_PROJECT_DIR})
		message(STATUS ">> OpenFOAM found in : " $ENV{WM_PROJECT_DIR})
	else()
		message(FATAL_ERROR "The OpenFOAM bashrc is not sourced or OpenFOAM is not installed")
	endif()
endif()

if(BUILD_WITH_PRECICE)
	find_package(precice REQUIRED)
endif()

if(BUILD_WITH_LAMMPS)
  find_package(PkgConfig REQUIRED)
  pkg_search_module(LAMMPS REQUIRED liblammps)
endif()

if("${MD_SIM}" STREQUAL "LS1_MARDYN")
  set(LS1_SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/ls1" CACHE PATH "ls1 source directory")
  set(LS1_BUILD_DIR_NAME "build" CACHE STRING "ls1 build directory name")
  load_cache(${LS1_SRC_DIR}/${LS1_BUILD_DIR_NAME} READ_WITH_PREFIX LS1_CACHE_ ENABLE_AUTOPAS ENABLE_MPI PRECISION ENABLE_VTK CMAKE_BUILD_TYPE)
  if(NOT LS1_CACHE_PRECISION)
    message(FATAL_ERROR "Could not read the ls1 build cache. Please ensure the path to the ls1 root folder and the name of the build folder are correct.")
  endif()
  if(LS1_CACHE_ENABLE_MPI AND NOT BUILD_WITH_MPI)
    message(FATAL_ERROR "ls1 built with MPI, but MaMiCo does not have MPI enabled! Aborting...")
  endif()
  if(NOT LS1_CACHE_ENABLE_MPI AND BUILD_WITH_MPI)
    message(FATAL_ERROR "ls1 not built with MPI, but MaMiCo has MPI enabled! Aborting...")
  endif()
  option(LS1_BUILT_WITH_AUTOPAS "Specify whether ls1 was compiled with AutoPas support." ${LS1_CACHE_ENABLE_AUTOPAS})
  find_library(MarDyn REQUIRED NAMES "libMarDyn.a" PATHS ${LS1_SRC_DIR}/${LS1_BUILD_DIR_NAME}/src )
  if(LS1_BUILT_WITH_AUTOPAS)
    find_library(autopas REQUIRED NAMES "libautopas.a" PATHS ${LS1_SRC_DIR}/${LS1_BUILD_DIR_NAME}/_deps/autopasfetch-build/src/autopas)
    if("${LS1_CACHE_CMAKE_BUILD_TYPE}" STREQUAL "Debug")
      find_library(spdlog REQUIRED NAMES "libspdlogd.a" PATHS ${LS1_SRC_DIR}/${LS1_BUILD_DIR_NAME}/_deps/spdlog-build )
    else()
      find_library(spdlog REQUIRED NAMES "libspdlog.a" PATHS ${LS1_SRC_DIR}/${LS1_BUILD_DIR_NAME}/_deps/spdlog-build)
    endif()
    find_library(harmony REQUIRED NAMES "libharmony.a" PATHS ${LS1_SRC_DIR}/${LS1_BUILD_DIR_NAME}/_deps/autopasfetch-build/libs/harmony/include/lib )
  endif()
  if(LS1_CACHE_ENABLE_VTK)
    find_library(VTK_LIB xerces-c HINTS $ENV{XERCES_LIBDIR})
    if(NOT VTK_LIB)
      message(FATAL_ERROR "xerces-c lib not found. Disable VTK support, if you do not need it.")
    endif()
  endif()
  message(STATUS ">> ls1 found and initialized")
endif()

if(BUILD_TESTING)
  find_package(CppUnit REQUIRED)
endif()

#
# simplemd interface
#
add_library(simplemd_interface INTERFACE)
target_compile_definitions(simplemd_interface INTERFACE MDDim3)
target_compile_options(simplemd_interface INTERFACE -fpic)
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
	target_compile_definitions(simplemd_interface INTERFACE MDDebug TarchDebug MDError)
  target_compile_options(simplemd_interface INTERFACE -Og -g3)
elseif(CMAKE_BUILD_TYPE STREQUAL "DebugOptimized")
  target_compile_definitions(simplemd_interface INTERFACE MDError)
  target_compile_options(simplemd_interface INTERFACE -O3 -g3)
elseif(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
  target_compile_options(simplemd_interface INTERFACE -O3 -g3)
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
  target_compile_options(simplemd_interface INTERFACE -Ofast)
else()
  message(FATAL_ERROR "Unsupported build type")
endif()

if(BUILD_WITH_MPI)
    target_compile_definitions(simplemd_interface INTERFACE MDParallel TarchParallel)
    target_link_libraries(simplemd_interface INTERFACE MPI::MPI_CXX)
endif()
if(BUILD_WITH_ADIOS2)
  if(BUILD_WITH_MPI)
    target_link_libraries(simplemd_interface INTERFACE adios2::cxx11_mpi MPI::MPI_CXX)
  else()
    target_link_libraries(simplemd_interface INTERFACE adios2::cxx11)
  endif()
endif()
# Check the CXX compiler ID. IDs are defined here: https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_COMPILER_ID.html
message(STATUS ">> CMAKE_CXX_COMPILER : " ${CMAKE_CXX_COMPILER})
message(STATUS ">> CMAKE_CXX_COMPILER_ID : " ${CMAKE_CXX_COMPILER_ID})
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
	target_compile_options(simplemd_interface INTERFACE -Wall -Werror -pedantic -pedantic-errors -Wstrict-aliasing -fstrict-aliasing -Wno-long-long -Wno-unknown-pragmas -Wconversion -Wno-non-virtual-dtor)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
	target_compile_options(simplemd_interface INTERFACE -fstrict-aliasing -qpack_semantic=gnu -ipo -align -ansi-alias) #-w
	target_link_options(simplemd_interface INTERFACE -ipo)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "FujitsuClang")
	target_compile_options(simplemd_interface INTERFACE -Wall -Werror -pedantic -pedantic-errors -Wstrict-aliasing -fstrict-aliasing -Wno-long-long -Wno-unknown-pragmas -Wconversion -Wno-non-virtual-dtor -Wno-sign-conversion -Wno-unused-private-field -mcpu=a64fx)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
	target_compile_options(simplemd_interface INTERFACE -Wall -pedantic -pedantic-errors -Wstrict-aliasing -fstrict-aliasing -Wno-long-long -Wno-unknown-pragmas -Wconversion -Wno-non-virtual-dtor -Wno-sign-conversion -Wno-unused-private-field -Wno-gnu-array-member-paren-init)
else()
	message(STATUS ">> Unsupported compiler.")
endif()

#
# simplemd library target
#
set(SIMPLE_MD_SOURCES
  "simplemd/BoundaryTreatment.cpp"
  "simplemd/MolecularDynamicsSimulation.cpp"
  "simplemd/ProfilePlotter.cpp"
  "simplemd/cell-mappings/DeleteMoleculesMapping.cpp"
	"simplemd/cell-mappings/LennardJonesForceMapping.cpp"
	"simplemd/cell-mappings/LennardJonesPotentialEnergyMapping.cpp"
	"simplemd/cell-mappings/PeriodicAndParallelBoundaryFillCellsMapping.cpp"
	"simplemd/cell-mappings/PeriodicBoundaryEmptyCellsMapping.cpp"
	"simplemd/cell-mappings/RDFMapping.cpp"
	"simplemd/configurations/DomainConfiguration.cpp"
	"simplemd/configurations/ExternalForceConfiguration.cpp"
	"simplemd/configurations/MolecularDynamicsConfiguration.cpp"
	"simplemd/configurations/MoleculeConfiguration.cpp"
	"simplemd/configurations/MPIConfiguration.cpp"
	"simplemd/configurations/ProfilePlotterConfiguration.cpp"
	"simplemd/configurations/RDFConfiguration.cpp"
	"simplemd/configurations/SimulationConfiguration.cpp"
	"simplemd/configurations/VTKConfiguration.cpp"
	"simplemd/molecule-mappings/UpdateLinkedCellListsMapping.cpp"
	"simplemd/molecule-mappings/VelocityStoermerVerletMapping.cpp"
	"simplemd/molecule-mappings/VTKMoleculeWriter.cpp"
	"simplemd/molecule-mappings/WriteCheckPointMapping.cpp"
	"simplemd/services/ExternalForceService.cpp"
	"simplemd/services/LinkedCellService.cpp"
	"simplemd/services/MolecularPropertiesService.cpp"
	"simplemd/services/MoleculeService.cpp"
	"simplemd/services/ParallelAndLocalBufferService.cpp"
	"simplemd/services/ParallelTopologyService.cpp"
  "tarch/tinyxml2/tinyxml2.cpp"
  "tarch/utils/RandomNumberService.cpp"
  "simplemd/configurations/Adios2Configuration.cpp"
)
if(BUILD_WITH_ADIOS2)
  list(APPEND SIMPLE_MD_SOURCES "simplemd/molecule-mappings/Adios2Writer.cpp")
endif()

add_library(libsimplemd ${SIMPLE_MD_SOURCES})
target_link_libraries(libsimplemd PRIVATE simplemd_interface)
set_target_properties(libsimplemd PROPERTIES CXX_STANDARD 17)
set_target_properties(libsimplemd PROPERTIES OUTPUT_NAME simplemd)

#
# simplemd executable target
#
add_executable(simplemd ${SIMPLE_MD_SOURCES} "simplemd/main.cpp")
target_link_libraries(simplemd PRIVATE simplemd_interface)
set_target_properties(simplemd PROPERTIES CXX_STANDARD 17)
set_target_properties(simplemd PROPERTIES OUTPUT_NAME simplemd)

#
# multi-simplemd target
#
add_executable(multi-simplemd ${SIMPLE_MD_SOURCES} "simplemd/main_multi.cpp")
target_link_libraries(multi-simplemd PRIVATE simplemd_interface)
set_target_properties(multi-simplemd PROPERTIES CXX_STANDARD 17)

#
# MaMiCo interface
#
add_library(mamico INTERFACE)
if(BUILD_WITH_EIGEN)
	target_compile_definitions(mamico INTERFACE BUILD_WITH_EIGEN)
	target_link_libraries(mamico INTERFACE Eigen3::Eigen)
endif()
if(BUILD_WITH_ADIOS2)
  add_compile_definitions(BUILD_WITH_ADIOS2)
endif()
target_compile_definitions(mamico INTERFACE SIMPLE_MD MDDim3)
target_compile_options(mamico INTERFACE -Wall -Wno-unknown-warning-option -Wno-unknown-pragmas -Wfatal-errors)
target_compile_definitions(mamico INTERFACE MDCoupledParallel TarchParallel MPICH_IGNORE_CXX_SEEK ENABLE_POST_MULTI_INSTANCE_FILTERING)
target_compile_options(mamico INTERFACE -Werror -Wno-int-in-bool-context)
if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "FujitsuClang")
  target_compile_options(mamico INTERFACE -Wno-uninitialized)
else()
  target_compile_options(mamico INTERFACE -Wno-maybe-uninitialized)
endif()
target_link_libraries(mamico INTERFACE MPI::MPI_CXX)
target_compile_definitions(mamico INTERFACE BUILD_WITH_PRECICE)
target_link_libraries(mamico INTERFACE precice::precice)
if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "FujitsuClang" )
	target_compile_options(mamico INTERFACE -Wno-undefined-var-template)
endif()
target_link_libraries(mamico INTERFACE libsimplemd)
if("${MD_SIM}" STREQUAL "LAMMPS_MD")
  target_link_libraries(mamico INTERFACE ${LAMMPS_LINK_LIBRARIES})
  target_include_directories(mamico INTERFACE ${LAMMPS_INCLUDE_DIRS})
  #target_compile_options(mamico INTERFACE ${LAMMPS_CFLAGS_OTHER})
elseif("${MD_SIM}" STREQUAL "LS1_MARDYN")
  target_compile_definitions(mamico INTERFACE MAMICO_COUPLING)
  if ("${LS1_CACHE_PRECISION}" STREQUAL "DOUBLE")
    target_compile_definitions(mamico INTERFACE MARDYN_DPDP)
  elseif("${LS1_CACHE_PRECISION}" STREQUAL "SINGLE")
    target_compile_definitions(mamico INTERFACE MARDYN_SPSP)
  elseif("${LS1_CACHE_PRECISION}" STREQUAL "MIXED")
    target_compile_definitions(mamico INTERFACE MARDYN_SPDP)
  endif()
  target_include_directories(mamico INTERFACE ${LS1_SRC_DIR}/src ${LS1_SRC_DIR}/libs/rapidxml)
  target_link_directories(mamico INTERFACE ${LS1_SRC_DIR}/${LS1_BUILD_DIR_NAME}/src)
  target_link_libraries(mamico INTERFACE MarDyn)
  if(BUILD_WITH_MPI)
    target_compile_definitions(mamico INTERFACE ENABLE_MPI)
  endif()
  if(LS1_BUILT_WITH_AUTOPAS)
    target_compile_definitions(mamico INTERFACE MARDYN_AUTOPAS)
    target_include_directories(mamico INTERFACE ${LS1_SRC_DIR}/${LS1_BUILD_DIR_NAME}/_deps/autopasfetch-src/src ${LS1_SRC_DIR}/${LS1_BUILD_DIR_NAME}/_deps/spdlog-src/include)
    target_link_directories(mamico INTERFACE ${LS1_SRC_DIR}/${LS1_BUILD_DIR_NAME}/_deps/autopasfetch-build/src/autopas ${LS1_SRC_DIR}/${LS1_BUILD_DIR_NAME}/_deps/autopasfetch-build/libs/harmony/include/lib ${LS1_SRC_DIR}/${LS1_BUILD_DIR_NAME}/_deps/spdlog-build)
    target_link_libraries(mamico INTERFACE autopas pthread spdlog harmony)
    #pthread is required by spdlog
  endif()
  if(LS1_CACHE_ENABLE_VTK)
    target_compile_definitions(mamico INTERFACE VTK)
    target_link_libraries(mamico INTERFACE ${VTK_LIB})
  endif()
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  target_compile_definitions(mamico INTERFACE MDCoupledDebug TarchDebug NLM_DEBUG DEBUG_FILTER_PIPELINE)
  target_compile_options(mamico INTERFACE -Og -g3)
elseif(CMAKE_BUILD_TYPE STREQUAL "DebugOptimized")
  target_compile_definitions(mamico INTERFACE MDCoupledError)
  target_compile_options(mamico INTERFACE -O3 -g3)
elseif(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
  target_compile_options(mamico INTERFACE -O3 -g3)
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
  target_compile_options(mamico INTERFACE -Ofast)
else()
  message(FATAL_ERROR "Unsupported build type")
endif()

#
# Couette target
#
set(COUETTE_SOURCES
	"coupling/indexing/IndexingService.cpp"
	"coupling/configurations/ParticleInsertionConfiguration.cpp"
  "coupling/solvers/CoupledMolecularDynamicsSimulation.cpp"
  "coupling/scenario/main_couette.cpp"
)
add_executable(couette ${COUETTE_SOURCES})
set_target_properties(couette PROPERTIES CXX_STANDARD 17)
if (CMAKE_CXX_COMPILER_ID STREQUAL "FujitsuClang")
	target_compile_options(couette PRIVATE -std=c++17)
endif()
target_link_libraries(couette PRIVATE mamico)

#
# test_lammps target
#
if("${MD_SIM}" STREQUAL "LAMMPS_MD")
  set(TEST_LAMMPS_SOURCES
    "tarch/tinyxml2/tinyxml2.cpp"
    "coupling/configurations/ParticleInsertionConfiguration.cpp"
    "tarch/utils/RandomNumberService.cpp"
    "coupling/test/integration/main_lammps.cpp"
    "coupling/indexing/IndexingService.cpp"
  )
  add_executable(test_lammps ${TEST_LAMMPS_SOURCES})
  set_target_properties(test_lammps PROPERTIES CXX_STANDARD 17)
  if (CMAKE_CXX_COMPILER_ID STREQUAL "FujitsuClang")
    target_compile_options(test_lammps PRIVATE -std=c++17)
  endif()
  target_link_libraries(test_lammps PRIVATE mamico)
endif()

#
# Main target
#
set(MAMICO_MAIN_SOURCES
  "coupling/indexing/IndexingService.cpp"
	"coupling/configurations/ParticleInsertionConfiguration.cpp"
  "coupling/solvers/CoupledMolecularDynamicsSimulation.cpp"
  "coupling/solvers/LBCouetteSolver.h"
  "coupling/solvers/IcoFoamBufferSetup.h"
  "coupling/solvers/CoupledMolecularDynamicsSimulation.h"
  "coupling/solvers/CouetteSolver.h"
  "coupling/solvers/CouetteSolverInterface.h"
  "coupling/solvers/LBCouetteSolverInterface.h"
  "coupling/solvers/NumericalSolver.h"
  "test/integration/TestEspressoMDMolecule.h"
  "test/integration/Test.h"
  "test/integration/MacroscopicCellsTest.h"
  "test/integration/TestEspressoAddDeleteMolecule.h"
  "test/integration/main_lbcouettesolver.cpp"
  "test/integration/TestLammpsCalculateForceEnergy.h"
  "test/integration/CompileTest.h"
  "test/integration/TestEspressoMDSolverInterface.h"
  "test/integration/main_lammps.cpp"
  "test/integration/TestLammpsMoleculeIterator.h"
  "test/integration/TestMarDynMolecule.h"
  "test/integration/TestLammpsAddDeleteAtom.h"
  "test/integration/main_mardyn.cpp"
  "test/integration/TestEspressoForceEnergyCalculation.h"
  "test/integration/TestLammps.h"
  "test/integration/TestLammpsInitialVelocity.h"
  "test/integration/TestLBCouetteSolver.h"
  "test/integration/cell_idx_iter_benchmark.cpp"
  "test/integration/main_espresso.cpp"
  "test/integration/TestEspresso.h"
  "test/integration/TestMarDyn.h"
  "test/integration/TestLammpsGhost.h"
  "test/integration/TestMarDynAddDeleteMolecule.h"
  "test/integration/TestEspressoMoleculeIterator.h"
  "test/integration/TestMarDynMoleculeIterator.h"
  "coupling/scenario/CouetteScenario.h"
  "coupling/scenario/main_couette.cpp"
  "coupling/scenario/Scenario.h"
  "coupling/MPIIndexConversion.h"
  "coupling/InstanceHandling.h"
  "coupling/MomentumInsertion.h"
  "simplemd/services/ExternalForceService.cpp"
  "simplemd/services/LinkedCellService.h"
  "simplemd/services/ExternalForceService.h"
  "simplemd/services/MoleculeService.h"
  "simplemd/services/ParallelTopologyService.cpp"
  "simplemd/services/MolecularPropertiesService.h"
  "simplemd/services/LinkedCellService.cpp"
  "simplemd/services/MolecularPropertiesService.cpp"
  "simplemd/services/ParallelTopologyService.h"
  "simplemd/services/ParallelAndLocalBufferService.cpp"
  "simplemd/services/ParallelAndLocalBufferService.h"
  "simplemd/services/MoleculeService.cpp"
  "simplemd/MolecularDynamicsDefinitions.h"
  "simplemd/ProfilePlotter.cpp"
  "simplemd/MolecularDynamicsSimulation.cpp"
  "simplemd/BoundaryTreatment.cpp"
  "simplemd/cell-mappings/ResetPotentialEnergyMapping.h"
  "simplemd/cell-mappings/ParallelBoundaryEmptyCellsMapping.h"
  "simplemd/cell-mappings/RDFMapping.h"
  "simplemd/cell-mappings/CollectMoleculesMapping.h"
  "simplemd/cell-mappings/LennardJonesForceMapping.cpp"
  "simplemd/cell-mappings/ProfilePlotterMapping.h"
  "simplemd/cell-mappings/LennardJonesPotentialEnergyMapping.cpp"
  "simplemd/cell-mappings/ComputeTemperatureMapping.h"
  "simplemd/cell-mappings/ComputeMeanVelocityMapping.h"
  "simplemd/cell-mappings/CopyMoleculesMapping.h"
  "simplemd/cell-mappings/VaryCheckpointMapping.h"
  "simplemd/cell-mappings/EmptyLinkedListsMapping.h"
  "simplemd/cell-mappings/DeleteMoleculesMapping.cpp"
  "simplemd/cell-mappings/PeriodicBoundaryEmptyCellsMapping.h"
  "simplemd/cell-mappings/LennardJonesForceMapping.h"
  "simplemd/cell-mappings/RDFMapping.cpp"
  "simplemd/cell-mappings/DeleteMoleculesMapping.h"
  "simplemd/cell-mappings/LennardJonesPotentialEnergyMapping.h"
  "simplemd/cell-mappings/PeriodicAndParallelBoundaryFillCellsMapping.cpp"
  "simplemd/cell-mappings/PeriodicBoundaryEmptyCellsMapping.cpp"
  "simplemd/cell-mappings/PeriodicAndParallelBoundaryFillCellsMapping.h"
  "simplemd/MolecularDynamicsSimulation.h"
  "simplemd/MolecularProperties.h"
  "simplemd/molecule-mappings/InitialPositionAndForceUpdate.h"
  "simplemd/molecule-mappings/SetMeanVelocityMapping.h"
  "simplemd/molecule-mappings/VTKMoleculeWriter.h"
  "simplemd/molecule-mappings/VelocityStoermerVerletMapping.h"
  "simplemd/molecule-mappings/VelocityStoermerVerletMapping.cpp"
  "simplemd/molecule-mappings/UpdateLinkedCellListsMapping.cpp"
  "simplemd/molecule-mappings/ComputeMeanVelocityMapping.h"
  "simplemd/molecule-mappings/UpdateLinkedCellListsMapping.h"
  "simplemd/molecule-mappings/WriteCheckPointMapping.cpp"
  "simplemd/molecule-mappings/WriteCheckPointMapping.h"
  "simplemd/molecule-mappings/VTKMoleculeWriter.cpp"
  "simplemd/molecule-mappings/Adios2Writer.cpp"
  "simplemd/molecule-mappings/Adios2Writer.h"
  "simplemd/BoundaryTreatment.h"
  "simplemd/LinkedCell.h"
  "simplemd/main_multi.cpp"
  "simplemd/configurations/CheckpointConfiguration.h"
  "simplemd/configurations/MoleculeConfiguration.h"
  "simplemd/configurations/MolecularDynamicsConfiguration.cpp"
  "simplemd/configurations/MolecularDynamicsConfiguration.h"
  "simplemd/configurations/DomainConfiguration.h"
  "simplemd/configurations/MoleculeConfiguration.cpp"
  "simplemd/configurations/SimulationConfiguration.h"
  "simplemd/configurations/DomainConfiguration.cpp"
  "simplemd/configurations/VTKConfiguration.h"
  "simplemd/configurations/ProfilePlotterConfiguration.cpp"
  "simplemd/configurations/RDFConfiguration.cpp"
  "simplemd/configurations/ProfilePlotterConfiguration.h"
  "simplemd/configurations/MPIConfiguration.cpp"
  "simplemd/configurations/ExternalForceConfiguration.h"
  "simplemd/configurations/SimulationConfiguration.cpp"
  "simplemd/configurations/MPIConfiguration.h"
  "simplemd/configurations/ExternalForceConfiguration.cpp"
  "simplemd/configurations/RDFConfiguration.h"
  "simplemd/configurations/VTKConfiguration.cpp"
  "simplemd/configurations/Adios2Configuration.cpp"
  "simplemd/configurations/Adios2Configuration.h"
  "simplemd/main.cpp"
  "simplemd/Molecule.h"
  "simplemd/ProfilePlotter.h"
  "simplemd/MolecularDynamicsUserInput.h"
)
foreach(clangformat_source ${CLANGFORMAT_SOURCES})
  get_filename_component(clangformat_source ${clangformat_source} ABSOLUTE)
  list(APPEND CLANGFORMAT_SOURCES_ABS ${clangformat_source})
endforeach()
add_custom_target(
  clangformat
  COMMAND /usr/bin/env clang-format --verbose
  -i
  ${CLANGFORMAT_SOURCES_ABS}
)

#
# compiletest target
#
get_filename_component(COMPILETEST_SCRIPT "cmake/CompileTest.sh" ABSOLUTE)
add_custom_target(
  compiletest
  COMMAND ${COMPILETEST_SCRIPT}
  ${COMPILETEST_COMPILERS}
)

#
# Benchmarks target
#
set(BENCHMARKS_SOURCES
	"coupling/indexing/IndexingService.cpp"
	"coupling/configurations/ParticleInsertionConfiguration.cpp"
  "test/integration/cell_idx_iter_benchmark.cpp"
)
add_executable(benchmarks ${BENCHMARKS_SOURCES})
set_target_properties(benchmarks PROPERTIES CXX_STANDARD 17)
target_link_libraries(benchmarks PRIVATE mamico)

#
# Testing
#
if (BUILD_TESTING)    
  file(GLOB_RECURSE TEST_SOURCES "test/*Test.cpp")
  add_executable(testmamico "test/unit/main.cpp" ${TEST_SOURCES})
  target_link_libraries(testmamico PRIVATE mamico)
  target_include_directories(testmamico PRIVATE ${CPPUNIT_INCLUDE_DIRS})
  target_link_directories(testmamico PRIVATE ${CPPUNIT_LIBRARY})
  target_link_libraries(testmamico PRIVATE cppunit)
  foreach(TEST ${TEST_SOURCES})
    get_filename_component(TEST_NAME ${TEST} NAME_WE)
    if (BUILD_WITH_MPI)
      add_test(NAME ${TEST_NAME} COMMAND mpirun -np 4 ./testmamico ${TEST_NAME})
    else()
      add_test(NAME ${TEST_NAME} COMMAND ./testmamico ${TEST_NAME})
    endif()
  endforeach()
endif()