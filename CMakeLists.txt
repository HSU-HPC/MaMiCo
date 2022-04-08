cmake_minimum_required(VERSION 3.13)
project(mamico VERSION 2.0 LANGUAGES CXX)

# CMake Policies
# CMP0074
if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()
# CMP0075
if(POLICY CMP0075)
  cmake_policy(SET CMP0075 NEW)
endif()

#
# Setup Options
#
# Build type
if(NOT CMAKE_BUILD_TYPE)
  message(WARNING "CMAKE_BUILD_TYPE was not specified. Debug build is assumed.")
  set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release)
endif()
message(STATUS ">> Build type: " ${CMAKE_BUILD_TYPE})
# Automatically add the current source and build directories to the include path
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# Enable the MPI library for parallel build, default is OFF
option(BUILD_WITH_MPI "Enable the MPI library for parallel build, default is OFF." OFF)
message(STATUS ">> Build with MPI: " ${BUILD_WITH_MPI})
# Dimension for MD build, default is MDDim3
set(MD_DIM MDDim3 CACHE STRING "Set to MDDim3 or MDDim2, default is MDDim3.")
message(STATUS ">> SimpleMD dimension variable : " ${MD_DIM})
mark_as_advanced(MD_DIM)
set_property(CACHE MD_DIM PROPERTY STRINGS MDDim3 MDDim2)
# Enable the Eigen library used (so far) only for POD filtering, default is off
option(BUILD_WITH_EIGEN "Enable the Eigen library used (so far) only for POD filtering, default is off." OFF)
message(STATUS ">> Build with Eigen : " ${BUILD_WITH_EIGEN})
# Enable MaMiCo python bindings
option(BUILD_WITH_PYBIND11 "Enable the MaMiCo python binding, default is off." OFF)
message(STATUS ">> Build with pybind11 : " ${BUILD_WITH_PYBIND11})
if (BUILD_WITH_PYBIND11 AND NOT BUILD_WITH_MPI)
	message(STATUS ">> Build with pybind11 was enabled but build with MPI was not. Enabling build with MPI.")
	set(BUILD_WITH_MPI "ON")
	message(STATUS ">> Build with MPI: " ${BUILD_WITH_MPI})
endif()

#
# Find dependencies
#
if(BUILD_WITH_EIGEN)
	find_package(Eigen3 REQUIRED)
endif()

if(BUILD_WITH_MPI)
	find_package(MPI REQUIRED)
endif()

if(BUILD_WITH_PYBIND11)
	find_package(Python COMPONENTS Interpreter Development REQUIRED)
	find_package(pybind11 REQUIRED)
endif()

#
# simpleMD target
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
)
add_library(simplemd ${SIMPLE_MD_SOURCES})
set_target_properties(simplemd PROPERTIES CXX_STANDARD 17)
target_compile_definitions(simplemd PRIVATE ${MD_DIM})
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
	target_compile_definitions(simplemd PRIVATE MDError)
endif()
if(BUILD_WITH_MPI)
    target_compile_definitions(simplemd PRIVATE MDParallel TarchParallel)
    target_link_libraries(simplemd PRIVATE MPI::MPI_CXX)
endif()
# Check the CXX compiler ID. IDs are defined here: https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_COMPILER_ID.html
message(STATUS ">> CMAKE_CXX_COMPILER : " ${CMAKE_CXX_COMPILER})
message(STATUS ">> CMAKE_CXX_COMPILER_ID : " ${CMAKE_CXX_COMPILER_ID})
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
	target_compile_options(simplemd PRIVATE -Wall -Werror -pedantic -pedantic-errors -Wstrict-aliasing -fstrict-aliasing -Wno-long-long -Wno-unknown-pragmas -Wconversion -Wno-non-virtual-dtor -O3)
	if (CMAKE_BUILD_TYPE STREQUAL "Debug")
		target_compile_options(simplemd PRIVATE -g3)
	endif()
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
	target_compile_options(simplemd PRIVATE -fstrict-aliasing -qpack_semantic=gnu -ipo)
	target_link_options(simplemd PRIVATE -ipo)
	if (CMAKE_BUILD_TYPE PRIVATE "Debug")
		target_compile_options(simplemd PRIVATE -O3)
	else()
		target_compile_options(simplemd PRIVATE -fast -w -Werror-all -align -ansi-alias)
	endif()
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Fujitsu")
	target_compile_options(simplemd PRIVATE -Nclang -Wall -Werror -pedantic -pedantic-errors -Wstrict-aliasing -fstrict-aliasing -Wno-long-long -Wno-unknown-pragmas -Wconversion -Wno-non-virtual-dtor -Wno-sign-conversion -Wno-unused-private-field -mcpu=a64fx)
	target_link_options(simplemd PRIVATE -ipo)
	if (CMAKE_BUILD_TYPE STREQUAL "Debug")
		target_compile_options(simplemd PRIVATE -g3 -O3)
	else()
		target_compile_options(simplemd PRIVATE -Ofast)
	endif()
else()
	message(STATUS ">> Unsupported compiler.")
endif()

#
# MaMiCo target
#
add_library(mamico INTERFACE)
set_target_properties(mamico PROPERTIES CXX_STANDARD 17)
target_link_libraries(mamico INTERFACE simplemd)
if(BUILD_WITH_EIGEN)
	target_compile_definitions(mamico INTERFACE BUILD_WITH_EIGEN)
	target_link_libraries(mamico INTERFACE Eigen3::Eigen)
endif()
target_compile_definitions(mamico INTERFACE SIMPLE_MD MDDim3)
target_compile_options(mamico INTERFACE -Wall -Wno-unknown-pragmas -Wfatal-errors -O3)
if(BUILD_WITH_MPI)
    target_compile_definitions(mamico INTERFACE MDCoupledParallel TarchParallel MPICH_IGNORE_CXX_SEEK ENABLE_POST_MULTI_INSTANCE_FILTERING)
    target_compile_options(mamico INTERFACE -Werror -Wno-int-in-bool-context -Wno-maybe-uninitialized)
    target_link_libraries(mamico INTERFACE MPI::MPI_CXX)
endif()

#
# Couette target
#
set(COUETTE_SOURCES
		"coupling/configurations/ParticleInsertionConfiguration.cpp"
        "coupling/solvers/CoupledMolecularDynamicsSimulation.cpp"
        "coupling/tests/main_couette.cpp"
)
add_executable(couette ${COUETTE_SOURCES})
set_target_properties(couette PROPERTIES CXX_STANDARD 17)
target_link_libraries(couette PUBLIC mamico)

#
# pybind target
#
if(BUILD_WITH_PYBIND11)
	set(PYBIND_SOURCES11
			"coupling/python-binding/mamico.cpp"
	)
	pybind11_add_module(pybind11
		${PYBIND_SOURCES11}
	)
	target_link_libraries(pybind11 PUBLIC mamico)
endif()

