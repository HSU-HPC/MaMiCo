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
# Parallel build, default is OFF
option(BUILD_WITH_MPI "Build with MPI." OFF)
message(STATUS ">> Build with MPI: " ${BUILD_WITH_MPI})
# Dimension for MD build, default is MDDim3
set(MD_DIM MDDim3 CACHE STRING "Set to MDDim3 (default) or MDDim2.")
message(STATUS ">> SimpleMD dimension variable : " ${MD_DIM})
mark_as_advanced(MD_DIM)
set_property(CACHE MD_DIM PROPERTY STRINGS MDDim3 MDDim2)

#
# Find dependencies
#
find_package(Eigen3 3.3 REQUIRED)
IF(BUILD_WITH_MPI)
	find_package(MPI REQUIRED)
	message(STATUS ">> MPI_INCLUDE_PATH : ${MPI_INCLUDE_PATH}")
	message(STATUS ">> MPI_LIBRARIES : ${MPI_LIBRARIES}")
	message(STATUS ">> MPI_CXX_INCLUDE_PATH : ${MPI_CXX_INCLUDE_PATH}")
	message(STATUS ">> MPI_CXX_INCLUDE_DIRS : ${MPI_CXX_INCLUDE_DIRS}")
	message(STATUS ">> MPI_CXX_LIBRARIES : ${MPI_CXX_LIBRARIES}")
	message(STATUS ">> MPI_CXX_FOUND : ${MPI_CXX_FOUND}")
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
set_property(TARGET simplemd PROPERTY CXX_STANDARD 17)
target_compile_definitions(simplemd PUBLIC ${MD_DIM})
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
	target_compile_definitions(simplemd PUBLIC MDError)
endif()
IF(BUILD_WITH_MPI)
    target_compile_definitions(simplemd PUBLIC MDParallel TarchParallel)
    target_link_libraries(simplemd PUBLIC MPI::MPI_CXX)
ENDIF()
# Check the CXX compiler ID. IDs are defined here: https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_COMPILER_ID.html
message(STATUS ">> CMAKE_CXX_COMPILER : " ${CMAKE_CXX_COMPILER})
message(STATUS ">> CMAKE_CXX_COMPILER_ID : " ${CMAKE_CXX_COMPILER_ID})
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
	target_compile_options(simplemd PUBLIC -Wall -Werror -pedantic -pedantic-errors -Wstrict-aliasing -fstrict-aliasing -Wno-long-long -Wno-unknown-pragmas -Wconversion -Wno-non-virtual-dtor -O3)
	if (CMAKE_BUILD_TYPE STREQUAL "Debug")
		target_compile_options(simplemd PUBLIC -g3)
	endif()
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
	target_compile_options(simplemd PUBLIC -fstrict-aliasing -qpack_semantic=gnu -ipo)
	target_link_options(simplemd PUBLIC -ipo)
	if (CMAKE_BUILD_TYPE STREQUAL "Debug")
		target_compile_options(simplemd PUBLIC -O3)
	else()
		target_compile_options(simplemd PUBLIC -fast -w -Werror-all -align -ansi-alias)
	endif()
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Fujitsu")
	target_compile_options(simplemd PUBLIC -Nclang -Wall -Werror -pedantic -pedantic-errors -Wstrict-aliasing -fstrict-aliasing -Wno-long-long -Wno-unknown-pragmas -Wconversion -Wno-non-virtual-dtor -Wno-sign-conversion -Wno-unused-private-field -mcpu=a64fx)
	target_link_options(simplemd PUBLIC -ipo)
	if (CMAKE_BUILD_TYPE STREQUAL "Debug")
		target_compile_options(simplemd PUBLIC -g3 -O3)
	else()
		target_compile_options(simplemd PUBLIC -Ofast)
	endif()
else()
	message(STATUS ">> Unsupported compiler. Check the CXX variable value.")
endif()



