cmake_minimum_required(VERSION 3.0)
project(mamico VERSION 2.0 LANGUAGES CXX)

# CMake Policies
# TODO: check policies
if(POLICY CMP0057)
  cmake_policy(SET CMP0057 NEW)
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
set(MD_SIM SIMPLE_MD CACHE STRING "Choose SIMPLE_MD")
set_property(CACHE MD_SIM PROPERTY STRINGS SIMPLE_MD)
message(STATUS ">> MD simulation used : " ${MD_SIM})
# Enable the Eigen library used (so far) only for POD filtering, default is off
option(BUILD_WITH_EIGEN "Enable the Eigen library used (so far) only for POD filtering, default is off." OFF)
message(STATUS ">> Build with Eigen : " ${BUILD_WITH_EIGEN})
#Enable Adios2 output for megamol
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

find_package(precice REQUIRED)

#
# simplemd interface
#
add_library(simplemd_interface INTERFACE)
target_compile_definitions(simplemd_interface INTERFACE ${MDDim})
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
target_compile_definitions(mamico INTERFACE ${MD_SIM} MDDim3)
target_compile_options(mamico INTERFACE -Wall -Wno-unknown-warning-option -Wno-unknown-pragmas -Wfatal-errors)
if(BUILD_WITH_MPI)
  target_compile_definitions(mamico INTERFACE MDCoupledParallel TarchParallel MPICH_IGNORE_CXX_SEEK ENABLE_POST_MULTI_INSTANCE_FILTERING)
  target_compile_options(mamico INTERFACE -Werror -Wno-int-in-bool-context)
  if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "FujitsuClang")
    target_compile_options(mamico INTERFACE -Wno-uninitialized)
  else()
    target_compile_options(mamico INTERFACE -Wno-maybe-uninitialized)
  endif()
  target_link_libraries(mamico INTERFACE MPI::MPI_CXX)
endif()
target_compile_definitions(mamico INTERFACE BUILD_WITH_PRECICE)
target_link_libraries(mamico INTERFACE precice::precice)
if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "FujitsuClang" )
	target_compile_options(mamico INTERFACE -Wno-undefined-var-template)
endif()
# we always link against libsimplemd, even if ${MD_SIM} not SIMPLE_MD, so that we can use the simpleMD configuration
target_link_libraries(mamico INTERFACE libsimplemd)

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
# Main target
#
set(MAMICO_MAIN_SOURCES
  "coupling/indexing/IndexingService.cpp"
	"coupling/configurations/ParticleInsertionConfiguration.cpp"
  "coupling/solvers/CoupledMolecularDynamicsSimulation.cpp"
  "coupling/scenarios/main.cpp"
)
add_executable(main ${MAMICO_MAIN_SOURCES})
set_target_properties(main PROPERTIES CXX_STANDARD 17)
target_link_libraries(main PUBLIC mamico)
