cmake_minimum_required(VERSION 3.0)
project(mamico VERSION 2.0 LANGUAGES CXX)

if(POLICY CMP0057)
  cmake_policy(SET CMP0057 NEW)
endif()

set(CMAKE_INCLUDE_CURRENT_DIR ON)

if(NOT CMAKE_BUILD_TYPE)
  message(WARNING "CMAKE_BUILD_TYPE was not specified. DebugOptimized build is assumed.")
  set(CMAKE_BUILD_TYPE DebugOptimized CACHE STRING "Build type." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug DebugOptimized RelWithDebInfo Release)
endif()
message(STATUS ">> Build type: " ${CMAKE_BUILD_TYPE})
option(TWO_WAY "Two way coupling." OFF)
message(STATUS ">> Build with two way coupling: " ${TWO_WAY})

#
# Find dependencies
#
find_package(MPI REQUIRED)
find_package(precice REQUIRED)

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

target_compile_definitions(simplemd_interface INTERFACE MDParallel TarchParallel)
target_link_libraries(simplemd_interface INTERFACE MPI::MPI_CXX)
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
target_include_directories(mamico INTERFACE "tarch/logging/")
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
if(TWO_WAY)
  target_compile_definitions(mamico INTERFACE TWO_WAY)
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
