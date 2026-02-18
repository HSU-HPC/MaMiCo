cmake_minimum_required(VERSION 3.0...3.28)
project(mamico VERSION 2.0 LANGUAGES CXX)

# Make our cmake modules visible to CMake
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

# CMake Policies
# TODO: check policies
if(POLICY CMP0057)
  cmake_policy(SET CMP0057 NEW)
endif()

# Testing and coverage
option(BUILD_TESTING "Build testing" OFF)
message(STATUS ">> Build with testing: " ${BUILD_TESTING})
if (BUILD_TESTING)
  include(CTest)
  option(ENABLE_COVERAGE "Enable code coverage, requires GCC and CMAKE_BUILD_TYPE=Debug." OFF)
  message(STATUS ">> Enable code coverage: " ${ENABLE_COVERAGE})
endif()
if (ENABLE_COVERAGE)
  if (NOT BUILD_TESTING OR NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(WARNING ">> Disabling code coverage.\n CMAKE_BUILD_TYPE=Debug and BUILD_TESTING=ON required.")
    set(ENABLE_COVERAGE OFF CACHE BOOL "Enable code coverage, requires GCC and CMAKE_BUILD_TYPE=Debug." FORCE)
    message(STATUS ">> Enable code coverage: " ${ENABLE_COVERAGE})
  endif()
  if (NOT CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    message(FATAL_ERROR ">> ENABLE_COVERAGE only works with GCC. Change your compiler or set ENABLE_COVERAGE=OFF")
  endif()
endif()

#
# Setup Options
#
# Build type
if(NOT CMAKE_BUILD_TYPE)
  message(WARNING "CMAKE_BUILD_TYPE was not specified. DebugOptimized build is assumed.")
  set(CMAKE_BUILD_TYPE DebugOptimized CACHE STRING "Build type." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS DebugMD Debug DebugOptimized RelWithDebInfo Release)
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
# Enable the Eigen library used for POD filtering, and for coupling with AutoPas, default is off
option(BUILD_WITH_EIGEN "Enable the Eigen library used for POD filtering, and for coupling with AutoPas, default is off." OFF)
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

# Build mode list to be used for compiletest
set(COMPILETEST_MODES "Debug DebugOptimized RelWithDebInfo Release" CACHE STRING "Build mode list to be used for compiletest")
message(STATUS ">> Build mode list to be used for compiletest : " ${COMPILETEST_MODES})

# Enable coupling with pint debug, default is off
option(BUILD_WITH_PINT_DEBUG "Enable parallel-in-time integration debugging, default is off." OFF)
message(STATUS ">> Build with PINT_DEBUG : " ${BUILD_WITH_PINT_DEBUG})

# Enable Adios2 output for megamol
option(BUILD_WITH_ADIOS2 "Enable output in the Adios2 format." OFF)
message(STATUS ">> Build with Adios2: " ${BUILD_WITH_ADIOS2})

#
# Find dependencies
#

if(BUILD_WITH_MPI)
  find_package(MPI REQUIRED)
endif()

set(KOKKOS_TARGET OPENMP CACHE STRING "Choose CUDA, OPENMP, or SERIAL")
set_property(CACHE KOKKOS_TARGET PROPERTY STRINGS CUDA OPENMP SERIAL)
message(STATUS ">> KOKKOS_TARGET: " ${KOKKOS_TARGET})

find_package(Kokkos CONFIG) # Try to find Kokkos externally
if(Kokkos_FOUND)
    message(STATUS "Found Kokkos: ${Kokkos_DIR} (version \"${Kokkos_VERSION}\")")
else()
    message(STATUS "Kokkos not found externally. Fetching via FetchContent.")
    include(FetchContent)
    set(Kokkos_ARCH_NATIVE ON CACHE STRING "Always build Kokkos for the CPU on which MaMiCo is being compiled" FORCE)

    # Backends
    set(Kokkos_ENABLE_SERIAL ON CACHE STRING "Always build Kokkos with Serial support" FORCE)
    if(KOKKOS_TARGET STREQUAL "SERIAL")
      set(Kokkos_ENABLE_OPENMP OFF CACHE STRING "Build Kokkos without OpenMP support when targeting Serial" FORCE)
    else()
      set(Kokkos_ENABLE_OPENMP ON CACHE STRING "Build Kokkos with OpenMP support when not targeting Serial" FORCE)
    endif()
    if(KOKKOS_TARGET STREQUAL "CUDA")
      set(Kokkos_ENABLE_CUDA ON CACHE STRING "Build MaMiCo/Kokkos with CUDA support" FORCE)
      set(Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE ON CACHE STRING "MaMiCo requires CUDA relocatable device code (RDC)" FORCE)
    else()
      set(Kokkos_ENABLE_CUDA OFF CACHE STRING "Build MaMiCo/Kokkos without CUDA support" FORCE)
    endif()

    # Debugging
    if(CMAKE_BUILD_TYPE MATCHES "^Debug" AND (NOT CMAKE_BUILD_TYPE STREQUAL "DebugOptimized"))
      set(KOKKOS_DEBUG "ON")
    else()
      set(KOKKOS_DEBUG "OFF")
    endif()
    set(Kokkos_ENABLE_DEBUG ${KOKKOS_DEBUG} CACHE STRING "Debug builds only" FORCE)
    set(Kokkos_ENABLE_DEBUG_BOUNDS_CHECK ${KOKKOS_DEBUG} CACHE STRING "Debug builds only" FORCE)
    if(CMAKE_BUILD_TYPE MATCHES "^Debug" AND (${KOKKOS_TARGET} STREQUAL "CUDA"))
      set(Kokkos_ENABLE_CUDA_CONSTEXPR ON CACHE STRING "Allow host constexpr code for CUDA" FORCE)
    endif()

    FetchContent_Declare(
        Kokkos
        URL https://github.com/kokkos/kokkos/archive/refs/tags/5.0.2.tar.gz
    )
    FetchContent_MakeAvailable(Kokkos)

    # Hide Kokkos related options in ccmake
    get_cmake_property(_vars CACHE_VARIABLES)
    foreach(_var ${_vars})
      if(_var MATCHES "^FETCHCONTENT_" OR _var MATCHES "^Kokkos_")
        mark_as_advanced(${_var})
      endif()
    endforeach()
    target_compile_options(kokkoscore PUBLIC -fPIC)
    target_link_options(kokkoscore PUBLIC -fPIC)
endif()

if(NOT "${MD_SIM}" STREQUAL "LS1_MARDYN" AND BUILD_WITH_ADIOS2) #ls1 can come bundled with own adios2, hence finding adios2 later for ls1
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
  set(LS1_BUILD_DIR_NAMES "build-seq;build-par;build" CACHE STRING "ls1 build directory names. First one working is used")
  list(LENGTH LS1_BUILD_DIR_NAMES NUM_BUILD_DIRS)
  # we have to use the last one anyway (no other options left in list), so skip that check
  MATH(EXPR NUM_BUILD_DIRS "${NUM_BUILD_DIRS}-1")
  set(BUILD_CANDIDATE_IDX 0)
  list(GET LS1_BUILD_DIR_NAMES ${BUILD_CANDIDATE_IDX} LS1_BUILD_DIR_NAME)
  while(BUILD_CANDIDATE_IDX LESS NUM_BUILD_DIRS)
    # check if this build dir works
    if(EXISTS "${LS1_SRC_DIR}/${LS1_BUILD_DIR_NAME}")
      load_cache(${LS1_SRC_DIR}/${LS1_BUILD_DIR_NAME} READ_WITH_PREFIX CANDIDATE_CACHE_ ENABLE_MPI PRECISION)
    endif()
    if(CANDIDATE_CACHE_PRECISION AND (CANDIDATE_CACHE_ENABLE_MPI STREQUAL BUILD_WITH_MPI))
      # this candidate build dir seems good, stop searching and try it
      break()
    endif()
    # this candidate build dir does not work, try next one in list
    MATH(EXPR BUILD_CANDIDATE_IDX "${BUILD_CANDIDATE_IDX}+1")
    list(GET LS1_BUILD_DIR_NAMES ${BUILD_CANDIDATE_IDX} LS1_BUILD_DIR_NAME)
  endwhile()
  message(STATUS ">> ls1 build directory used: " ${LS1_SRC_DIR}/${LS1_BUILD_DIR_NAME})
  load_cache(${LS1_SRC_DIR}/${LS1_BUILD_DIR_NAME} READ_WITH_PREFIX LS1_CACHE_ ENABLE_AUTOPAS ENABLE_MPI PRECISION ENABLE_VTK ENABLE_ADIOS2 FIND_PACKAGE_ADIOS2 Eigen3_ForceBundled AUTOPAS_ENABLE_ENERGY_MEASUREMENTS CMAKE_BUILD_TYPE )
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
  find_library(MarDyn REQUIRED NAMES "MarDyn" PATHS ${LS1_SRC_DIR}/${LS1_BUILD_DIR_NAME}/src )
  if(LS1_BUILT_WITH_AUTOPAS)
    if(LS1_CACHE_AUTOPAS_ENABLE_ENERGY_MEASUREMENTS)
      message(FATAL_ERROR "Unsupported flag! Please rebuild ls1 with AUTOPAS_ENABLE_ENERGY_MEASUREMENTS turned off!")
    endif()
    if(NOT BUILD_WITH_EIGEN)
      message(WARNING "AutoPas requires Eigen3 support, but BUILD_WITH_EIGEN not enabled!\nProceeding with Eigen3 enabled only for ls1")
    endif()
    if(NOT LS1_CACHE_Eigen3_ForceBundled)
      find_package(Eigen3 REQUIRED NO_MODULE)
      message(STATUS "Using Eigen3 available in environment for ls1")
    else()
      set(eigen3_SOURCE_DIR "${LS1_SRC_DIR}/${LS1_BUILD_DIR_NAME}/_deps/eigen3-src")
      add_library(Eigen3::Eigen OBJECT IMPORTED GLOBAL)
      target_include_directories(Eigen3::Eigen SYSTEM INTERFACE "${eigen3_SOURCE_DIR}")
      message(STATUS "Using Eigen3 fetched by ls1")
      if(BUILD_WITH_EIGEN)
        message(WARNING "MaMiCo will use Eigen3 bundled with ls1. Double check the library versions/paths!")
      endif()
    endif()
    find_library(autopas REQUIRED NAMES "autopas" PATHS ${LS1_SRC_DIR}/${LS1_BUILD_DIR_NAME}/_deps/autopasfetch-build/src/autopas)
    find_library(molecularDynamicsLibrary REQUIRED NAMES "molecularDynamicsLibrary" PATHS ${LS1_SRC_DIR}/${LS1_BUILD_DIR_NAME}/_deps/autopasfetch-build/applicationLibrary)
    if("${LS1_CACHE_CMAKE_BUILD_TYPE}" STREQUAL "Debug")
      find_library(spdlog REQUIRED NAMES "spdlogd" PATHS ${LS1_SRC_DIR}/${LS1_BUILD_DIR_NAME}/_deps/spdlog-build )
    else()
      find_library(spdlog REQUIRED NAMES "spdlog" PATHS ${LS1_SRC_DIR}/${LS1_BUILD_DIR_NAME}/_deps/spdlog-build)
    endif()
    if(LS1_CACHE_AUTOPAS_ENABLE_RULES_BASED_AND_FUZZY_TUNING)
      find_library(antlr4-runtime REQUIRED NAMES "antlr4-runtime" PATHS ${LS1_SRC_DIR}/${LS1_BUILD_DIR_NAME}/_deps/autopasfetch-build/_deps/antlr4cppPrefix/install/lib)
    endif()
    if(LS1_CACHE_AUTOPAS_ENABLE_HARMONY)
      find_library(harmony REQUIRED NAMES "harmony" PATHS ${LS1_SRC_DIR}/${LS1_BUILD_DIR_NAME}/_deps/autopasfetch-build/libs/harmony/include/lib )
    endif()
  endif()
  if(LS1_CACHE_ENABLE_VTK)
    find_library(VTK_LIB xerces-c HINTS $ENV{XERCES_LIBDIR})
    if(NOT VTK_LIB)
      message(FATAL_ERROR "xerces-c lib not found. Disable VTK support, if you do not need it.")
    endif()
  endif()
  if(BUILD_WITH_ADIOS2 AND NOT LS1_CACHE_ENABLE_ADIOS2)
    message(FATAL_ERROR "Please rebuild ls1 with ENABLE_ADIOS2 turned on, if adios2 is needed, else turn BUILD_WITH_ADIOS2 off")
  endif()
  if(LS1_CACHE_ENABLE_ADIOS2)
    if(NOT BUILD_WITH_ADIOS2)
      message(WARNING "ls1 built with adios2, but BUILD_WITH_ADIOS2 not enabled by MaMiCo!\nProceeding with adios2 enabled only for ls1")
    endif()
    if(LS1_CACHE_FIND_PACKAGE_ADIOS2)
      enable_language(C)
      find_package(ADIOS2 REQUIRED)
      message(STATUS "Using adios2 available in environment for ls1")
    else()
      enable_language(C)
      find_package(ADIOS2 REQUIRED PATHS ${LS1_SRC_DIR}/${LS1_BUILD_DIR_NAME}/_deps/adios2fetch-build NO_DEFAULT_PATH)
      message(STATUS "Using adios2 fetched by ls1")
      if(BUILD_WITH_ADIOS2)
        message(WARNING "Both ls1 and MaMiCo will use adios2. Double check the library versions/paths!")
      endif()
    endif()
  endif()
  message(STATUS ">> ls1 found and initialized")
endif()

if(BUILD_WITH_EIGEN AND NOT LS1_BUILT_WITH_AUTOPAS) #ls1 with autopas enabled can come bundled with own eigen3, hence finding eigen3 earlier for ls1
  find_package(Eigen3 REQUIRED NO_MODULE)
endif()

if(BUILD_TESTING)
  find_package(CppUnit REQUIRED)
endif()


#
# shared macros across targets
#
add_library(shared_defines INTERFACE)
target_compile_definitions(shared_defines INTERFACE ${MD_SIM})
target_compile_definitions(shared_defines INTERFACE KOKKOS_TARGET_${KOKKOS_TARGET})
if(CMAKE_BUILD_TYPE STREQUAL "DebugMD")
  target_compile_definitions(shared_defines INTERFACE MDDebug TarchDebug MDError)
elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
  target_compile_definitions(shared_defines INTERFACE MDError TarchDebug)
elseif(CMAKE_BUILD_TYPE STREQUAL "DebugOptimized")
  target_compile_definitions(shared_defines INTERFACE MDError TarchDebug)
endif()
if(BUILD_WITH_MPI)
    target_compile_definitions(shared_defines INTERFACE MDParallel TarchParallel)
endif()

#
# simplemd interface
#
add_library(simplemd_interface INTERFACE)
target_link_libraries(simplemd_interface INTERFACE shared_defines)
target_compile_definitions(simplemd_interface INTERFACE ${MDDim})
target_compile_options(simplemd_interface INTERFACE -fpic)
if(NOT KOKKOS_TARGET STREQUAL "SERIAL")
  target_compile_options(simplemd_interface INTERFACE -fopenmp)
endif()

if(KOKKOS_TARGET STREQUAL "CUDA")
    set(FLAG_OPT_DEBUG "-O0")
else()
    set(FLAG_OPT_DEBUG "-Og")
endif()

if(CMAKE_BUILD_TYPE STREQUAL "DebugMD")
  target_compile_options(simplemd_interface INTERFACE ${FLAG_OPT_DEBUG} -g3)
  if(KOKKOS_TARGET STREQUAL "CUDA")
    target_compile_options(simplemd_interface INTERFACE --fmad=false)
  endif()
elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
  target_compile_options(simplemd_interface INTERFACE ${FLAG_OPT_DEBUG} -g3)
  if(KOKKOS_TARGET STREQUAL "CUDA")
    target_compile_options(simplemd_interface INTERFACE --fmad=false)
  endif()
elseif(CMAKE_BUILD_TYPE STREQUAL "DebugOptimized")
  # Note: -fexpensive-optimizations breaks SimpleMDBench
  target_compile_options(simplemd_interface INTERFACE -O3 -fno-expensive-optimizations -g3)
  if(KOKKOS_TARGET STREQUAL "CUDA")
    target_compile_options(simplemd_interface INTERFACE --fmad=false)
  endif()
elseif(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
  target_compile_options(simplemd_interface INTERFACE -O3 -g3)
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
  if(KOKKOS_TARGET STREQUAL "CUDA")
    target_compile_options(simplemd_interface INTERFACE -O3)
  else()
    target_compile_options(simplemd_interface INTERFACE -Ofast)
  endif()
else()
  message(FATAL_ERROR "Unsupported build type")
endif()

if(BUILD_WITH_MPI)
    target_link_libraries(simplemd_interface INTERFACE MPI::MPI_CXX)
endif()
target_link_libraries(simplemd_interface INTERFACE Kokkos::kokkos)
if(BUILD_WITH_ADIOS2)
  target_compile_definitions(simplemd_interface INTERFACE BUILD_WITH_ADIOS2)
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
  target_compile_options(simplemd_interface INTERFACE -Wall -Werror -pedantic -pedantic-errors -Wstrict-aliasing -fstrict-aliasing -Wno-long-long -Wno-unknown-pragmas -Wno-non-virtual-dtor)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  target_compile_options(simplemd_interface INTERFACE -fstrict-aliasing -qpack_semantic=gnu -ipo -align -ansi-alias) #-w
  target_link_options(simplemd_interface INTERFACE -ipo)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "FujitsuClang")
  target_compile_options(simplemd_interface INTERFACE -Wall -Werror -pedantic -pedantic-errors -Wstrict-aliasing -fstrict-aliasing -Wno-long-long -Wno-unknown-pragmas -Wno-non-virtual-dtor -Wno-sign-conversion -Wno-unused-private-field -mcpu=a64fx)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  target_compile_options(simplemd_interface INTERFACE -Wall -pedantic -pedantic-errors -Wstrict-aliasing -fstrict-aliasing -Wno-long-long -Wno-unknown-pragmas -Wno-non-virtual-dtor -Wno-sign-conversion -Wno-unused-private-field -Wno-gnu-array-member-paren-init)
else()
  message(STATUS ">> Unsupported compiler.")
endif()

#
# git version info target
#
set(GIT_HEADER ${CMAKE_BINARY_DIR}/generated/mamico_git_version_info.h)
add_custom_command(
    OUTPUT ${GIT_HEADER}
    COMMAND ${CMAKE_COMMAND}
            -DOUTPUT_FILE=${GIT_HEADER}
            -P ${CMAKE_SOURCE_DIR}/cmake/GenerateGitHeader.cmake
    DEPENDS ${CMAKE_SOURCE_DIR}/.git/HEAD
            ${CMAKE_SOURCE_DIR}/.git/index
    COMMENT "Generating git version header"
)

add_custom_target(mamico_git_version_info DEPENDS ${GIT_HEADER})

#
# simplemd library target
#
set(SIMPLE_MD_SOURCES
  "simplemd/BoundaryTreatment.cpp"
  "simplemd/MolecularDynamicsSimulation.cpp"
  "simplemd/ProfilePlotter.cpp"
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
  "simplemd/configurations/DomainDecompConfiguration.cpp"
  "simplemd/molecule-mappings/VelocityStoermerVerletMapping.cpp"
  "simplemd/molecule-mappings/VTKMoleculeWriter.cpp"
  "simplemd/molecule-mappings/WriteCheckPointMapping.cpp"
  "simplemd/services/ExternalForceService.cpp"
  "simplemd/services/MolecularPropertiesService.cpp"
  "simplemd/services/MoleculeService.cpp"
  "simplemd/services/ParallelAndLocalBufferService.cpp"
  "simplemd/services/ParallelTopologyService.cpp"
  "tarch/tinyxml2/tinyxml2.cpp"
  "tarch/utils/RandomNumberService.cpp"
  "simplemd/configurations/Adios2Configuration.cpp"
  "simplemd/MoleculeContainer.cpp"
)
if(BUILD_WITH_ADIOS2)
  list(APPEND SIMPLE_MD_SOURCES "simplemd/molecule-mappings/Adios2Writer.cpp")
endif()


add_library(libsimplemd ${SIMPLE_MD_SOURCES})
target_link_libraries(libsimplemd INTERFACE shared_defines)
target_link_libraries(libsimplemd PRIVATE simplemd_interface)
set_target_properties(libsimplemd PROPERTIES CXX_STANDARD 17)
set_target_properties(libsimplemd PROPERTIES OUTPUT_NAME simplemd)

#
# simplemd executable target
#
add_executable(simplemd ${SIMPLE_MD_SOURCES} "simplemd/main.cpp")
target_link_libraries(simplemd INTERFACE shared_defines)
target_link_libraries(simplemd PRIVATE simplemd_interface)
set_target_properties(simplemd PROPERTIES CXX_STANDARD 17)
set_target_properties(simplemd PROPERTIES OUTPUT_NAME simplemd)

#
# multi-simplemd target
#
add_executable(multi-simplemd ${SIMPLE_MD_SOURCES} "simplemd/main_multi.cpp")
target_link_libraries(multi-simplemd INTERFACE shared_defines)
target_link_libraries(multi-simplemd PRIVATE simplemd_interface)
set_target_properties(multi-simplemd PROPERTIES CXX_STANDARD 17)

#
# MaMiCo interface
#
add_library(mamico INTERFACE)
add_dependencies(mamico mamico_git_version_info)
target_include_directories(mamico INTERFACE ${CMAKE_BINARY_DIR}/generated)
target_link_libraries(mamico INTERFACE shared_defines)
target_link_libraries(mamico INTERFACE Kokkos::kokkos)
if(BUILD_WITH_MPI)
  option(MAMICO_USE_COLLECTIVE_MPI "Use broadcast/reduce between MD instances, default is OFF." OFF)
  if(MAMICO_USE_COLLECTIVE_MPI)
    target_compile_definitions(mamico INTERFACE MAMICO_USE_COLLECTIVE_MPI)
  endif()
endif()
if(BUILD_WITH_EIGEN)
  target_compile_definitions(mamico INTERFACE BUILD_WITH_EIGEN)
  target_link_libraries(mamico INTERFACE Eigen3::Eigen)
endif()
if(BUILD_WITH_ADIOS2)
  target_compile_definitions(mamico INTERFACE BUILD_WITH_ADIOS2)
endif()
target_compile_definitions(mamico INTERFACE MDDim3)
target_compile_options(mamico INTERFACE -Wall -Wno-unknown-warning-option -Wno-unknown-pragmas -Wfatal-errors)
if(NOT KOKKOS_TARGET STREQUAL "SERIAL")
  target_compile_options(mamico INTERFACE -fopenmp)
endif()
if(BUILD_WITH_MPI)
  target_compile_definitions(mamico INTERFACE MDCoupledParallel MPICH_IGNORE_CXX_SEEK)
  target_compile_options(mamico INTERFACE -Werror -Wno-int-in-bool-context)
  if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "FujitsuClang")
    target_compile_options(mamico INTERFACE -Wno-uninitialized)
  else()
    target_compile_options(mamico INTERFACE -Wno-maybe-uninitialized)
  endif()
  target_link_libraries(mamico INTERFACE MPI::MPI_CXX)
endif()
if(BUILD_WITH_OPENFOAM)
  target_compile_definitions(mamico INTERFACE BUILD_WITH_OPENFOAM WM_LABEL_SIZE=32 WM_DP NoRepository)
  target_include_directories(mamico INTERFACE $ENV{FOAM_SRC}/OpenFOAM/lnInclude $ENV{FOAM_SRC}/OSspecific/POSIX/lnInclude $ENV{FOAM_SRC}/finiteVolume/lnInclude $ENV{FOAM_SRC}/meshTools/lnInclude)
  target_link_directories(mamico INTERFACE $ENV{FOAM_LIBBIN} $ENV{FOAM_LIBBIN}/dummy $ENV{FOAM_LIBBIN}/sys-openmpi)
  target_link_libraries(mamico INTERFACE OpenFOAM finiteVolume meshTools Pstream surfMesh fileFormats dl m)
endif()
if(BUILD_WITH_PRECICE)
  target_compile_definitions(mamico INTERFACE BUILD_WITH_PRECICE)
  target_link_libraries(mamico INTERFACE precice::precice)
endif()
if(BUILD_WITH_PINT_DEBUG)
  target_compile_definitions(mamico INTERFACE PINT_DEBUG)
endif()
if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "FujitsuClang" )
  target_compile_options(mamico INTERFACE -Wno-undefined-var-template)
endif()
# we always link against libsimplemd, even if ${MD_SIM} not SIMPLE_MD, so that we can use the simpleMD configuration
target_link_libraries(mamico INTERFACE libsimplemd)
if("${MD_SIM}" STREQUAL "LAMMPS_MD" OR "${MD_SIM}" STREQUAL "LAMMPS_DPD")
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
    target_include_directories(mamico INTERFACE ${LS1_SRC_DIR}/${LS1_BUILD_DIR_NAME}/_deps/autopasfetch-src/src ${LS1_SRC_DIR}/${LS1_BUILD_DIR_NAME}/_deps/spdlog-src/include ${LS1_SRC_DIR}/${LS1_BUILD_DIR_NAME}/_deps/autopasfetch-src/applicationLibrary/molecularDynamics)
    target_link_directories(mamico INTERFACE ${LS1_SRC_DIR}/${LS1_BUILD_DIR_NAME}/_deps/autopasfetch-build/src/autopas ${LS1_SRC_DIR}/${LS1_BUILD_DIR_NAME}/_deps/spdlog-build ${LS1_SRC_DIR}/${LS1_BUILD_DIR_NAME}/_deps/autopasfetch-build/applicationLibrary ${LS1_SRC_DIR}/${LS1_BUILD_DIR_NAME}/_deps/autopasfetch-build/_deps/antlr4cppPrefix/install/lib)
    target_link_libraries(mamico INTERFACE autopas pthread spdlog molecularDynamicsLibrary)
    if(LS1_CACHE_AUTOPAS_ENABLE_HARMONY)
      target_link_directories(mamico INTERFACE ${LS1_SRC_DIR}/${LS1_BUILD_DIR_NAME}/_deps/autopasfetch-build/libs/harmony/include/lib)
      target_link_libraries(mamico INTERFACE harmony)
    endif()
    if(LS1_CACHE_AUTOPAS_ENABLE_RULES_BASED_AND_FUZZY_TUNING)
      target_link_libraries(mamico INTERFACE antlr4-runtime)
    endif()
    #pthread is required by spdlog
    target_compile_definitions(mamico INTERFACE BUILD_WITH_EIGEN)
    target_link_libraries(mamico INTERFACE Eigen3::Eigen)
  endif()
  if(LS1_CACHE_ENABLE_VTK)
    target_compile_definitions(mamico INTERFACE VTK)
    target_link_libraries(mamico INTERFACE ${VTK_LIB})
  endif()
  if(LS1_CACHE_ENABLE_ADIOS2)
    if (NOT DEFINED ADIOS2_HAVE_MPI)
        set(ADIOS2_HAVE_MPI OFF)
    endif()
    if ((ADIOS2_HAVE_MPI OR BUILD_WITH_MPI) AND (NOT (ADIOS2_HAVE_MPI AND BUILD_WITH_MPI))) # handish implementation of xor (cmake's EQUAL does not work)
        message(FATAL_ERROR "
        ADIOS2_HAVE_MPI set to \"${ADIOS2_HAVE_MPI}\" and mamico BUILD_WITH_MPI is set to \"${BUILD_WITH_MPI}\".
        For a parallel build, MPI has to be enabled in the packaged ADIOS2 as well as mamico (ADIOS2_ENABLE_MPI=ON, BUILD_WITH_MPI=ON).
        For a sequential build, deactivate MPI for both.")
    endif()
    target_compile_definitions(mamico INTERFACE ENABLE_ADIOS2)
    target_link_libraries(mamico INTERFACE adios2::adios2)
  endif()
endif()
if(ENABLE_COVERAGE)
  target_compile_options(mamico INTERFACE --coverage)
  target_link_libraries(mamico INTERFACE gcov)
endif()

if(CMAKE_BUILD_TYPE STREQUAL "DebugMD")
  target_compile_definitions(mamico INTERFACE MDCoupledDebug MDCoupledError NLM_DEBUG DEBUG_FILTER_PIPELINE)
  target_compile_options(mamico INTERFACE ${FLAG_OPT_DEBUG} -g3)
elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
  if (NOT ENABLE_COVERAGE)
    target_compile_definitions(mamico INTERFACE MDCoupledDebug MDCoupledError NLM_DEBUG DEBUG_FILTER_PIPELINE PINT_DEBUG)
  endif()
  target_compile_options(mamico INTERFACE ${FLAG_OPT_DEBUG} -g3)
elseif(CMAKE_BUILD_TYPE STREQUAL "DebugOptimized")
  target_compile_definitions(mamico INTERFACE MDCoupledError)
  target_compile_options(mamico INTERFACE -O3 -g3)
elseif(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
  target_compile_options(mamico INTERFACE -O3 -g3)
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
  if(KOKKOS_TARGET STREQUAL "CUDA")
    target_compile_options(mamico INTERFACE -O3)
  else()
    target_compile_options(mamico INTERFACE -Ofast)
  endif()
else()
  message(FATAL_ERROR "Unsupported build type")
endif()

# List of sources for all Mamico targets
set(MAMICO_SOURCES
  "coupling/indexing/IndexingService.cpp"
  "coupling/configurations/ParticleInsertionConfiguration.cpp"
  "coupling/solvers/LBCouetteSolverState.cpp"
)
#
# Couette target
#
set(COUETTE_SOURCES
  ${MAMICO_SOURCES}
  "coupling/scenario/main_couette.cpp"
)
add_executable(couette ${COUETTE_SOURCES})
target_link_libraries(couette INTERFACE shared_defines)
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
    ${MAMICO_SOURCES}
    "tarch/tinyxml2/tinyxml2.cpp"
    "tarch/utils/RandomNumberService.cpp"
    "test/integration/main_lammps.cpp"
  )
  add_executable(test_lammps ${TEST_LAMMPS_SOURCES})
  target_link_libraries(test_lammps INTERFACE shared_defines)
  set_target_properties(test_lammps PROPERTIES CXX_STANDARD 17)
  if (CMAKE_CXX_COMPILER_ID STREQUAL "FujitsuClang")
    target_compile_options(test_lammps PRIVATE -std=c++17)
  endif()
  target_link_libraries(test_lammps PRIVATE mamico)
endif()

#
# pybind target
#
if(BUILD_WITH_PYBIND11)
  set(PYBIND_SOURCES11
    ${MAMICO_SOURCES}
    "coupling/python-binding/mamico.cpp"
  )
  pybind11_add_module(mamico_pybind
    ${PYBIND_SOURCES11}
  )
  target_link_libraries(mamico_pybind INTERFACE shared_defines)
  set_target_properties(mamico_pybind PROPERTIES OUTPUT_NAME mamico)
  target_compile_options(mamico_pybind PRIVATE -fPIC)
  set_target_properties(mamico_pybind PROPERTIES CXX_STANDARD 17)
  target_link_libraries(mamico_pybind PRIVATE mamico)
endif()

#
# Clangformat target
#
set(CLANGFORMAT_DIRS
  coupling
  simplemd
  tarch
  test
)
set(CLANGFORMAT_EXTS
  *.c
  *.h
  *.hpp
  *.cpp
  *.cxx
  *.cpph
)
set(CLANGFORMAT_SOURCES)
foreach(dir IN LISTS CLANGFORMAT_DIRS)
  foreach(ext IN LISTS CLANGFORMAT_EXTS)
    file(GLOB_RECURSE _tmp
      CONFIGURE_DEPENDS
      FOLLOW_SYMLINKS
      "${CMAKE_SOURCE_DIR}/${dir}/${ext}"
    )
    list(APPEND CLANGFORMAT_SOURCES ${_tmp})
  endforeach()
endforeach()
list(REMOVE_DUPLICATES CLANGFORMAT_SOURCES)
foreach(clangformat_source ${CLANGFORMAT_SOURCES})
  get_filename_component(clangformat_source ${clangformat_source} ABSOLUTE)
  list(APPEND CLANGFORMAT_SOURCES_ABS ${clangformat_source})
endforeach()

# Require specific clang-format version
execute_process(
  COMMAND clang-format --version
  OUTPUT_VARIABLE CLANGFORMAT_VERSION_OUTPUT
  ERROR_VARIABLE CLANGFORMAT_VERSION_ERROR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
if(CLANGFORMAT_VERSION_OUTPUT STREQUAL "")
  set(CLANGFORMAT_VERSION "<not installed>")
else()
  string(REGEX MATCH "[0-9]+\\.[0-9]+\\.[0-9]+" CLANGFORMAT_VERSION ${CLANGFORMAT_VERSION_OUTPUT})
  message(STATUS "Found clang-format version: ${CLANGFORMAT_VERSION}")
endif()
if(NOT CLANGFORMAT_VERSION MATCHES "^18\\.1\\.8")
  message(WARNING "clang-format version must be 18.1.8 but found ${CLANGFORMAT_VERSION} (Target \"clangformat\" will not be generated)")
else()
  add_custom_target(
    clangformat
    COMMAND /usr/bin/env clang-format --verbose
    -i
    ${CLANGFORMAT_SOURCES_ABS}
  )
endif()

#
# compiletest target
#
get_filename_component(COMPILETEST_SCRIPT "cmake/CompileTest.sh" ABSOLUTE)
add_custom_target(
  compiletest
  COMMAND ${COMPILETEST_SCRIPT}
  ${COMPILETEST_COMPILERS}
  ${COMPILETEST_MODES}
)

#
# Benchmarks target
#
set(BENCHMARKS_SOURCES
  ${MAMICO_SOURCES}
  "test/integration/benchmarks.cpp"
)
add_executable(benchmarks ${BENCHMARKS_SOURCES})
target_link_libraries(benchmarks INTERFACE shared_defines)
set_target_properties(benchmarks PROPERTIES CXX_STANDARD 17)
target_link_libraries(benchmarks PRIVATE mamico)
target_compile_definitions(benchmarks PRIVATE TarchDebug)

#
# Coverage
#
if(ENABLE_COVERAGE)
  # find required tools
  find_program(LCOV lcov REQUIRED)
  find_program(GENHTML genhtml REQUIRED)
  find_package(Python3 REQUIRED)

  # add coverage target
  add_custom_target(
    coverage
    # gather data
    COMMAND
      ${LCOV} --directory . --capture --exclude */test/* --exclude '/usr/*' --exclude '*/simplemd/*' --output-file coverage.info
    # generate report
    COMMAND ${GENHTML} --rc check_data_consistency=0 --demangle-cpp -o coverage coverage.info
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

    # add my-coverage target
    add_custom_target(
      my-coverage
      COMMAND ${PYTHON3} ../test/check_my_coverage.py --build-directory . --skip-make-coverage
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

    add_dependencies(my-coverage coverage)
endif()

#
# Testing
#
if (BUILD_TESTING)    
  file(GLOB_RECURSE TEST_SOURCES_TOP "test/*Test.cpp")
  list(FILTER TEST_SOURCES_TOP EXCLUDE REGEX "test/unit/coupling/interface/impl/") # remove tests corresponding to specific MD solvers
  
  if("${MD_SIM}" STREQUAL "SIMPLE_MD")
    file(GLOB_RECURSE TEST_SOURCES_MD "test/unit/coupling/interface/impl/SimpleMD/*Test.cpp")
  elseif("${MD_SIM}" STREQUAL "LS1_MARDYN")
    file(GLOB_RECURSE TEST_SOURCES_MD "test/unit/coupling/interface/impl/ls1/*Test.cpp")
    list(FILTER TEST_SOURCES_TOP EXCLUDE REGEX "test/unit/coupling/InstanceHandlingTest.cpp")
  elseif("${MD_SIM}" STREQUAL "LAMMPS_MD" OR "${MD_SIM}" STREQUAL "LAMMPS_DPD")
    file(GLOB_RECURSE TEST_SOURCES_MD "test/unit/coupling/interface/impl/LAMMPS/*Test.cpp")
  else()
    message(FATAL_ERROR "Could not determine the MD solver while compiling tests")
  endif()

  set(TEST_SOURCES ${TEST_SOURCES_TOP} ${TEST_SOURCES_MD})
  if(GET_TEST_NAMES)
    add_executable(getTestNames
      ${MAMICO_SOURCES}
      "test/unit/getTestNames.cpp"
      ${TEST_SOURCES}
    )
    target_link_libraries(getTestNames INTERFACE shared_defines)
    set_target_properties(getTestNames PROPERTIES CXX_STANDARD 17)
    target_link_libraries(getTestNames PRIVATE mamico)
    target_include_directories(getTestNames PRIVATE ${CPPUNIT_INCLUDE_DIRS})
    target_link_directories(getTestNames PRIVATE ${CPPUNIT_LIBRARY})
    target_link_libraries(getTestNames PRIVATE cppunit)
    unset(GET_TEST_NAMES CACHE)
    return()
  endif()

  add_executable(testmamico 
    ${MAMICO_SOURCES} 
    "test/unit/main.cpp"
    ${TEST_SOURCES}
  )
  target_link_libraries(testmamico INTERFACE shared_defines)
  set_target_properties(testmamico PROPERTIES CXX_STANDARD 17)
  target_compile_definitions(testmamico PRIVATE MDCoupledError MDError)
  target_link_libraries(testmamico PRIVATE mamico)
  target_include_directories(testmamico PRIVATE ${CPPUNIT_INCLUDE_DIRS})
  target_link_directories(testmamico PRIVATE ${CPPUNIT_LIBRARY})
  target_link_libraries(testmamico PRIVATE cppunit)

  message(STATUS ">> Running cmake with GET_TEST_NAMES ON ...")
  message(STATUS "(((((((((((((((((((((((((((((((((( SUBPROCESSES_START ")
  execute_process(COMMAND cmake -DGET_TEST_NAMES=ON -DBUILD_TESTING=ON -DENABLE_COVERAGE=OFF -DBUILD_WITH_MPI=${BUILD_WITH_MPI} .)
  message(STATUS ">> Building getTestNames executable ...")
  execute_process(COMMAND make -j4 getTestNames)
  message(STATUS "SUBPROCESSES_END ))))))))))))))))))))))))))))))))))")
  execute_process(COMMAND ./getTestNames OUTPUT_VARIABLE TEST_NAMES)
  separate_arguments(TEST_NAMES UNIX_COMMAND ${TEST_NAMES})

  foreach(TEST_NAME ${TEST_NAMES})
    if (BUILD_WITH_MPI)
    # We have to use the --oversubscribe OpenMPI option here as openmpi is installed in the test container (see test/docker/Dockerfile) 
    # and the  the github runner VM has less than 4 VCPUS.
      add_test(NAME ${TEST_NAME} COMMAND mpirun --oversubscribe -np 4 ./testmamico ${TEST_NAME})
    else()
      add_test(NAME ${TEST_NAME} COMMAND ./testmamico ${TEST_NAME})
    endif()
  endforeach()
endif()
