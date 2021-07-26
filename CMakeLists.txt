cmake_minimum_required(VERSION 3.13)
project ("mamico")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Debug
			CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel ..." FORCE)
ENDIF()
#### Options ####
# build with MPI
option(PARALLEL ON)
# Set MD Simulation to use
set(MD_SIM SIMPLE_MD CACHE STRING "Choose one of SIMPLE_MD, LAMMPS_MD, LAMMPS_DPD")
set(MD_DIM MDDim3 CACHE STRING "Set to MDDim3 (default) or MDDim2")
#### Loag required packages ####
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(MPI REQUIRED)
#### MaMiCo Target ####
add_library(mamico INTERFACE)
target_link_libraries(mamico INTERFACE Eigen3::Eigen)
target_include_directories(
		mamico
		INTERFACE   "./")
target_compile_definitions(mamico INTERFACE ${MD_SIM} MDCoupled)
target_compile_options(mamico INTERFACE -pedantic -Werror -Wall -Wno-unknown-pragmas)
#### *.CPP Globs ####
file(GLOB MAMICO_SOURCES
		"coupling/configurations/ParticleInsertionConfiguration.cpp"
		"coupling/solvers/CoupledMolecularDynamicsSimulation.cpp")
file(GLOB SIMPLE_MD_SRC
		"simplemd/BoundaryTreatment.cpp"
		"simplemd/MolecularDynamicsSimulation.cpp"
		"simplemd/ProfilePlotter.cpp"
		"simplemd/*/*.cpp"
		"tarch/tinyxml2/*.cpp"
		"tarch/configuration/*.cpp"
		"tarch/la/*cpp"
		"tarch/utils/*cpp")
#### Simple MD Target ####
add_library(simpleMD ${SIMPLE_MD_SRC})
target_include_directories(simpleMD PUBLIC "./")
target_compile_definitions(simpleMD PUBLIC ${MD_DIM})
target_compile_options(simpleMD PUBLIC -pedantic -Werror -Wall -Wno-unknown-pragmas)
target_link_libraries(mamico INTERFACE Eigen3::Eigen)
#### Default Compiler Flags ####
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic -Werror -Wall -Wno-unknown-pragmas")
if("${MD_SIM}" STREQUAL "SIMPLE_MD")
	set(LINKED_MD simpleMD)
elseif("${MD_SIM}" STREQUAL "LAMMPS_MD")
	set(USER-MAMICO_SOURCES_DIR "coupling/interface/impl/LAMMPS/USER-MAMICO")
	set(PARALLEL ON)
	set(MD_DIM "MDDim3")
	include(FetchContent)
	FetchContent_Declare(
			lammps
			GIT_REPOSITORY https://github.com/lammps/lammps.git
			GIT_TAG        stable_3Mar2020
			SOURCE_SUBDIR  cmake
	)
	if(NOT lammps_POPULATED)
		FetchContent_populate(lammps)
		set(BUILD_EXE OFF CACHE INTERNAL "lammps build exe")
		set(BUILD_LIB ON CACHE INTERNAL "lammps build lib")
		#set(BUILD_OMP ON CACHE INTERNAL "lammps enable openmp")
		set(PACKAGES asphere body class2 colloid dipole fld gpu granular kim
				kspace manybody mc meam misc molecule opt peri poems reax replica
				rigid shock srd voronoi xtc
				)
		set(PACKUSER user-atc user-awpmd user-cg-cmm user-colvars
				user-cuda user-eff user-lb user-misc user-omp user-molfile
				user-phonon user-reaxc user-sph
				)
		set(PACKLIB gpu kim meam poems reax voronoi
				user-atc user-awpmd user-colvars user-cuda user-molfile)
		foreach(PKG ${PACKAGES} ${PACKUSER} ${PACKLIB})
			string(TOUPPER "${PKG}" PKG)
			message(STATUS ${PKG})
			set(PKG_${PKG} ON)
		endforeach()
		add_subdirectory(${lammps_SOURCE_DIR}/cmake ${lammps_BINARY_DIR})
	endif()
	RegisterStyles(${USER-MAMICO_SOURCES_DIR})
	target_include_directories(mamico INTERFACE
			"${USER-MAMICO_SOURCES_DIR}"
			"${lammps_SOURCE_DIR}/src"
			)
	RegisterPackages(${USER-MAMICO_SOURCES_DIR})
	target_sources(lammps PRIVATE "${USER-MAMICO_SOURCES_DIR}/fix_mamico.cpp")
	target_include_directories(lammps PRIVATE "./")
	target_link_libraries(lammps Eigen3::Eigen)
	set(LINKED_MD lammps)
	option(PKG_USER-MAMICO ON)
	include("${USER-MAMICO_SOURCES_DIR}/PKG_USER-MAMICO.cmake")
endif()
IF(PARALLEL)

	target_compile_definitions(mamico INTERFACE MDCoupledParallel )
	target_compile_definitions(simpleMD PUBLIC MDParallel TarchParallel)

	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_FLAGS}")

	target_link_libraries(simpleMD ${MPI_CXX_LIBRARIES})
	target_include_directories(simpleMD PUBLIC "${MPI_CXX_INCLUDE_DIRS}")
	target_include_directories(mamico INTERFACE "${MPI_CXX_INCLUDE_DIRS}")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_CXX_LINK_FLAGS}")
ENDIF()
#### TARGETS ####
## Dynamic MD Test
add_executable(test "coupling/tests/main_couette.cpp"
		"${MAMICO_SOURCES}"
		"${SIMPLE_MD_SRC}")
target_link_libraries(test mamico)
target_link_libraries(test ${LINKED_MD})
target_compile_definitions(test PUBLIC ${MD_DIM} ${MD_SIM})
target_compile_options(test PUBLIC -pedantic -Werror -Wall -Wno-unknown-pragmas)
if(PARALLEL)
	target_compile_definitions(test PUBLIC TarchParallel MPICH_IGNORE_CXX_SEEK)
	target_link_libraries(test ${MPI_CXX_LIBRARIES})
endif()