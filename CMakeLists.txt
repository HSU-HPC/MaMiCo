CMAKE_MINIMUM_REQUIRED (VERSION 3.19) 
PROJECT ( mamico ) 

set ( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${mamico_SOURCE_DIR})

if ( NOT CMAKE_BUILD_TYPE )
    set ( CMAKE_BUILD_TYPE Release CACHE STRING "Build Types: Debug Release"  FORCE )
endif()
SET_PROPERTY( CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release)

#TODO: move to CompilerTest.sh
#Set compiler
set(CMAKE_CXX_COMPILER "/usr/bin/clang++")


#Compile SimpleMD using scons script
add_custom_target(
	TARGET MAMICO_PYBIND
	#TODO: make this customizable? (dim, parallel...)
	COMMAND cd ${mamico_SOURCE_DIR} && scons target=libsimplemd dim=3 build=release parallel=yes -j4 || exit $?
	)
#TODO: Add simplemd to libraries (add_libraries?)

#TODO: fix this
include_directories ( ${CMAKE_CURRENT_SOURCE_DIR} )

#TODO: use personal_settings or rather CMakeCache
include_directories ( /usr/include/python3.8 )
include_directories ( /usr/lib/openmpi/include )
include_directories ( /usr/include/eigen3 )
                 
#define target
add_library(MAMICO_PYBIND 
	${mamico_SOURCE_DIR}/coupling/solvers/CoupledMolecularDynamicsSimulation.cpp
	${mamico_SOURCE_DIR}/coupling/configurations/ParticleInsertionConfiguration.cpp
	${mamico_SOURCE_DIR}/coupling/python-binding/mamico.cpp
)

set_property(
	TARGET MAMICO_PYBIND
	PROPERTY COMPILE_DEFINITIONS MDDim3 MDError SIMPLE_MD
	)

#enforce c++17 in order to support the use of constexxr
target_compile_features(MAMICO_PYBIND PRIVATE cxx_std_17)

#if(MAMICO_PARALLEL)
	set_property(
		TARGET MAMICO_PYBIND
   		APPEND PROPERTY COMPILE_DEFINITIONS MDParallel Parallel MDCoupledParallel TarchParallel
   	)
#endif()


#add_subdirectory( coupling )
#add_subdirectory( tarch )
#add_subdirectory( simplemd )

