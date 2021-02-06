CMAKE_MINIMUM_REQUIRED (VERSION 3.19) 
PROJECT (MAMICO_CTEST LANGUAGES CXX)

find_package(Python COMPONENTS Interpreter Development REQUIRED)
find_package(pybind11 REQUIRED)

set ( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR})
include_directories ( ${CMAKE_CURRENT_SOURCE_DIR} )

#--- Define target ---

#Case: standard C++ executable
#add_library(mamico_ctest
#	${CMAKE_CURRENT_SOURCE_DIR}/coupling/solvers/CoupledMolecularDynamicsSimulation.cpp
#	${CMAKE_CURRENT_SOURCE_DIR}/coupling/configurations/ParticleInsertionConfiguration.cpp
#	[...]main.cpp
#	Add Simplemd to libraries?
#)

#Case: pybind11
pybind11_add_module( mamico_ctest
	${CMAKE_CURRENT_SOURCE_DIR}/coupling/python-binding/mamico.cpp
)

#enforce c++17 in order to support the use of constexpr
target_compile_features(mamico_ctest PRIVATE cxx_std_17)

#Add SimpleMD's scons build script to target ---
add_custom_target(
	TARGET mamico_ctest
	#TODO: make this customizable? (dim, parallel...)
	COMMAND cd ${CMAKE_CURRENT_SOURCE_DIR} && scons target=libsimplemd dim=3 build=release parallel=yes -j4 || exit $?
)

#set simplemd macros
set_property(
	TARGET mamico_ctest
	PROPERTY COMPILE_DEFINITIONS MDDim3 MDError SIMPLE_MD
	)

#For pybind, this is always true TODO: Support standard C++ versions, e.g. Couette
set(MAMICO_PARALLEL TRUE)
if(MAMICO_PARALLEL)
	set_property(
		TARGET mamico_ctest
   		APPEND PROPERTY COMPILE_DEFINITIONS MDParallel Parallel MDCoupledParallel TarchParallel
   	)
endif()

#If thesee include directories do not suffice, you can add more manually via -DCMAKE_CXX_FLAGS
target_include_directories (mamico_ctest PUBLIC /usr/include/ )
target_include_directories (mamico_ctest PUBLIC /usr/lib/ )


