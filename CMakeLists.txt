CMAKE_MINIMUM_REQUIRED (VERSION 2.6) 
PROJECT ( mamico ) 
käs(e
set ( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${mamico_SOURCE_DIR})

if ( NOT CMAKE_BUILD_TYPE )
    set ( CMAKE_BUILD_TYPE Release CACHE STRING "Build Types: Debug Release"  FORCE )
endif()
SET_PROPERTY( CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release)
               
include_directories ( ${CMAKE_CURRENT_SOURCE_DIR} )
                 
FILE(GLOB_RECURSE MAMICOFILES "*.cpp")
FILE(GLOB_RECURSE to_remove "coupling/tests/*.cpp")
FILE(GLOB_RECURSE to_remove2 "coupling/interface/impl/LAMMPS/*.cpp")
FILE(GLOB_RECURSE to_remove3 "coupling/interface/impl/Espresso/*.cpp")
FILE(GLOB_RECURSE to_remove4 "coupling/interface/impl/MarDyn/*.cpp")
FILE(GLOB_RECURSE to_remove5 "coupling/solvers/exahype/*.cpp")
list(REMOVE_ITEM MAMICOFILES ${to_remove} ${to_remove2} ${to_remove3} ${to_remove4} ${to_remove5})
add_library(MAMICO??? ${MAMICOFILES})

set_property(
   	TARGET MAMICO
   	PROPERTY COMPILE_DEFINITIONS MDDim3 MDError
)

if(MAMICO_PARALLEL)
	set_property(
   		TARGET MAMICO
   		APPEND PROPERTY COMPILE_DEFINITIONS MDParallel Parallel MDCoupledParallel
   	)
endif()


add_subdirectory( coupling )
add_subdirectory( tarch )
add_subdirectory( simplemd )
