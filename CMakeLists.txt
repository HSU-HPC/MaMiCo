cmake_minimum_required(VERSION 3.13)

project ("mamico")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(MPI)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

add_library(mamico INTERFACE)

target_include_directories(
    mamico
    INTERFACE   "./")

target_compile_definitions(mamico INTERFACE "SIMPLE_MD")

file(GLOB SOURCES
    "coupling/cell-mappings/*cpp"
    "coupling/configurations/*cpp"
    "coupling/datastructures/*cpp"
    "coupling/filtering/*cpp"
    "coupling/interface/*cpp"
    "coupling/noisereduction/*cpp"
    "coupling/paralleltopology/*cpp"
    "coupling/sendrecv/*cpp"
    "coupling/services/*cpp"
    "coupling/solvers/CoupledMolecularDynamicsSimulation.cpp"
    "coupling/transferstrategies/*cpp"
    "simplemd/BoundaryTreatment.cpp"
    "simplemd/MolecularDynamicsSimulation.cpp"
    "simplemd/ProfilePlotter.cpp"
    "simplemd/*/*.cpp"
    "tarch/tinyxml2/*.cpp"
    "tarch/configuration/*.cpp"
    "tarch/la/*cpp"
    "tarch/utils/*cpp")

add_library(simpleMD "simplemd/main_multi.cpp"
                        "simplemd/BoundaryTreatment.cpp"
                        "simplemd/MolecularDynamicsSimulation.cpp"
                        "simplemd/ProfilePlotter.cpp"
                        ${SOURCES})

target_include_directories(simpleMD PUBLIC "./")
#target_link_directories(simpleMD PUBLIC "./")
target_compile_definitions(simpleMD PUBLIC "MDDim3")
target_link_libraries(simpleMD Eigen3::Eigen)

add_executable(dynamicMDTest "coupling/tests/main_dynamic_md.cpp")
                             #"${SOURCES}")

target_include_directories(dynamicMDTest PUBLIC "./" 
                                                "coupling/interface/" 
                                                "coupling/interface/impl/SimpleMD"
                                                "${MPI_CXX_INCLUDE_DIRS}")

target_link_libraries(dynamicMDTest Eigen3::Eigen)
target_link_libraries(dynamicMDTest ${MPI_CXX_LIBRARIES})
target_link_libraries(dynamicMDTest mamico)
target_link_libraries(dynamicMDTest simpleMD)

target_compile_definitions(dynamicMDTest PUBLIC "-DMDDim3 -pedantic -Werror -Wall \
                                                 -Wno-unknown-pragmas -DSIMPLE_MD \
                                                 -DMDCoupledParallel -DTarchParallel \
                                                 -DMPICH_IGNORE_CXX_SEEK -O0 -g3 ${MPI_CXX_FLAGS}")