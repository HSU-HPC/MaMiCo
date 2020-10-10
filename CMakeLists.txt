cmake_minimum_required(VERSION 3.13)

project ("mamico")

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_library(mamico INTERFACE)

target_include_directories(
    mamico
    INTERFACE   "coupling/services"
                "tarch/configuration"
)


file(GLOB SOURCES
    "simplemd/*/*.cpp"
    "tarch/tinyxml2/*.cpp"
    "tarch/configuration/*.cpp"
    "tarch/la/*cpp"
    "tarch/utils/*cpp")

add_executable(simpleMD "simplemd/main_multi.cpp" 
                        "simplemd/BoundaryTreatment.cpp"
                        "simplemd/MolecularDynamicsSimulation.cpp"
                        "simplemd/ProfilePlotter.cpp"
                        ${SOURCES})


target_include_directories(simpleMD PUBLIC "./")
target_link_directories(simpleMD PUBLIC "./")

target_compile_definitions(simpleMD PUBLIC "MDDim3")