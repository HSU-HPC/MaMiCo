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
# Enable the MPI library for parallel build, default is OFF
option(BUILD_WITH_MPI "Enable the MPI library for parallel build, default is OFF." OFF)
message(STATUS ">> Build with MPI: " ${BUILD_WITH_MPI})
# Dimension for MD build, default is MDDim3
set(MD_DIM MDDim3 CACHE STRING "Set to MDDim3 or MDDim2, default is MDDim3.")
message(STATUS ">> SimpleMD dimension variable : " ${MD_DIM})
mark_as_advanced(MD_DIM)
set_property(CACHE MD_DIM PROPERTY STRINGS MDDim3 MDDim2)
# Enable the Eigen library used (so far) only for POD filtering, default is off
option(BUILD_WITH_EIGEN "Enable the Eigen library used (so far) only for POD filtering, default is off." OFF)
message(STATUS ">> Build with Eigen : " ${BUILD_WITH_EIGEN})
# Enable MaMiCo python bindings
option(BUILD_WITH_PYBIND11 "Enable the MaMiCo python binding, default is off." OFF)
message(STATUS ">> Build with pybind11 : " ${BUILD_WITH_PYBIND11})
if (BUILD_WITH_PYBIND11 AND NOT BUILD_WITH_MPI)
	message(STATUS ">> Build with pybind11 was enabled but build with MPI was not. Enabling build with MPI.")
	set(BUILD_WITH_MPI "ON")
	message(STATUS ">> Build with MPI: " ${BUILD_WITH_MPI})
endif()
# Enable coupling the OpenFOAM icoFOAM API, default is off
option(BUILD_WITH_ICOFOAM "Enable coupling the OpenFOAM icoFOAM API, default is off." OFF)
message(STATUS ">> Build with OpenFOAM icoFOAM : " ${BUILD_WITH_ICOFOAM})
# Enable coupling the OpenFOAM rhoCentralFOAM API, default is off
option(BUILD_WITH_RHOCENTRALFOAM "Enable coupling the OpenFOAM rhoCentralFOAM API, default is off." OFF)
message(STATUS ">> Build with OpenFOAM rhoCentralFOAM : " ${BUILD_WITH_RHOCENTRALFOAM})
# Enable coupling with the preCICE API, default is off
option(BUILD_WITH_PRECICE "Enable coupling with the preCICE API, default is off." OFF)
message(STATUS ">> Build with preCICE : " ${BUILD_WITH_PRECICE})

#
# Find dependencies
#
if(BUILD_WITH_EIGEN)
	find_package(Eigen3 REQUIRED)
endif()

if(BUILD_WITH_MPI)
	find_package(MPI REQUIRED)
endif()

if(BUILD_WITH_PYBIND11)
	find_package(Python COMPONENTS Interpreter Development REQUIRED)
	find_package(pybind11 REQUIRED)
endif()

if(BUILD_WITH_ICOFOAM OR BUILD_WITH_RHOCENTRALFOAM)
	# Check valid OpenFOAM
	if(DEFINED ENV{WM_PROJECT_DIR})
		MESSAGE(STATUS ">> OpenFOAM found in : " $ENV{WM_PROJECT_DIR})
	else()
		message(FATAL_ERROR "The OpenFOAM bashrc is not sourced or OpenFOAM is not installed")
	endif()
endif()

if(BUILD_WITH_PRECICE)
	find_package(precice REQUIRED)
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
set_target_properties(simplemd PROPERTIES CXX_STANDARD 17)
target_compile_definitions(simplemd PRIVATE ${MD_DIM})
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
	target_compile_definitions(simplemd PRIVATE MDError)
endif()
if(BUILD_WITH_MPI)
    target_compile_definitions(simplemd PRIVATE MDParallel TarchParallel)
    target_link_libraries(simplemd PRIVATE MPI::MPI_CXX)
endif()
# Check the CXX compiler ID. IDs are defined here: https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_COMPILER_ID.html
message(STATUS ">> CMAKE_CXX_COMPILER : " ${CMAKE_CXX_COMPILER})
message(STATUS ">> CMAKE_CXX_COMPILER_ID : " ${CMAKE_CXX_COMPILER_ID})
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
	target_compile_options(simplemd PRIVATE -Wall -Werror -pedantic -pedantic-errors -Wstrict-aliasing -fstrict-aliasing -Wno-long-long -Wno-unknown-pragmas -Wconversion -Wno-non-virtual-dtor -O3)
	if (CMAKE_BUILD_TYPE STREQUAL "Debug")
		target_compile_options(simplemd PRIVATE -g3)
	endif()
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
	target_compile_options(simplemd PRIVATE -fstrict-aliasing -qpack_semantic=gnu -ipo)
	target_link_options(simplemd PRIVATE -ipo)
	if (CMAKE_BUILD_TYPE PRIVATE "Debug")
		target_compile_options(simplemd PRIVATE -O3)
	else()
		target_compile_options(simplemd PRIVATE -fast -w -Werror-all -align -ansi-alias)
	endif()
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Fujitsu")
	target_compile_options(simplemd PRIVATE -Nclang -Wall -Werror -pedantic -pedantic-errors -Wstrict-aliasing -fstrict-aliasing -Wno-long-long -Wno-unknown-pragmas -Wconversion -Wno-non-virtual-dtor -Wno-sign-conversion -Wno-unused-private-field -mcpu=a64fx)
	target_link_options(simplemd PRIVATE -ipo)
	if (CMAKE_BUILD_TYPE STREQUAL "Debug")
		target_compile_options(simplemd PRIVATE -g3 -O3)
	else()
		target_compile_options(simplemd PRIVATE -Ofast)
	endif()
else()
	message(STATUS ">> Unsupported compiler.")
endif()

#
# MaMiCo target
#
add_library(mamico INTERFACE)
set_target_properties(mamico PROPERTIES CXX_STANDARD 17)
target_link_libraries(mamico INTERFACE simplemd)
if(BUILD_WITH_EIGEN)
	target_compile_definitions(mamico INTERFACE BUILD_WITH_EIGEN)
	target_link_libraries(mamico INTERFACE Eigen3::Eigen)
endif()
target_compile_definitions(mamico INTERFACE SIMPLE_MD MDDim3)
target_compile_options(mamico INTERFACE -Wall -Wno-unknown-pragmas -Wfatal-errors -O3)
if(BUILD_WITH_MPI)
    target_compile_definitions(mamico INTERFACE MDCoupledParallel TarchParallel MPICH_IGNORE_CXX_SEEK ENABLE_POST_MULTI_INSTANCE_FILTERING)
    #-Werror
    target_compile_options(mamico INTERFACE -Wno-int-in-bool-context -Wno-maybe-uninitialized)
    target_link_libraries(mamico INTERFACE MPI::MPI_CXX)
endif()
if(BUILD_WITH_ICOFOAM)
	target_compile_definitions(mamico INTERFACE BUILD_WITH_OPENFOAM WM_LABEL_SIZE=32 WM_DP NoRepository)
	target_include_directories(mamico INTERFACE $ENV{FOAM_SRC}/OpenFOAM/lnInclude $ENV{FOAM_SRC}/OSspecific/POSIX/lnInclude $ENV{FOAM_SRC}/finiteVolume/lnInclude $ENV{FOAM_SRC}/meshTools/lnInclude)
	target_link_directories(mamico INTERFACE $ENV{FOAM_LIBBIN} $ENV{FOAM_LIBBIN}/dummy $ENV{FOAM_LIBBIN}/sys-openmpi)
	target_link_libraries(mamico INTERFACE OpenFOAM finiteVolume meshTools Pstream surfMesh fileFormats dl m)
endif()
if(BUILD_WITH_RHOCENTRALFOAM)
	target_compile_definitions(mamico INTERFACE BUILD_WITH_OPENFOAM WM_LABEL_SIZE=32 WM_DP NoRepository)
	target_include_directories(mamico INTERFACE $ENV{FOAM_SRC}/OpenFOAM/lnInclude $ENV{FOAM_SRC}/OSspecific/POSIX/lnInclude
  $ENV{FOAM_SRC}/finiteVolume/lnInclude $ENV{FOAM_SRC}/meshTools/lnInclude $ENV{FOAM_SRC}/dynamicFvMesh/lnInclude $ENV{FOAM_SRC}/thermophysicalModels/basic/lnInclude
  $ENV{FOAM_SRC}/transportModels/compressible/lnInclude $ENV{FOAM_SRC}/TurbulenceModels/turbulenceModels/lnInclude $ENV{FOAM_SRC}/TurbulenceModels/compressible/lnInclude
  $ENV{FOAM_APP}/solvers/compressible/rhoCentralFoam/BCs/rho $ENV{FOAM_APP}/solvers/compressible/rhoCentralFoam)
	target_link_directories(mamico INTERFACE $ENV{FOAM_LIBBIN} $ENV{FOAM_LIBBIN}/dummy $ENV{FOAM_LIBBIN}/sys-openmpi)
	target_link_libraries(mamico INTERFACE OpenFOAM finiteVolume meshTools Pstream surfMesh fileFormats dl m fvOptions rhoCentralFoam
  compressibleTransportModels compressibleTurbulenceModels turbulenceModels fluidThermophysicalModels topoChangerFvMesh dynamicFvMesh specie)
endif()
if(BUILD_WITH_PRECICE)
	target_compile_definitions(mamico INTERFACE BUILD_WITH_PRECICE)
	target_link_libraries(mamico INTERFACE precice::precice)
endif()

#
# Couette target
#
set(COUETTE_SOURCES
	"coupling/configurations/ParticleInsertionConfiguration.cpp"
  "coupling/solvers/CoupledMolecularDynamicsSimulation.cpp"
  "coupling/tests/main_couette.cpp"
)
add_executable(couette ${COUETTE_SOURCES})
set_target_properties(couette PROPERTIES CXX_STANDARD 17)
target_link_libraries(couette PUBLIC mamico)

#
# Evaporation target
#
set(EVAPORATION_SOURCES
	"coupling/configurations/ParticleInsertionConfiguration.cpp"
  "coupling/solvers/CoupledMolecularDynamicsSimulation.cpp"
  "coupling/tests/main_evaporation.cpp"
)
add_executable(evaporation ${EVAPORATION_SOURCES})
set_target_properties(evaporation PROPERTIES CXX_STANDARD 17)
target_link_libraries(evaporation PUBLIC mamico)

#
# pybind target
#
if(BUILD_WITH_PYBIND11)
	set(PYBIND_SOURCES11
			"coupling/python-binding/mamico.cpp"
	)
	pybind11_add_module(pybind11
		${PYBIND_SOURCES11}
	)
	target_link_libraries(pybind11 PUBLIC mamico)
endif()

#
# Clangformat target
#
set(CLANGFORMAT_SOURCES
  "tarch/utils/Uncopyable.h"
  "tarch/utils/MultiMDService.cpph"
  "tarch/utils/RandomNumberService.h"
  "tarch/utils/RandomNumberService.cpp"
  "tarch/utils/MultiMDService.h"
  "tarch/configuration/ParseConfiguration.h"
  "tarch/configuration/Configuration.h"
  "tarch/tests/main_readconfig.cpp"
  "tarch/la/VectorOperations.cpph"
  "tarch/la/Vector.h"
  "tarch/la/ScalarOperations.h"
  "tarch/la/Matrix.h"
  "tarch/TarchDefinitions.h"
  "coupling/services/MultiMDCellService.h"
  "coupling/services/MacroscopicCellServiceMacroOnly.h"
  "coupling/services/MacroscopicCellService.cpph"
  "coupling/services/MacroscopicCellTraversalWrappers.cpph"
  "coupling/services/MacroscopicCellService.h"
  "coupling/ErrorEstimation.h"
  "coupling/ParticleInsertion.h"
  "coupling/IndexConversionMD2Macro.cpph"
  "coupling/SetGivenVelocity4MomentumInsertion.h"
  "coupling/MomentumController.h"
  "coupling/filtering/FilterPipeline.cpph"
  "coupling/filtering/FilterPipeline.h"
  "coupling/filtering/SequentialFilter.h"
  "coupling/filtering/sequencing/AsymmetricalFilterJunction.cpph"
  "coupling/filtering/sequencing/FilterJunction.h"
  "coupling/filtering/sequencing/FilterSequence.cpph"
  "coupling/filtering/sequencing/FilterJunction.cpph"
  "coupling/filtering/sequencing/FilterSequence.h"
  "coupling/filtering/sequencing/AsymmetricalFilterJunction.h"
  "coupling/filtering/SequentialFilter.cpph"
  "coupling/filtering/interfaces/FilterInterfaceReadOnly.h"
  "coupling/filtering/interfaces/FilterInterface.h"
  "coupling/filtering/interfaces/AsymmetricalJunctorInterface.h"
  "coupling/filtering/interfaces/JunctorInterface.h"
  "coupling/filtering/filters/NLM.cpph"
  "coupling/filtering/filters/POD.cpph"
  "coupling/filtering/filters/Strouhal.h"
  "coupling/filtering/filters/NLM.h"
  "coupling/filtering/filters/Copy.h"
  "coupling/filtering/filters/Gauss.h"
  "coupling/filtering/filters/Constant.h"
  "coupling/filtering/filters/ReadFromFile.cpph"
  "coupling/filtering/filters/Gauss.cpph"
  "coupling/filtering/filters/FilterFromFunction.h"
  "coupling/filtering/filters/ReadFromFile.h"
  "coupling/filtering/filters/WriteToFile.h"
  "coupling/filtering/filters/Datastructures.h"
  "coupling/filtering/filters/WriteToFileJunctor.h"
  "coupling/filtering/filters/POD.h"
  "coupling/filtering/filters/Strouhal.cpph"
  "coupling/filtering/filters/WriteToFile.cpph"
  "coupling/UsherParticleInsertion.h"
  "coupling/NieVelocityImposition.h"
  "coupling/MultiMDMediator.cpph"
  "coupling/python-binding/mamico.cpp"
  "coupling/python-binding/conversion.h"
  "coupling/NoParticleInsertion.h"
  "coupling/MacroscopicCellPlotter.h"
  "coupling/NoMomentumInsertion.h"
  "coupling/cell-mappings/ComputeKineticEnergyMapping.h"
  "coupling/cell-mappings/MoleculeExtractor.h"
  "coupling/cell-mappings/ComputeMomentumMapping.h"
  "coupling/cell-mappings/VTKMoleculePlotter.h"
  "coupling/cell-mappings/SetKineticEnergyMapping.h"
  "coupling/cell-mappings/ComputeAvgForceAndVelocity.h"
  "coupling/cell-mappings/SetMomentumMapping.h"
  "coupling/cell-mappings/ComputeTemperatureMapping.h"
  "coupling/cell-mappings/ComputeMeanPotentialEnergyMapping.h"
  "coupling/cell-mappings/DeleteParticleMapping.h"
  "coupling/cell-mappings/ComputeVelocityBinMapping.h"
  "coupling/cell-mappings/ComputeMassMapping.h"
  "coupling/cell-mappings/VelocityGradientRelaxationMapping.h"
  "coupling/cell-mappings/PerturbateVelocityMapping.h"
  "coupling/cell-mappings/ComputeTemperatureBinMapping.h"
  "coupling/cell-mappings/SetTemperatureMapping.h"
  "coupling/cell-mappings/NieVelocityImpositionMapping.h"
  "coupling/cell-mappings/ZhouBoundaryForce.h"
  "coupling/ZhouBoundaryForceController.h"
  "coupling/IndexConversionMD2Macro.h"
  "coupling/indexing/CellIndex.cpph"
  "coupling/indexing/IndexingService.cpp"
  "coupling/indexing/Operations.h"
  "coupling/indexing/Testing.cpph"
  "coupling/indexing/Testing.h"
  "coupling/indexing/IndexingService.h"
  "coupling/indexing/CellIndex.h"
  "coupling/sendrecv/FromMacro2MDSendOnly.h"
  "coupling/sendrecv/FromMD2MacroRecvOnly.h"
  "coupling/sendrecv/FromMacro2MD.h"
  "coupling/sendrecv/FromMD2MacroRecvOnly.cpph"
  "coupling/sendrecv/SendReceiveBuffer.h"
  "coupling/sendrecv/FromMacro2MDSendOnly.cpph"
  "coupling/sendrecv/DataExchangeFromMD2Macro.h"
  "coupling/sendrecv/FromMD2Macro.cpph"
  "coupling/sendrecv/DataExchange.h"
  "coupling/sendrecv/DataExchangeFromMacro2MD.h"
  "coupling/sendrecv/SendReceiveBuffer.cpph"
  "coupling/sendrecv/FromMacro2MD.cpph"
  "coupling/sendrecv/FromMD2Macro.h"
  "coupling/IndexConversion.cpph"
  "coupling/UsherParticleInsertion.cpph"
  "coupling/NoBoundaryForce.h"
  "coupling/CouplingMDDefinitions.h"
  "coupling/AdditiveMomentumInsertion.h"
  "coupling/IndexConversion.h"
  "coupling/transferstrategies/AveragingTransferStrategy.h"
  "coupling/transferstrategies/TransferStrategy4NieCoupling.h"
  "coupling/transferstrategies/TransferStrategy4SchwarzCoupling.h"
  "coupling/transferstrategies/TransferStrategy4SchwarzCoupling.cpph"
  "coupling/transferstrategies/TransferStrategy.h"
  "coupling/transferstrategies/DifferenceTransferStrategy.h"
  "coupling/transferstrategies/TransferStrategy4NieCoupling.cpph"
  "coupling/transferstrategies/DirectTransferStrategy.h"
  "coupling/paralleltopology/ParallelTopology.h"
  "coupling/paralleltopology/ZYXTopology.h"
  "coupling/paralleltopology/XYZTopology.h"
  "coupling/paralleltopology/ParallelTopologyFactory.h"
  "coupling/datastructures/MacroscopicCells.h"
  "coupling/datastructures/MacroscopicCell.h"
  "coupling/datastructures/MacroscopicCells.cpph"
  "coupling/datastructures/Molecule.h"
  "coupling/MultiMDMediator.h"
  "coupling/configurations/ParticleInsertionConfiguration.h"
  "coupling/configurations/TransferStrategyConfiguration.h"
  "coupling/configurations/BoundaryForceConfiguration.h"
  "coupling/configurations/CouetteConfiguration.h"
  "coupling/configurations/MaMiCoConfiguration.cpph"
  "coupling/configurations/MomentumInsertionConfiguration.h"
  "coupling/configurations/ParallelTopologyConfiguration.h"
  "coupling/configurations/MacroscopicCellConfiguration.h"
  "coupling/configurations/ParticleInsertionConfiguration.cpp"
  "coupling/configurations/MacroscopicCellConfiguration.cpph"
  "coupling/configurations/MaMiCoConfiguration.h"
  "coupling/configurations/ThermostatConfiguration.h"
  "coupling/VelocityGradientRelaxation.h"
  "coupling/BoundaryForceController.h"
  "coupling/interface/MDSimulationFactory.h"
  "coupling/interface/MoleculeIterator.h"
  "coupling/interface/MacroscopicSolverInterface.h"
  "coupling/interface/impl/SimpleMD/SimpleMDMoleculeIterator.h"
  "coupling/interface/impl/SimpleMD/SimpleMDSolverInterface.h"
  "coupling/interface/impl/SimpleMD/SimpleMDMolecule.h"
  "coupling/interface/impl/macroscopictestsolvers/TestMacroscopicSolverInterface.h"
  "coupling/interface/impl/macroscopictestsolvers/TestMacroscopicSolver.h"
  "coupling/interface/impl/macroscopictestsolvers/MacroscopicSolverFactory.h"
  "coupling/interface/impl/macroscopictestsolvers/VoidMacroscopicSolverInterface.h"
  "coupling/interface/impl/macroscopictestsolvers/Helper.h"
  "coupling/interface/impl/macroscopictestsolvers/UniformFlowSolverInterface.h"
  "coupling/interface/impl/macroscopictestsolvers/VoidMacroscopicSolver.h"
  "coupling/interface/MamicoInterfaceProvider.h"
  "coupling/interface/Molecule.h"
  "coupling/interface/MDSolverInterface.h"
  "coupling/KineticEnergyController.h"
  "coupling/solvers/FiniteDifferenceSolver.h"
  "coupling/solvers/FDCouetteSolver.h"
  "coupling/solvers/IcoFoamInterface.h"
  "coupling/solvers/CoupledMolecularDynamicsSimulation.cpp"
  "coupling/solvers/LBCouetteSolver.h"
  "coupling/solvers/IcoFoamBufferSetup.h"
  "coupling/solvers/CoupledMolecularDynamicsSimulation.h"
  "coupling/solvers/CouetteSolver.h"
  "coupling/solvers/CouetteSolverInterface.h"
  "coupling/solvers/LBCouetteSolverInterface.h"
  "coupling/solvers/NumericalSolver.h"
  "coupling/tests/TestEspressoMDMolecule.h"
  "coupling/tests/Test.h"
  "coupling/tests/MacroscopicCellsTest.h"
  "coupling/tests/TestEspressoAddDeleteMolecule.h"
  "coupling/tests/main_lbcouettesolver.cpp"
  "coupling/tests/TestLammpsCalculateForceEnergy.h"
  "coupling/tests/main_couette.cpp"
  "coupling/tests/EvaporationTest.h"
  "coupling/tests/CouetteTest.h"
  "coupling/tests/FromMD2MacroTest.h"
  "coupling/tests/CompileTest.h"
  "coupling/tests/FromMacro2MDTest.h"
  "coupling/tests/TestEspressoMDSolverInterface.h"
  "coupling/tests/main_lammps.cpp"
  "coupling/tests/TestLammpsMoleculeIterator.h"
  "coupling/tests/TestMarDynMolecule.h"
  "coupling/tests/TestLammpsAddDeleteAtom.h"
  "coupling/tests/IndexConversionTest.h"
  "coupling/tests/main_mardyn.cpp"
  "coupling/tests/TestCell.h"
  "coupling/tests/TestEspressoForceEnergyCalculation.h"
  "coupling/tests/main_evaporation.cpp"
  "coupling/tests/TestLammps.h"
  "coupling/tests/TestLammpsInitialVelocity.h"
  "coupling/tests/main.cpp"
  "coupling/tests/TestLBCouetteSolver.h"
  "coupling/tests/cell_idx_iter_benchmark.cpp"
  "coupling/tests/TestDataExchangeFromMD2Macro.h"
  "coupling/tests/main_espresso.cpp"
  "coupling/tests/TestEspresso.h"
  "coupling/tests/TestMarDyn.h"
  "coupling/tests/TestLammpsGhost.h"
  "coupling/tests/TestMarDynAddDeleteMolecule.h"
  "coupling/tests/TestEspressoMoleculeIterator.h"
  "coupling/tests/TestMarDynMoleculeIterator.h"
  "coupling/MPIIndexConversion.h"
  "coupling/InstanceHandling.h"
  "coupling/MomentumInsertion.h"
  "simplemd/services/ExternalForceService.cpp"
  "simplemd/services/LinkedCellService.h"
  "simplemd/services/ExternalForceService.h"
  "simplemd/services/MoleculeService.h"
  "simplemd/services/ParallelTopologyService.cpp"
  "simplemd/services/MolecularPropertiesService.h"
  "simplemd/services/LinkedCellService.cpp"
  "simplemd/services/MolecularPropertiesService.cpp"
  "simplemd/services/ParallelTopologyService.h"
  "simplemd/services/ParallelAndLocalBufferService.cpp"
  "simplemd/services/ParallelAndLocalBufferService.h"
  "simplemd/services/MoleculeService.cpp"
  "simplemd/MolecularDynamicsDefinitions.h"
  "simplemd/ProfilePlotter.cpp"
  "simplemd/MolecularDynamicsSimulation.cpp"
  "simplemd/BoundaryTreatment.cpp"
  "simplemd/cell-mappings/ResetPotentialEnergyMapping.h"
  "simplemd/cell-mappings/ParallelBoundaryEmptyCellsMapping.h"
  "simplemd/cell-mappings/RDFMapping.h"
  "simplemd/cell-mappings/CollectMoleculesMapping.h"
  "simplemd/cell-mappings/LennardJonesForceMapping.cpp"
  "simplemd/cell-mappings/ProfilePlotterMapping.h"
  "simplemd/cell-mappings/LennardJonesPotentialEnergyMapping.cpp"
  "simplemd/cell-mappings/ComputeTemperatureMapping.h"
  "simplemd/cell-mappings/ComputeMeanVelocityMapping.h"
  "simplemd/cell-mappings/CopyMoleculesMapping.h"
  "simplemd/cell-mappings/VaryCheckpointMapping.h"
  "simplemd/cell-mappings/EmptyLinkedListsMapping.h"
  "simplemd/cell-mappings/DeleteMoleculesMapping.cpp"
  "simplemd/cell-mappings/PeriodicBoundaryEmptyCellsMapping.h"
  "simplemd/cell-mappings/LennardJonesForceMapping.h"
  "simplemd/cell-mappings/RDFMapping.cpp"
  "simplemd/cell-mappings/DeleteMoleculesMapping.h"
  "simplemd/cell-mappings/LennardJonesPotentialEnergyMapping.h"
  "simplemd/cell-mappings/PeriodicAndParallelBoundaryFillCellsMapping.cpp"
  "simplemd/cell-mappings/PeriodicBoundaryEmptyCellsMapping.cpp"
  "simplemd/cell-mappings/PeriodicAndParallelBoundaryFillCellsMapping.h"
  "simplemd/MolecularDynamicsSimulation.h"
  "simplemd/MolecularProperties.h"
  "simplemd/molecule-mappings/InitialPositionAndForceUpdate.h"
  "simplemd/molecule-mappings/SetMeanVelocityMapping.h"
  "simplemd/molecule-mappings/VTKMoleculeWriter.h"
  "simplemd/molecule-mappings/VelocityStoermerVerletMapping.h"
  "simplemd/molecule-mappings/VelocityStoermerVerletMapping.cpp"
  "simplemd/molecule-mappings/UpdateLinkedCellListsMapping.cpp"
  "simplemd/molecule-mappings/ComputeMeanVelocityMapping.h"
  "simplemd/molecule-mappings/UpdateLinkedCellListsMapping.h"
  "simplemd/molecule-mappings/WriteCheckPointMapping.cpp"
  "simplemd/molecule-mappings/WriteCheckPointMapping.h"
  "simplemd/molecule-mappings/VTKMoleculeWriter.cpp"
  "simplemd/BoundaryTreatment.h"
  "simplemd/LinkedCell.h"
  "simplemd/main_multi.cpp"
  "simplemd/configurations/CheckpointConfiguration.h"
  "simplemd/configurations/MoleculeConfiguration.h"
  "simplemd/configurations/MolecularDynamicsConfiguration.cpp"
  "simplemd/configurations/MolecularDynamicsConfiguration.h"
  "simplemd/configurations/DomainConfiguration.h"
  "simplemd/configurations/MoleculeConfiguration.cpp"
  "simplemd/configurations/SimulationConfiguration.h"
  "simplemd/configurations/DomainConfiguration.cpp"
  "simplemd/configurations/VTKConfiguration.h"
  "simplemd/configurations/ProfilePlotterConfiguration.cpp"
  "simplemd/configurations/RDFConfiguration.cpp"
  "simplemd/configurations/ProfilePlotterConfiguration.h"
  "simplemd/configurations/MPIConfiguration.cpp"
  "simplemd/configurations/ExternalForceConfiguration.h"
  "simplemd/configurations/SimulationConfiguration.cpp"
  "simplemd/configurations/MPIConfiguration.h"
  "simplemd/configurations/ExternalForceConfiguration.cpp"
  "simplemd/configurations/RDFConfiguration.h"
  "simplemd/configurations/VTKConfiguration.cpp"
  "simplemd/main.cpp"
  "simplemd/Molecule.h"
  "simplemd/ProfilePlotter.h"
  "simplemd/MolecularDynamicsUserInput.h"
)
foreach(clangformat_source ${CLANGFORMAT_SOURCES})
  get_filename_component(clangformat_source ${clangformat_source} ABSOLUTE)
  list(APPEND CLANGFORMAT_SOURCES_ABS ${clangformat_source})
endforeach()
add_custom_target(
  clangformat
  COMMAND /usr/bin/env clang-format --verbose
  -i
  ${CLANGFORMAT_SOURCES_ABS}
)

#
# compiletest target
#
set(COMPILETEST_COMPILERS "g++ clang++" CACHE STRING "Compiler list to be used for compiletest")
get_filename_component(COMPILETEST_SCRIPT "tarch/tests/CompileTest.sh" ABSOLUTE)
add_custom_target(
  compiletest
  COMMAND ${COMPILETEST_SCRIPT}
  ${COMPILETEST_COMPILERS}
)
