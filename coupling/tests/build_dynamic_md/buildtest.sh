#!/bin/bash


### local settings like path variables
SETTINGS=../../../personal_settings

if test -f "$SETTINGS"; then
	source ${SETTINGS}
else
	echo "ERROR! No personal settings file found at $SETTINGS ."
	exit -1
fi

### build directory for library of SIMPLE_MD (currently specified for gnu compiler (intel variant: .../icc/..)
SIMPLEMD_PARALLEL_PATH=${MAMICO_PATH}/build/libsimplemd/release/dim3/parallel_yes/gcc/gprof_no/
SIMPLEMD_SEQUENTIAL_PATH=${MAMICO_PATH}/build/libsimplemd/release/dim3/parallel_no/gcc/gprof_no/
### name of lib for SIMPLE_MD
LIBSIMPLEMD=simplemd

### FROM HERE ON, ALL VARIABLES ARE AUTOMATICALLY DETERMINED ###

### build path for DynamicMDTest
BUILD_PATH=${MAMICO_PATH}/coupling/tests/build_dynamic_md;

parallel=$1;
mdSim=$2;

if [ "${parallel}" == "parallel" ] || [ "${parallel}" == "sequential" ]
then
    echo "Build mode: ${parallel}"
else
    echo "ERROR! ./buildtest parallel/sequential SIMPLE_MD/LAMMPS_MD/LAMMPS_DPD"
    exit -1
fi

if [ "${mdSim}" == "SIMPLE_MD" ] || [ "${mdSim}" == "LAMMPS_MD" ] || [ "${mdSim}" == "LAMMPS_DPD" ]
then
    echo "MD Simulation: ${mdSim}"
else
    echo "ERROR! ./test parallel/sequential SIMPLE_MD/LAMMPS_MD/LAMMPS_DPD"
    exit -1
fi



rm ${BUILD_PATH}/test;
rm ${BUILD_PATH}/*.o;

compiler=""
libraries=""
objects=""
includes="-I${MAMICO_PATH}"

### specify flags, includes, libraries,compiler for parallel or sequential build
if [ "${parallel}" == "parallel" ]
then
    # note: we need to set MDDim3 for ALL Simulations since we use the configuration classes from SimpleMD
    FLAGS="-D${mdSim} -DMDDim3 -std=c++1z -pedantic -Werror -Wno-unknown-pragmas -Wno-int-in-bool-context -Wall -DMDCoupledParallel -DTarchParallel -DMPICH_IGNORE_CXX_SEEK -O3 -g3"
    # -DMDCoupledDebug"
    includes="${includes} -I${MPI_INCLUDE_PATH} -I${LIB_EIGEN_PATH}"
    libraries="-L${MPI_LIB_PATH} -l${LIB_MPI}"
    compiler="mpicxx"
else
    FLAGS="-D${mdSim} -DMDDim3 -std=c++1z -pedantic -Wall -Wno-unknown-pragmas -O3"
    # -Werror
    includes="${includes} -I${LIB_EIGEN_PATH}"
    compiler="g++"
fi
### 

### builds, objects, libraries for coupling -> we require several parts from simplemd
cd ${MAMICO_PATH} || exit



# specific built for SIMPLE_MD
if [ "${mdSim}" == "SIMPLE_MD" ]
then
  if [ "${parallel}" == "parallel" ]
  then
    scons target=libsimplemd dim=3 build=release parallel=yes compiler=gcc machine=hww-cluster -j4
    libraries="${libraries} -L${SIMPLEMD_PARALLEL_PATH} -l${LIBSIMPLEMD}"
    FLAGS="${FLAGS} -DMDParallel"
  else
    scons target=libsimplemd dim=3 build=release parallel=no compiler=gcc machine=hww-cluster -j4
    libraries="${libraries} -L${SIMPLEMD_SEQUENTIAL_PATH} -l${LIBSIMPLEMD}"
  fi
  cd ${BUILD_PATH} || exit 
    ${compiler} ${MAMICO_PATH}/coupling/solvers/CoupledMolecularDynamicsSimulation.cpp ${FLAGS} ${includes} -c -o ${BUILD_PATH}/CoupledMolecularDynamicsSimulation.o
    objects="${objects} ${BUILD_PATH}/CoupledMolecularDynamicsSimulation.o"
elif [ "${mdSim}" == "LAMMPS_MD" ]
then
  if [ "${parallel}" == "sequential" ]
  then
    echo "ERROR: LAMMPS only works for option 'parallel'"
    exit -1
  fi
  
  cp -a ${MAMICO_PATH}/coupling/interface/impl/LAMMPS/USER-MAMICO ${LAMMPS_PATH}/src/

  export LIB_EIGEN_PATH
  export MAMICO_PATH

  # build lammps
  cd ${LAMMPS_PATH}/src || exit
  make yes-user-mamico -j4 || exit
  make g++_openmpi -j4 || exit
  make makelib
  make -f Makefile.lib openmpi
  ln -s liblammps_g++_openmpi.a lib${LIB_LAMMPS}.a
  includes="${includes} -I${LAMMPS_PATH}/src"
  libraries="${libraries} -L${LAMMPS_PATH}/src -Wl,-Bstatic -l${LIB_LAMMPS} -Wl,-Bdynamic"
  # remove Werror flag
  FLAGS=`echo ${FLAGS} | sed 's/-Werror//'`
elif [ "${mdSim}" == "LAMMPS_DPD" ]
then
  if [ "${parallel}" == "sequential" ]
  then
    echo "ERROR: LAMMPS only works for option 'parallel'"
    exit -1
  fi
  # build lammps
  cd ${LAMMPS_PATH} || exit
  make yes-user-mamico -j4
  make openmpi -j4
  make makelib -j4
  make -f Makefile.lib openmpi
  ln -s liblammps_g++_openmpi.a lib${LIB_LAMMPS}.a
  includes="${includes} -I${LAMMPS_PATH}/src"
  libraries="${libraries} -L${LAMMPS_PATH}/src -Wl,-Bstatic -l${LIB_LAMMPS} -Wl,-Ddynamic"
  # remove Werror flag
  FLAGS=`echo ${FLAGS} | sed 's/-Werror//'`
fi

### builds, linking, objects for coupled simulation with MaMiCo
cd ${BUILD_PATH} || exit
${compiler} ${MAMICO_PATH}/tarch/tinyxml2/tinyxml2.cpp ${includes} ${libraries} -c -o ${BUILD_PATH}/tinyxml2.o
${compiler} ${MAMICO_PATH}/tarch/utils/RandomNumberService.cpp ${includes} ${libraries} -c -o ${BUILD_PATH}/RandomNumberService.o
${compiler} ${MAMICO_PATH}/coupling/configurations/ParticleInsertionConfiguration.cpp ${FLAGS} ${includes} -c -o ${BUILD_PATH}/ParticleInsertionConfiguration.o
${compiler} ${MAMICO_PATH}/simplemd/configurations/VTKConfiguration.cpp ${FLAGS} ${includes} -c -o ${BUILD_PATH}/VTKConfiguration.o
${compiler} ${MAMICO_PATH}/simplemd/configurations/MolecularDynamicsConfiguration.cpp ${FLAGS} ${includes} -c -o ${BUILD_PATH}/MolecularDynamicsConfiguration.o
${compiler} ${MAMICO_PATH}/simplemd/configurations/DomainConfiguration.cpp ${FLAGS} ${includes} -c -o ${BUILD_PATH}/DomainConfiguration.o
${compiler} ${MAMICO_PATH}/simplemd/configurations/ExternalForceConfiguration.cpp ${FLAGS} ${includes} -c -o ${BUILD_PATH}/ExternalForceConfiguration.o
${compiler} ${MAMICO_PATH}/simplemd/configurations/MoleculeConfiguration.cpp ${FLAGS} ${includes} -c -o ${BUILD_PATH}/MoleculeConfiguration.o
${compiler} ${MAMICO_PATH}/simplemd/configurations/MPIConfiguration.cpp ${FLAGS} -DMDParallel ${includes} -c -o ${BUILD_PATH}/MPIConfiguration.o
${compiler} ${MAMICO_PATH}/simplemd/configurations/ProfilePlotterConfiguration.cpp ${FLAGS} ${includes} -c -o ${BUILD_PATH}/ProfilePlotterConfiguration.o
${compiler} ${MAMICO_PATH}/simplemd/configurations/RDFConfiguration.cpp ${FLAGS} ${includes} -c -o ${BUILD_PATH}/RDFConfiguration.o
${compiler} ${MAMICO_PATH}/simplemd/configurations/SimulationConfiguration.cpp ${FLAGS} ${includes} -c -o ${BUILD_PATH}/SimulationConfiguration.o
${compiler} ${MAMICO_PATH}/coupling/tests/main_dynamic_md.cpp ${FLAGS} ${includes} -c -o ${BUILD_PATH}/main_dynamic_md.o
objects="${objects} ${BUILD_PATH}/ParticleInsertionConfiguration.o ${BUILD_PATH}/main_dynamic_md.o ${BUILD_PATH}/tinyxml2.o ${BUILD_PATH}/RandomNumberService.o ${BUILD_PATH}/VTKConfiguration.o ${BUILD_PATH}/MolecularDynamicsConfiguration.o ${BUILD_PATH}/DomainConfiguration.o ${BUILD_PATH}/ExternalForceConfiguration.o ${BUILD_PATH}/MoleculeConfiguration.o ${BUILD_PATH}/MPIConfiguration.o ${BUILD_PATH}/ProfilePlotterConfiguration.o ${BUILD_PATH}/RDFConfiguration.o ${BUILD_PATH}/SimulationConfiguration.o"

echo "${compiler} ${objects} ${libraries} -o ${BUILD_PATH}/dynamicMDTest" 
${compiler} ${objects} ${libraries} -o ${BUILD_PATH}/dynamicMDTest_${parallel}

