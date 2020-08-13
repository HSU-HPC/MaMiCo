#!/bin/bash

### include and library paths for MPI and name of MPI library
MPI_INCLUDE_PATH=/usr/lib/openmpi/include
MPI_LIB_PATH=/usr/lib/openmpi/lib
LIB_MPI=mpi

LIB_EIGEN_PATH=/usr/include/eigen3

### home directory of MAMICO
MAMICO_PATH=/home/$USER/Documents/MaMiCo/mamico-dev

### build directory for library of SIMPLE_MD (currently specified for gnu compiler (intel variant: .../icc/..)
SIMPLEMD_PARALLEL_PATH=${MAMICO_PATH}/build/libsimplemd/release/dim3/parallel_yes/gcc/gprof_no/
SIMPLEMD_SEQUENTIAL_PATH=${MAMICO_PATH}/build/libsimplemd/release/dim3/parallel_no/gcc/gprof_no/
### name of lib for SIMPLE_MD
LIBSIMPLEMD=simplemd

### FROM HERE ON, ALL VARIABLES ARE AUTOMATICALLY DETERMINED ###

### build path for CouetteTest
BUILD_PATH=${MAMICO_PATH}/coupling/tests/build_couette;

parallel=$1;

if [ "${parallel}" == "parallel" ] || [ "${parallel}" == "sequential" ]
then
    echo "Build mode: ${parallel}"
else
    echo "ERROR! ./buildtest parallel/sequential"
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
    FLAGS="-DSIMPLE_MD -DMDDim3 -std=c++1z -pedantic -Werror -Wno-unknown-pragmas -Wall -DMDCoupledParallel -DTarchParallel -DMPICH_IGNORE_CXX_SEEK -O3"
    # -DMDCoupledDebug"
    includes="${includes} -I${MPI_INCLUDE_PATH} -I${LIB_EIGEN_PATH}"
    libraries="-L${MPI_LIB_PATH} -l${LIB_MPI}"
    compiler="mpicxx"
else
    FLAGS="-DSIMPLE_MD -DMDDim3 -std=c++1z -pedantic -Wall -Wno-unknown-pragmas -O3"
    # -Werror
    includes="${includes} -I${LIB_EIGEN_PATH}"
    compiler="g++"
fi
### 

### builds, objects, libraries for coupling -> we require several parts from simplemd
cd ${MAMICO_PATH}
if [ "${parallel}" == "parallel" ]
then
        scons target=libsimplemd dim=3 build=release parallel=yes -j4
        libraries="${libraries} -L${SIMPLEMD_PARALLEL_PATH} -l${LIBSIMPLEMD}"
        FLAGS="${FLAGS} -DMDParallel"
else
        scons target=libsimplemd dim=3 build=release parallel=no -j4
        libraries="${libraries} -L${SIMPLEMD_SEQUENTIAL_PATH} -l${LIBSIMPLEMD}"
fi

cd ${BUILD_PATH}
${compiler} ${MAMICO_PATH}/coupling/solvers/CoupledMolecularDynamicsSimulation.cpp ${FLAGS} ${includes} -c -o ${BUILD_PATH}/CoupledMolecularDynamicsSimulation.o
objects="${objects} ${BUILD_PATH}/CoupledMolecularDynamicsSimulation.o"

### builds, linking, objects for coupled simulation with MaMiCo
cd ${BUILD_PATH}
${compiler} ${MAMICO_PATH}/coupling/configurations/ParticleInsertionConfiguration.cpp ${FLAGS} ${includes} -c -o ${BUILD_PATH}/ParticleInsertionConfiguration.o
${compiler} ${MAMICO_PATH}/coupling/tests/main_couette.cpp ${FLAGS} ${includes} -c -o ${BUILD_PATH}/main_couette.o
objects="${objects} ${BUILD_PATH}/ParticleInsertionConfiguration.o ${BUILD_PATH}/main_couette.o"

${compiler} ${objects} ${libraries} -o ${BUILD_PATH}/test

