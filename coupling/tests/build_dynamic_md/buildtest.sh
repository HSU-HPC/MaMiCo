#!/bin/bash

### include and library paths for MPI and name of MPI library
MPI_INCLUDE_PATH=$(pkg-config --cflags-only-I ompi)
MPI_LIB_PATH=$(pkg-config --libs-only-L ompi)
LIB_MPI=$(pkg-config --libs-only-l ompi)

LIB_EIGEN_PATH=$(pkg-config --cflags-only-I eigen3)

### home directory of MAMICO
MAMICO_PATH=/home/niklas/Dokumente/Git/mamico-dev

### build directory for library of SIMPLE_MD (currently specified for gnu compiler (intel variant: .../icc/..)
SIMPLEMD_PARALLEL_PATH=${MAMICO_PATH}/build/libsimplemd/release/dim3/parallel_yes/gcc/gprof_no/
### name of lib for SIMPLE_MD
LIBSIMPLEMD=simplemd

### FROM HERE ON, ALL VARIABLES ARE AUTOMATICALLY DETERMINED ###

### build path for DynamicMDTest
BUILD_PATH=${MAMICO_PATH}/coupling/tests/build_dynamic_md;

dynamic=$1;

if [ "${dynamic}" == "sudden" ] || [ "${dynamic}" == "successive" ]
then
    echo "Build mode: ${dynamic}"
else
    echo "ERROR! ./test sudden/successive"
    exit 1
fi



rm ${BUILD_PATH}/test;
rm ${BUILD_PATH}/*.o;

compiler=""
libaries=""
objects=""
includes="-I${MAMICO_PATH}"

### specify flags, includes, libraries,compiler for parallel build
# note: we need to set MDDim3 for ALL Simulations since we use the configuration classes from SimpleMD
FLAGS="-DSIMPLE_MD -DMDDim3 -std=c++1z -pedantic -Werror -Wno-unknown-pragmas -Wall -DMDCoupledParallel -DTarchParallel -DMPICH_IGNORE_CXX_SEEK -O0 -g3"
# -DMDCoupledDebug"
includes="${includes} ${MPI_INCLUDE_PATH} ${LIB_EIGEN_PATH}"
libraries="${MPI_LIB_PATH} ${LIB_MPI}"
compiler="mpicxx"

### 

### builds, objects, libraries for coupling -> we require several parts from simplemd
cd ${MAMICO_PATH} || exit
scons target=libsimplemd dim=3 build=release parallel=yes -j4
libraries="${libraries} -L${SIMPLEMD_PARALLEL_PATH} -l${LIBSIMPLEMD}"
FLAGS="${FLAGS} -DMDParallel"


if [ "${dynamic}" == "sudden" ]; then
    FLAGS="${FLAGS} -DCOUPLING_DYNAMIC_MD_SUDDEN"
else
    FLAGS="${FLAGS} -DCOUPLING_DYNAMIC_MD_SUCCESSIVE"
fi

cd ${BUILD_PATH} || exit
${compiler} ${MAMICO_PATH}/coupling/solvers/CoupledMolecularDynamicsSimulation.cpp ${FLAGS} ${includes} -c -o ${BUILD_PATH}/CoupledMolecularDynamicsSimulation.o
objects="${objects} ${BUILD_PATH}/CoupledMolecularDynamicsSimulation.o"

### builds, linking, objects for coupled simulation with MaMiCo
cd ${BUILD_PATH} || exit
${compiler} ${MAMICO_PATH}/coupling/configurations/ParticleInsertionConfiguration.cpp ${FLAGS} ${includes} -c -o ${BUILD_PATH}/ParticleInsertionConfiguration.o
${compiler} ${MAMICO_PATH}/coupling/tests/main_dynamic_md.cpp ${FLAGS} ${includes} -c -o ${BUILD_PATH}/main_dynamic_md.o
objects="${objects} ${BUILD_PATH}/ParticleInsertionConfiguration.o ${BUILD_PATH}/main_dynamic_md.o"

echo "${compiler} ${objects} ${libraries} -o ${BUILD_PATH}/test" 
${compiler} ${objects} ${libraries} -o ${BUILD_PATH}/test

