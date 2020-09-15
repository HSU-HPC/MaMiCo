#!/bin/bash

# !!! Minimum required: openmpi version  >= 3.0 !!!

### include and library paths for MPI and name of MPI library
MPI_INCLUDE_PATH=/usr/lib/x86_64-linux-gnu/openmpi/include
MPI_LIB_PATH=/usr/lib/x86_64-linux-gnu/openmpi/lib
LIB_MPI=mpi

LIB_EIGEN_PATH=/usr/local/include/Eigen

### home directory of MAMICO
MAMICO_PATH=/home/$USER/Documents/MaMiCo/mamico-dev

### build directory for library of SIMPLE_MD (currently specified for gnu compiler (intel variant: .../icc/..)
SIMPLEMD_PATH=${MAMICO_PATH}/build/libsimplemd/release/dim3/parallel_yes/gcc/gprof_no/
### name of lib for SIMPLE_MD
LIBSIMPLEMD=simplemd

echo "Compiling mamico python bindings..."

### build path for CouetteTest
BUILD_PATH=${MAMICO_PATH}/coupling/python-binding;

TARGET=${BUILD_PATH}/mamico`python3-config --extension-suffix`;

rm -f ${TARGET}; 
rm -f ${BUILD_PATH}/*.o; 

### specify flags, includes, libraries,compiler for parallel or sequential build
# note: we need to set MDDim3 for ALL Simulations since we use the configuration classes from SimpleMD
FLAGS="-DSIMPLE_MD -DMDDim3 -std=c++1z -pedantic -Werror -Wno-unknown-pragmas -Wall -DMDCoupledParallel -fPIC -DTarchParallel -DMPICH_IGNORE_CXX_SEEK -O3"
# -DMDCoupledDebug"
includes="-I${MAMICO_PATH} -I${MPI_INCLUDE_PATH} -I${LIB_EIGEN_PATH} `python3 -m pybind11 --includes`"
libraries="-L${MPI_LIB_PATH} -l${LIB_MPI}"
compiler="mpicxx"

### builds, objects, libraries for coupling -> we require several parts from simplemd
cd ${MAMICO_PATH} || exit $?
/usr/bin/env python $(which scons) target=libsimplemd dim=3 build=release parallel=yes -j4 || exit $?
libraries="${libraries} -L${SIMPLEMD_PATH} -l${LIBSIMPLEMD}"
FLAGS="${FLAGS} -DMDParallel"

cd ${BUILD_PATH} || exit $?
${compiler} ${MAMICO_PATH}/coupling/solvers/CoupledMolecularDynamicsSimulation.cpp ${FLAGS} ${includes} -c -o ${BUILD_PATH}/CoupledMolecularDynamicsSimulation.o || exit $?
objects="${BUILD_PATH}/CoupledMolecularDynamicsSimulation.o"

### builds, linking, objects for coupled simulation with MaMiCo
${compiler} ${MAMICO_PATH}/coupling/configurations/ParticleInsertionConfiguration.cpp ${FLAGS} ${includes} -c -o ${BUILD_PATH}/ParticleInsertionConfiguration.o || exit $?
${compiler} ${MAMICO_PATH}/coupling/python-binding/mamico.cpp ${FLAGS} ${includes} -c -o ${BUILD_PATH}/mamico.o || exit $?
objects="${objects} ${BUILD_PATH}/ParticleInsertionConfiguration.o ${BUILD_PATH}/mamico.o"

${compiler} ${objects} ${libraries} -shared -o ${TARGET} || exit $?

echo "Successfully finished compiling mamico python bindings to " ${TARGET}
