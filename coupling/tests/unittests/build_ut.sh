#!/bin/bash

### local settings like path variables
SETTINGS=../../../personal_settings

if test -f "$SETTINGS"; then
	source ${SETTINGS}
else
	echo "ERROR! No personal settings file found at $SETTINGS ."
	exit -1
fi

SIMPLEMD_PATH=${MAMICO_PATH}/build/libsimplemd/release/dim3/parallel_yes/gcc/gprof_no/
LIBSIMPLEMD=simplemd

BUILD_PATH=${MAMICO_PATH}/coupling/tests/unittests/build

rm ${BUILD_PATH}/main;
rm ${BUILD_PATH}/*.o;

#set flags, includes, libraries
FLAGS="-DSIMPLE_MD -DMDDim3 -std=c++20 -Wno-unknown-pragmas -Wno-int-in-bool-context -Wall -DMDParallel -DMDCoupledParallel -DTarchParallel -DMPICH_IGNORE_CXX_SEEK -O3"
includes="-I${MAMICO_PATH} -I${MPI_INCLUDE_PATH}"
libraries="-L${MPI_LIB_PATH} -l${LIB_MPI}"
compiler="mpicxx"

#build simplemd
cd ${MAMICO_PATH}
scons target=libsimplemd dim=3 build=release parallel=yes -j4
libraries="${libraries} -L${SIMPLEMD_PATH} -l${LIBSIMPLEMD}"
cd ${BUILD_PATH}

${compiler} ${MAMICO_PATH}/coupling/solvers/CoupledMolecularDynamicsSimulation.cpp ${FLAGS} ${includes} -c -o ${BUILD_PATH}/CoupledMolecularDynamicsSimulation.o
${compiler} ${MAMICO_PATH}/coupling/configurations/ParticleInsertionConfiguration.cpp ${FLAGS} ${includes} -c -o ${BUILD_PATH}/ParticleInsertionConfiguration.o
${compiler} ${MAMICO_PATH}/coupling/tests/unittests/main.cpp ${FLAGS} ${includes} -c -o ${BUILD_PATH}/main_ut.o -L/usr/lib/python3.8/site-packages/pybind11 #$(python3.8-config --ldflags --embed)

cd ..
objects="${BUILD_PATH}/ParticleInsertionConfiguration.o ${BUILD_PATH}/main_ut.o ${BUILD_PATH}/CoupledMolecularDynamicsSimulation.o"

${compiler} ${objects} ${libraries} -o main
