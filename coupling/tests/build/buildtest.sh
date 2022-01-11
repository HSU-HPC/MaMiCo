#!/bin/bash

### local settings like path variables
SETTINGS=../../../personal_settings

if test -f "$SETTINGS"; then
	source ../../../personal_settings
else
	echo "ERROR! No personal settings file found at $SETTINGS ."
	exit -1
fi

BUILD_PATH=${MAMICO_PATH}/coupling/tests/build
SIMPLEMD_PARALLEL_PATH=${MAMICO_PATH}/build/libsimplemd/release/dim3/parallel_yes/gcc/gprof_no/
SIMPLEMD_SEQUENTIAL_PATH=${MAMICO_PATH}/build/libsimplemd/release/dim3/parallel_no/gcc/gprof_no/
LIBSIMPLEMD=simplemd

rm ${BUILD_PATH}/test;
rm ${BUILD_PATH}/*.o;
parallel=$1;
FLAGS="-DMDDim3 -DMDCoupledDebug -std=c++17 -pedantic -Werror -Wall -DTarchDebug"


# build SimpleMD library for testing purposes
cd ${MAMICO_PATH};
if [ "${parallel}" == "parallel" ]
then
    scons target=libsimplemd dim=3 build=release parallel=yes -j4
else
    scons target=libsimplemd dim=3 build=release parallel=no -j4
fi

cd ${BUILD_PATH}
if [ "${parallel}" == "parallel" ]
then
    mpicxx ${MAMICO_PATH}/coupling/configurations/ParticleInsertionConfiguration.cpp ${FLAGS} -DMDCoupledParallel -I${MPI_INCLUDE_PATH} -I${MAMICO_PATH} -c -o ${BUILD_PATH}/ParticleInsertionConfiguration.o
    mpicxx ${MAMICO_PATH}/coupling/solvers/CoupledMolecularDynamicsSimulation.cpp ${FLAGS} -DMDCoupledParallel -I${MPI_INCLUDE_PATH} -I${MAMICO_PATH} -c -o ${BUILD_PATH}/CoupledMolecularDynamicsSimulation.o
    mpicxx ${MAMICO_PATH}/coupling/tests/main.cpp ${FLAGS} -DMDCoupledParallel -I${MPI_INCLUDE_PATH} -I${MAMICO_PATH} -c -o ${BUILD_PATH}/main.o
    mpicxx ${BUILD_PATH}/main.o ${BUILD_PATH}/ParticleInsertionConfiguration.o ${BUILD_PATH}/CoupledMolecularDynamicsSimulation.o -L${SIMPLEMD_PARALLEL_PATH} -l${LIBSIMPLEMD} -o ${BUILD_PATH}/test
else
    g++ ${MAMICO_PATH}/coupling/configurations/ParticleInsertionConfiguration.cpp ${FLAGS} -I${MAMICO_PATH} -c -o ${BUILD_PATH}/ParticleInsertionConfiguration.o
    g++ ${MAMICO_PATH}/coupling/solvers/CoupledMolecularDynamicsSimulation.cpp ${FLAGS} -I${MAMICO_PATH} -c -o ${BUILD_PATH}/CoupledMolecularDynamicsSimulation.o
    g++ ${MAMICO_PATH}/coupling/tests/main.cpp ${FLAGS} -I${MAMICO_PATH} -c -o ${BUILD_PATH}/main.o
    g++ ${BUILD_PATH}/main.o ${BUILD_PATH}/ParticleInsertionConfiguration.o ${BUILD_PATH}/CoupledMolecularDynamicsSimulation.o -L${SIMPLEMD_SEQUENTIAL_PATH} -l${LIBSIMPLEMD} -o ${BUILD_PATH}/test
fi


