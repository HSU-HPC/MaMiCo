#!/bin/bash

### path variables; adapt these to your installation
MARDYN_PATH=/home/neumanph/workspace/mardyn_mamico/src;
MARDYN_PATH_EXTERNAL=${MARDYN_PATH}/External;
MPI_PATH_INCLUDE=/usr/lib/openmpi/include;
MPI_PATH_LIB=/usr/lib/openmpi/lib;
MPI_LIB=mpi;
MAMICO_PATH=/home/neumanph/workspace/mamico/mamico_cpc_v1.1;


# further variables; leave those unchanged (except for the compiler)
BUILD_PATH=${MAMICO_PATH}/coupling/tests/build_mardyn;
COMPILER=mpic++;
MARDYN_LIB=mardyn;
DEFINES="-DMDDim3 -DMDCoupledError -DMDCoupledParallel -DENABLE_MPI -DMPICH_IGNORE_CXX_SEEK -std=c++11"

rm ${BUILD_PATH}/test_mardyn;
rm ${BUILD_PATH}/*.o;

echo "Build..."
# tarch, and coupling classes
${COMPILER} ${MAMICO_PATH}/tarch/utils/RandomNumberService.cpp ${DEFINES} -I${MPI_PATH_INCLUDE} -I${MAMICO_PATH} -c -o ${BUILD_PATH}/RandomNumberService.o
${COMPILER} ${MAMICO_PATH}/coupling/configurations/ParticleInsertionConfiguration.cpp ${DEFINES} -I${MPI_PATH_INCLUDE} -I${MAMICO_PATH} -c -o ${BUILD_PATH}/ParticleInsertionConfiguration.o
${COMPILER} ${MAMICO_PATH}/tarch/tinyxml2/tinyxml2.cpp ${DEFINES} -I{MPI_PATH_INCLUDE} -I${MAMICO_PATH} -c -o ${BUILD_PATH}/tinyxml2.o
${COMPILER} ${MAMICO_PATH}/coupling/solvers/DummySolverInterfaceService.cpp ${DEFINES} -I${MPI_PATH_INCLUDE} -I${TARCH_PATH} -I${MAMICO_PATH} -c -o ${BUILD_PATH}/DummySolverInterfaceService.o
#executable
${COMPILER} -g ${MAMICO_PATH}/coupling/tests/main_mardyn.cpp ${DEFINES} -I${MPI_PATH_INCLUDE} -I${MAMICO_PATH} -I${MARDYN_PATH} -I${MARDYN_PATH_EXTERNAL} -c -o ${BUILD_PATH}/main_mardyn.o
${COMPILER} -o ${BUILD_PATH}/test_mardyn ${BUILD_PATH}/main_mardyn.o  ${BUILD_PATH}/ParticleInsertionConfiguration.o ${BUILD_PATH}/tinyxml2.o ${BUILD_PATH}/RandomNumberService.o ${BUILD_PATH}/DummySolverInterfaceService.o -L${MARDYN_PATH} -l${MARDYN_LIB}
