#!/bin/bash

### path variables for this script
DEFINES="-DMDCoupledDebug -DMDCoupledParallel -std=c++11 -DTarchDebug";
MPI_PATH_INCLUDE=/usr/lib/openmpi/include;
MPI_PATH_LIB=/usr/lib/openmpi/lib;
LIB_MPI=mpi;
MAMICO_PATH=/home/neumanph/workspace/mamico/mamico_cpc_v1.1;
BUILD_PATH=${MAMICO_PATH}/coupling/tests/build_lammps;
LAMMPS_PATH=/home/neumanph/programs/lammps/src;
LIB_LAMMPS=lammps_openmpi_mamico;
COMPILER=mpic++;

rm ${BUILD_PATH}/test_lammps;
rm ${BUILD_PATH}/*.o;

${COMPILER} ${MAMICO_PATH}/tarch/tinyxml2/tinyxml2.cpp ${DEFINES} -c -o ${BUILD_PATH}/tinyxml2.o
${COMPILER} ${MAMICO_PATH}/coupling/configurations/ParticleInsertionConfiguration.cpp ${DEFINES} -I${MPI_PATH_INCLUDE} -I${MAMICO_PATH} -I${LAMMPS_PATH} -c -o ${BUILD_PATH}/ParticleInsertionConfiguration.o
${COMPILER} ${MAMICO_PATH}/tarch/utils/RandomNumberService.cpp ${DEFINES} -I${MPI_PATH_INCLUDE} -I${MAMICO_PATH} -I${LAMMPS_PATH} -c -o ${BUILD_PATH}/RandomNumberService.o
${COMPILER} ${MAMICO_PATH}/coupling/tests/main_lammps.cpp ${DEFINES} -I${MPI_PATH_INCLUDE} -I${MAMICO_PATH} -I${LAMMPS_PATH} -c -o ${BUILD_PATH}/main_lammps.o
${COMPILER} -o ${BUILD_PATH}/test_lammps ${BUILD_PATH}/tinyxml2.o ${BUILD_PATH}/RandomNumberService.o ${BUILD_PATH}/ParticleInsertionConfiguration.o ${BUILD_PATH}/main_lammps.o -L${MPI_PATH_LIB} -l${LIB_MPI} -L${LAMMPS_PATH} -l${LIB_LAMMPS}

