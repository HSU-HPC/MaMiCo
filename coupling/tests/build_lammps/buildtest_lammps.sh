#!/bin/bash

### local settings like path variables
SETTINGS=../../../personal_settings

if test -f "$SETTINGS"; then
	source ../../../personal_settings
else
	echo "ERROR! No personal settings file found at $SETTINGS ."
	exit -1
fi

### path variables for this script
DEFINES="-DMDCoupledDebug -DMDCoupledParallel -std=c++11 -DTarchDebug";
BUILD_PATH=${MAMICO_PATH}/coupling/tests/build_lammps;
COMPILER=mpic++;

rm ${BUILD_PATH}/test_lammps;
rm ${BUILD_PATH}/*.o;

${COMPILER} ${MAMICO_PATH}/tarch/tinyxml2/tinyxml2.cpp ${DEFINES} -c -o ${BUILD_PATH}/tinyxml2.o
${COMPILER} ${MAMICO_PATH}/coupling/configurations/ParticleInsertionConfiguration.cpp ${DEFINES} -I${MPI_PATH_INCLUDE} -I${MAMICO_PATH} -I${LAMMPS_PATH} -c -o ${BUILD_PATH}/ParticleInsertionConfiguration.o
${COMPILER} ${MAMICO_PATH}/tarch/utils/RandomNumberService.cpp ${DEFINES} -I${MPI_PATH_INCLUDE} -I${MAMICO_PATH} -I${LAMMPS_PATH} -c -o ${BUILD_PATH}/RandomNumberService.o
${COMPILER} ${MAMICO_PATH}/coupling/tests/main_lammps.cpp ${DEFINES} -I${MPI_PATH_INCLUDE} -I${MAMICO_PATH} -I${LAMMPS_PATH} -c -o ${BUILD_PATH}/main_lammps.o
${COMPILER} -o ${BUILD_PATH}/test_lammps ${BUILD_PATH}/tinyxml2.o ${BUILD_PATH}/RandomNumberService.o ${BUILD_PATH}/ParticleInsertionConfiguration.o ${BUILD_PATH}/main_lammps.o -L${MPI_PATH_LIB} -l${LIB_MPI} -L${LAMMPS_PATH} -l${LIB_LAMMPS}

