#!/bin/bash

### local settings like path variables
SETTINGS=../../../personal_settings

if test -f "$SETTINGS"; then
	source ../../../personal_settings
else
	echo "ERROR! No personal settings file found at $SETTINGS ."
	exit -1
fi

BUILD_PATH=${MAMICO_PATH}/coupling/tests/build_lbcouettesolver

rm ${BUILD_PATH}/test;
rm ${BUILD_PATH}/*.o;
FLAGS="-std=c++11 -pedantic -Werror -Wall -O3 -fopenmp -lgomp -lmpi -DMDCoupledParallel -fno-strict-overflow"
#FLAGS="-std=c++11 -pedantic -Wall -O3 -DMDCoupledParallel"


mpicxx ${FLAGS} ${MAMICO_PATH}/coupling/tests/main_lbcouettesolver.cpp -I${MAMICO_PATH} -I${MPI_INCLUDE_PATH} -L${MPI_LIB_PATH} -lpthread -lmpi -o ${BUILD_PATH}/test


