#!/bin/bash

### path variables for this script
MAMICO_PATH=/home/neumanph/workspace/mamico/mamico_cpc_v1.1;
BUILD_PATH=${MAMICO_PATH}/coupling/tests/build_lbcouettesolver

rm ${BUILD_PATH}/test;
rm ${BUILD_PATH}/*.o;
FLAGS="-std=c++11 -pedantic -Werror -Wall -O3 -fopenmp -lgomp -lmpi -DMDCoupledParallel -fno-strict-overflow"
#FLAGS="-std=c++11 -pedantic -Wall -O3 -DMDCoupledParallel"


mpicxx ${FLAGS} ${MAMICO_PATH}/coupling/tests/main_lbcouettesolver.cpp -I${MAMICO_PATH} -I/usr/lib/openmpi/include -L/usr/lib/openmpi/lib -lpthread -lmpi -o ${BUILD_PATH}/test


