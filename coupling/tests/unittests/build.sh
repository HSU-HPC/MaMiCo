#!/bin/bash

SETTINGS=../../../personal_settings

source ${SETTINGS}

#Assumes simplemd to be compiled before
#To do so, go to MAMICO_PATH and execute:
#	scons target=libsimplemd dim=3 build=release parallel=no -j4
SIMPLEMD_PATH=${MAMICO_PATH}/build/libsimplemd/release/dim3/parallel_no/gcc/gprof_no/

flags="-DSIMPLE_MD -DMDDim3 --std=c++20 -O3"
includes="-I${MPI_INCLUDE_PATH} -I${MAMICO_PATH}"
libraries="-L${MPI_LIB_PATH} -l${LIB_MPI} -L${SIMPLEMD_PATH} -lsimplemd"

mpic++ ${flags} ${includes} ${libraries} main.cpp
