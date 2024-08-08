#!/bin/bash

FILENAME=infotest_mardyn.dat

# build test executable
./build_mardyn/buildtest_mardyn.sh 

cd build_mardyn
rm ${FILENAME}

echo "Run tests; write output to ${FILENAME} ..."
mpirun -np 1 ./test_mardyn  &> ${FILENAME}

cd ..
