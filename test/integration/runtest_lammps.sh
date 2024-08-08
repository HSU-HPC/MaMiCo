#!/bin/bash

# build test executable
echo "Build test executable..."
FILE_STEM="infotest_lammps";
cd build_lammps
./buildtest_lammps.sh

# run tests for 1,4,8,16,64 procs
for i in 1 4 8 16 64; do
  echo "Run test for ${i} procs";
  FILE="${FILE_STEM}${i}proc.dat";
  mpirun -np $i ./test_lammps &> ${FILE};
  echo "Errors reported:";
  cat ${FILE} | grep ERROR;
done

