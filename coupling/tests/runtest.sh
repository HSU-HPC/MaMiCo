#!/bin/bash

# run sequential test
cd build
rm output_test_*
echo "Run tests and write output into output_test_sequential.out"
./buildtest.sh
./test &> output_test_sequential.out
# run parallel test
./buildtest.sh parallel
for i in 1 4 8 9 27; do
  echo "Run parallel test with $i procs and write output into output_test_${i}proc.out";
  mpirun -np ${i} ./test &> output_test_${i}proc.out;
done
