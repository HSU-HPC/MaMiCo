#!/bin/bash

FILENAME=infotest_espresso.dat

# build test executable
./build_espresso/buildtest_espresso.sh 

#mpirun -np 1 valgrind --leak-check=full ./test_espresso &> infotest_espresso.dat;
cd build_espresso
rm ${FILENAME}

echo "Run tests; write output to ${FILENAME} ..."
for (( i=0; i<6; i++)); do
  echo "Run test $i" >> ${FILENAME}
  mpirun -np 1 ./test_espresso $i >> ${FILENAME}
done
echo "Run test 6" >> ${FILENAME}
mpirun -np 4 ./test_espresso 6 >> ${FILENAME}
cd ..
