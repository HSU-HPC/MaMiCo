#!/bin/bash

# check input arguments
if [ "$1" == "sequential" ] || [ "$1" == "parallel" ]
then
  echo "Simulation mode: $1"
else
  echo "ERROR: argument 1: sequential/parallel"
  exit -1
fi

if [ "$2" == "SIMPLE_MD" ] || [ "$2" == "LAMMPS_MD" ] || [ "$2" == "LAMMPS_DPD" ]
then
  echo "Simulation type: $2"
else
  echo "ERROR: argument 2: SIMPLE_MD/LAMMPS_MD/LAMMPS_DPD"
  exit -1
fi


# run NieTest: change into build folder and remove previous output
cd build_nie;
rm *csv *vtk *log.lammps* *out

# copy correct files to xml configs and header describing MD/DPD setup
if [ "$2" == "LAMMPS_DPD" ]
then
  cp NieTestDPD.h ../NieTest.h
  cp test_nie_mamico_dpd.xml test_nie.xml
  cp test_nie_dpd.xml test_nie_simplemd.xml
else
  cp NieTestMD.h ../NieTest.h
  cp test_nie_mamico_md.xml test_nie.xml
  cp test_nie_md.xml test_nie_simplemd.xml
fi

# build executable and break on failure
./buildtest.sh $1 $2 || { echo "ERROR!" ; exit -1; }

echo "Buildtest successful, execute simulation"

# execute coupled simulation using either 2 MPI ranks or a single process
if [ "$1" == "parallel" ]
then
  mpirun -np 2 ./test &> $2.out
else
  ./test &> $2.out
fi
