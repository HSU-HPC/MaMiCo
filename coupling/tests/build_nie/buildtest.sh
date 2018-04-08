#!/bin/bash

### include and library paths for MPI and name of MPI library
MPI_INCLUDE_PATH=/usr/lib/openmpi/include
MPI_LIB_PATH=/usr/lib/openmpi/lib
LIB_MPI=mpi

### home directory of MAMICO
MAMICO_PATH=/home/piet/mamico_v1.1/

### build directory for library of SIMPLE_MD (currently specified for gnu compiler (intel variant: .../icc/..)
SIMPLEMD_PARALLEL_PATH=${MAMICO_PATH}/build/libsimplemd/release/dim3/parallel_yes/gcc/gprof_no/
SIMPLEMD_SEQUENTIAL_PATH=${MAMICO_PATH}/build/libsimplemd/release/dim3/parallel_no/gcc/gprof_no/
### name of lib for SIMPLE_MD
LIBSIMPLEMD=simplemd

### LAMMPS
LAMMPS_PATH=/home/neumanph/programs/lammps/src/
LIBLAMMPS=lammps_openmpi_mamico


### FROM HERE ON, ALL VARIABLES ARE AUTOMATICALLY DETERMINED ###

### build path for Nie-test
BUILD_PATH=${MAMICO_PATH}/coupling/tests/build_nie;

parallel=$1;
mdSim=$2;

if [ "${parallel}" == "parallel" ] || [ "${parallel}" == "sequential" ]
then
    echo "Build mode: parallel"
else
    echo "ERROR! ./test parallel/sequential SIMPLE_MD/LAMMPS_MD/LAMMPS_DPD"
    exit -1
fi
if [ "${mdSim}" == "SIMPLE_MD" ] || [ "${mdSim}" == "LAMMPS_MD" ] || [ "${mdSim}" == "LAMMPS_DPD" ]
then
    echo "MD Simulation: ${mdSim}"
else
    echo "ERROR! ./test parallel/sequential SIMPLE_MD/LAMMPS_MD/LAMMPS_DPD"
    exit -1
fi

rm ${BUILD_PATH}/test;
rm ${BUILD_PATH}/*.o;

compiler=""
libaries=""
objects=""
includes="-I${MAMICO_PATH}"

### specify flags, includes, libraries,compiler for parallel or sequential build
if [ "${parallel}" == "parallel" ]
then
    # note: we need to set MDDim3 for ALL Simulations since we use the configuration classes from SimpleMD
    FLAGS="-D${mdSim} -DMDDim3 -std=c++0x -pedantic -Werror -Wno-unknown-pragmas -Wall -DMDCoupledParallel -DTarchParallel -DMPICH_IGNORE_CXX_SEEK -O3"
    # -DMDCoupledDebug"
    includes="${includes} -I${MPI_INCLUDE_PATH}"
    libraries="-L${MPI_LIB_PATH} -l${LIB_MPI}"
    compiler="mpicxx"
else
    FLAGS="-D${mdSim} -DMDDim3 -std=c++0x -pedantic -Werror -Wall -Wno-unknown-pragmas -O3"
    compiler="g++"
fi


### 

### builds, objects, libraries for coupling -> we require several parts from simplemd
cd ${MAMICO_PATH}
if [ "${parallel}" == "parallel" ]
then
        scons target=libsimplemd dim=3 build=release parallel=yes -j4
        libraries="${libraries} -L${SIMPLEMD_PARALLEL_PATH} -l${LIBSIMPLEMD}"
        FLAGS="${FLAGS} -DMDParallel"
else
        scons target=libsimplemd dim=3 build=release parallel=no -j4
        libraries="${libraries} -L${SIMPLEMD_SEQUENTIAL_PATH} -l${LIBSIMPLEMD}"
fi

# specific built for SIMPLE_MD
if [ "${mdSim}" == "SIMPLE_MD" ]
then
  cd ${BUILD_PATH}
  ${compiler} ${MAMICO_PATH}/coupling/solvers/CoupledMolecularDynamicsSimulation.cpp ${FLAGS} ${includes} -c -o ${BUILD_PATH}/CoupledMolecularDynamicsSimulation.o
  objects="${objects} ${BUILD_PATH}/CoupledMolecularDynamicsSimulation.o"
elif [ "${mdSim}" == "LAMMPS_MD" ]
then
  if [ "${parallel}" == "sequential" ]
  then
    echo "ERROR: LAMMPS only works for option 'parallel'"
    exit -1
  fi
  # build lammps
  cd ${LAMMPS_PATH}
  make yes-user-mamico
  make openmpi
  make makelib
  make -f Makefile.lib openmpi
  cp liblammps_openmpi.a lib${LIBLAMMPS}.a
  includes="${includes} -I${LAMMPS_PATH}"
  libraries="${libraries} -L${LAMMPS_PATH} -l${LIBLAMMPS}"
  # remove Werror flag
  FLAGS=`echo ${FLAGS} | sed 's/-Werror//'`
elif [ "${mdSim}" == "LAMMPS_DPD" ]
then
  if [ "${parallel}" == "sequential" ]
  then
    echo "ERROR: LAMMPS only works for option 'parallel'"
    exit -1
  fi
  # build lammps
  cd ${LAMMPS_PATH}
  make yes-user-mamico
  make openmpi
  make makelib
  make -f Makefile.lib openmpi
  cp liblammps_openmpi.a lib${LIBLAMMPS}.a
  includes="${includes} -I${LAMMPS_PATH}"
  libraries="${libraries} -L${LAMMPS_PATH} -l${LIBLAMMPS}"
  # remove Werror flag
  FLAGS=`echo ${FLAGS} | sed 's/-Werror//'`
fi

### builds, linking, objects for coupled simulation with MaMiCo
cd ${BUILD_PATH}
${compiler} ${MAMICO_PATH}/coupling/configurations/ParticleInsertionConfiguration.cpp ${FLAGS} ${includes} -c -o ${BUILD_PATH}/ParticleInsertionConfiguration.o
${compiler} ${MAMICO_PATH}/coupling/tests/main_nie.cpp ${FLAGS} ${includes} -c -o ${BUILD_PATH}/main_nie.o
objects="${objects} ${BUILD_PATH}/ParticleInsertionConfiguration.o ${BUILD_PATH}/main_nie.o"

${compiler} ${objects} ${libraries} -o ${BUILD_PATH}/test

