#!/bin/bash

source personal_settings


### build directory for library of SIMPLE_MD (currently specified for gnu compiler (intel variant: .../icc/..)
SIMPLEMD_PARALLEL_PATH=${MAMICO_PATH}/build/libsimplemd/release/dim3/parallel_yes/gcc/gprof_no/
### name of lib for SIMPLE_MD
LIBSIMPLEMD=simplemd

### FROM HERE ON, ALL VARIABLES ARE AUTOMATICALLY DETERMINED ###

### build path for DynamicMDTest
BUILD_PATH=${MAMICO_PATH}/coupling/tests/build_dynamic_md;

mdSim=$1;

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

### specify flags, includes, libraries,compiler for parallel build
# note: we need to set MDDim3 for ALL Simulations since we use the configuration classes from SimpleMD
FLAGS="-D${mdSim} -DMDDim3 -std=c++1z -pedantic -Werror -Wno-unknown-pragmas -Wall -DMDCoupledParallel -DTarchParallel -DMPICH_IGNORE_CXX_SEEK -O0 -g3"
# -DMDCoupledDebug"
includes="${includes} -I${MPI_INCLUDE_PATH} -I${LIB_EIGEN_PATH}"
libraries="-L${MPI_LIB_PATH} -l${LIB_MPI}"
compiler="mpicxx"

### 

### builds, objects, libraries for coupling -> we require several parts from simplemd
cd ${MAMICO_PATH} || exit
scons target=libsimplemd dim=3 build=release parallel=yes compiler=gcc machine=hww-cluster -j4
libraries="${libraries} -L${SIMPLEMD_PARALLEL_PATH} -l${LIBSIMPLEMD}"
FLAGS="${FLAGS} -DMDParallel"


# specific built for SIMPLE_MD
if [ "${mdSim}" == "SIMPLE_MD" ]
then
  cd ${BUILD_PATH} || exit 
    ${compiler} ${MAMICO_PATH}/coupling/solvers/CoupledMolecularDynamicsSimulation.cpp ${FLAGS} ${includes} -c -o ${BUILD_PATH}/CoupledMolecularDynamicsSimulation.o
    objects="${objects} ${BUILD_PATH}/CoupledMolecularDynamicsSimulation.o"
elif [ "${mdSim}" == "LAMMPS_MD" ]
then
  if [ "${parallel}" == "sequential" ]
  then
    echo "ERROR: LAMMPS only works for option 'parallel'"
    exit -1
  fi
  
  cp -a ${MAMICO_PATH}/coupling/interface/impl/LAMMPS/USER-MAMICO ${LAMMPS_PATH}/src/

  export MAMICO_PATH
  export LIB_EIGEN_PATH

  # build lammps
  cd ${LAMMPS_PATH}/src || exit
  make yes-user-mamico
  make g++_openmpi || exit
  #make makelib
  #make -f Makefile.lib openmpi
  ln -s liblammps_g++_openmpi.a lib${LIB_LAMMPS}.a
  includes="${includes} -I${LAMMPS_PATH}/src"
  libraries="${libraries} -L${LAMMPS_PATH}/src -l${LIB_LAMMPS}"
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
  cd ${LAMMPS_PATH} || exit
  make yes-user-mamico
  make openmpi
  make makelib
  make -f Makefile.lib openmpi
  ln -s liblammps_g++_openmpi.a lib${LIB_LAMMPS}.a
  includes="${includes} -I${LAMMPS_PATH}/src"
  libraries="${libraries} -L${LAMMPS_PATH}/src -Wl,-Bstatic -l${LIB_LAMMPS} -Wl,-Ddynamic"
  # remove Werror flag
  FLAGS=`echo ${FLAGS} | sed 's/-Werror//'`
fi

### builds, linking, objects for coupled simulation with MaMiCo
cd ${BUILD_PATH} || exit
${compiler} ${MAMICO_PATH}/coupling/configurations/ParticleInsertionConfiguration.cpp ${FLAGS} ${includes} -c -o ${BUILD_PATH}/ParticleInsertionConfiguration.o
${compiler} ${MAMICO_PATH}/coupling/tests/main_dynamic_md.cpp ${FLAGS} ${includes} -c -o ${BUILD_PATH}/main_dynamic_md.o
objects="${objects} ${BUILD_PATH}/ParticleInsertionConfiguration.o ${BUILD_PATH}/main_dynamic_md.o"

echo "${compiler} ${objects} ${libraries} -o ${BUILD_PATH}/test" 
${compiler} ${objects} ${libraries} -o ${BUILD_PATH}/test

