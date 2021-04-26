#!/bin/bash

### include and library paths for MPI and name of MPI library
MPI_INCLUDE_PATH=/usr/lib/x86_64-linux-gnu/openmpi/include
MPI_LIB_PATH=/usr/lib/x86_64-linux-gnu/openmpi/lib
LIB_MPI=mpi

LIB_EIGEN_PATH=/usr/include/eigen3/

### Include OpenFOAMpath
FOAM_PATH=/opt/openfoam7/

### home directory of MAMICO
MAMICO_PATH=/home/helene/Dokumente/mamico-dev
### local settings like path variables
SETTINGS=../../../personal_settings

if test -f "$SETTINGS"; then
	source ${SETTINGS}
else
	echo "ERROR! No personal settings file found at $SETTINGS ."
	exit -1
fi

### build directory for library of SIMPLE_MD (currently specified for gnu compiler (intel variant: .../icc/..)
SIMPLEMD_PARALLEL_PATH=${MAMICO_PATH}/build/libsimplemd/release/dim3/parallel_yes/gcc/gprof_no/
SIMPLEMD_SEQUENTIAL_PATH=${MAMICO_PATH}/build/libsimplemd/release/dim3/parallel_no/gcc/gprof_no/
### name of lib for SIMPLE_MD
LIBSIMPLEMD=simplemd

### FROM HERE ON, ALL VARIABLES ARE AUTOMATICALLY DETERMINED ###

### build path for CouetteTest
BUILD_PATH=${MAMICO_PATH}/coupling/tests/build_couette;

parallel=$1;

if [ "${parallel}" == "parallel" ] || [ "${parallel}" == "sequential" ]
then
    echo "Build mode: ${parallel}"
else
    echo "ERROR! ./buildtest parallel/sequential"
    exit -1
fi

rm ${BUILD_PATH}/test;
rm ${BUILD_PATH}/*.o;
rm ${BUILD_PATH}/*.txt;
rm ${BUILD_PATH}/*.vtk;
rm ${BUILD_PATH}/*.csv;

compiler=""
libraries=""
objects=""
includes="-I${MAMICO_PATH}"

### specify flags, includes, libraries,compiler for parallel or sequential build
if [ "${parallel}" == "parallel" ]
then
    # note: we need to set MDDim3 for ALL Simulations since we use the configuration classes from SimpleMD
    FLAGS="-DSIMPLE_MD -DMDDim3 -std=c++1z -Werror -Wno-unknown-pragmas -Wno-int-in-bool-context -Wall -DMDCoupledParallel -DTarchParallel -DMPICH_IGNORE_CXX_SEEK -O3"
    # -DMDCoupledDebug"
    includes="${includes} -I${MPI_INCLUDE_PATH} -I${LIB_EIGEN_PATH}"
    libraries="-L${MPI_LIB_PATH} -l${LIB_MPI}"
    compiler="mpicxx"
else
    FLAGS="-DSIMPLE_MD -DMDDim3 -std=c++1z -Wall -Wno-unknown-pragmas -O3"
    # -Werror
    includes="${includes} -I${LIB_EIGEN_PATH}"
    compiler="g++"
fi
###

### specity flags, includes and libraries for a run with OpenFOAM
if [ -v FOAM_PATH ]
then
  echo "MaMiCo is compiled including OpenFOAM library"
  FLAGS="${FLAGS} -DBUILD_WITH_OPENFOAM -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -Wno-attributes -DNoRepository -ftemplate-depth-100 -g -pg"
  includes="${includes} -I${FOAM_PATH}src/finiteVolume/lnInclude -I${FOAM_PATH}src/meshTools/lnInclude -IlnInclude -I. -I${FOAM_PATH}src/OpenFOAM/lnInclude -I${FOAM_PATH}src/OSspecific/POSIX/lnInclude"
  libraries="${libraries} -fPIC -fuse-ld=bfd -Xlinker --add-needed -Xlinker --no-as-needed -L${FOAM_PATH}platforms/linux64GccDPInt32Opt/lib -L${FOAM_PATH}platforms/linux64GccDPInt32Opt/lib/dummy -lfiniteVolume -lmeshTools -lOpenFOAM -ltriSurface -lPstream -lsurfMesh -lfileFormats -ldl -lm"
else
  FLAGS="${FLAGS} -pedantic"
fi
###

### builds, objects, libraries for coupling -> we require several parts from simplemd
cd ${MAMICO_PATH}
if [ "${parallel}" == "parallel" ]
then
        /usr/bin/env python $(which scons) target=libsimplemd dim=3 build=release parallel=yes -j4
        libraries="${libraries} -L${SIMPLEMD_PARALLEL_PATH} -l${LIBSIMPLEMD}"
        FLAGS="${FLAGS} -DMDParallel"
else
	/usr/bin/env python $(which scons) target=libsimplemd dim=3 build=release parallel=no -j4
        libraries="${libraries} -L${SIMPLEMD_SEQUENTIAL_PATH} -l${LIBSIMPLEMD}"
fi

cd ${BUILD_PATH}
${compiler} ${MAMICO_PATH}/coupling/solvers/CoupledMolecularDynamicsSimulation.cpp ${FLAGS} ${includes} -c -o ${BUILD_PATH}/CoupledMolecularDynamicsSimulation.o
objects="${objects} ${BUILD_PATH}/CoupledMolecularDynamicsSimulation.o"

### builds, linking, objects for coupled simulation with MaMiCo
cd ${BUILD_PATH}
${compiler} ${MAMICO_PATH}/coupling/configurations/ParticleInsertionConfiguration.cpp ${FLAGS} ${includes} -c -o ${BUILD_PATH}/ParticleInsertionConfiguration.o
${compiler} ${MAMICO_PATH}/coupling/tests/main_couette.cpp ${FLAGS} ${includes} -c -o ${BUILD_PATH}/main_couette.o
objects="${objects} ${BUILD_PATH}/ParticleInsertionConfiguration.o ${BUILD_PATH}/main_couette.o"

${compiler} ${objects} ${libraries} -o ${BUILD_PATH}/test
