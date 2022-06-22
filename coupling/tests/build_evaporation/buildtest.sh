#!/bin/bash

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
BUILD_PATH=${MAMICO_PATH}/coupling/tests/build_evaporation;

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

compiler=""
libraries=""
#libraries="-ldmalloc -ldmallocxx" uncomment this if you use dmalloc for debugging purposes
objects=""
includes="-I${MAMICO_PATH}"

### specify flags, includes, libraries,compiler for parallel or sequential build
if [ "${parallel}" == "parallel" ]
then
    # note: we need to set MDDim3 for ALL Simulations since we use the configuration classes from SimpleMD
    FLAGS="-DSIMPLE_MD -DMDDim3 -std=c++1z -Werror -Wno-unknown-pragmas -Wno-int-in-bool-context -Wno-maybe-uninitialized -Wall -DMDCoupledParallel -DTarchParallel -DMPICH_IGNORE_CXX_SEEK -DENABLE_POST_MULTI_INSTANCE_FILTERING -O3"
    # -DMDCoupledDebug"
    includes="${includes} -I${MPI_INCLUDE_PATH} -I${LIB_EIGEN_PATH}"
    libraries="-L${MPI_LIB_PATH} -l${LIB_MPI}"
    compiler="mpicxx"
else
    FLAGS="-DSIMPLE_MD -DMDDim3 -std=c++1z -Wall -Wno-unknown-pragmas -DENABLE_POST_MULTI_INSTANCE_FILTERING -O3"
    # -Werror
    includes="${includes} -I${LIB_EIGEN_PATH}"
    compiler="g++"
fi
###

### specity flags, includes and libraries for a run with OpenFOAM
if [ -v FOAM_PATH ]
then
  echo "MaMiCo is compiled including OpenFOAM library"
  FLAGS="${FLAGS} -DBUILD_WITH_OPENFOAM -m64 -DOPENFOAM=2112 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -Wno-attributes -DNoRepository -ftemplate-depth-100"
  includes="${includes} -I${FOAM_PATH}src/finiteVolume/lnInclude -IlnInclude -I. -I${FOAM_PATH}src/meshTools/lnInclude -I${FOAM_PATH}src/OpenFOAM/lnInclude -I${FOAM_PATH}src/OSspecific/POSIX/lnInclude -I${FOAM_PATH}src/dynamicFvMesh/lnInclude -I${FOAM_PATH}src/thermophysicalModels/basic/lnInclude -I${FOAM_PATH}src/transportModels/compressible/lnInclude -I${FOAM_PATH}src/TurbulenceModels/turbulenceModels/lnInclude -I${FOAM_PATH}src/TurbulenceModels/compressible/lnInclude -I${FOAM_PATH}applications/solvers/compressible/rhoCentralFoam/BCs/rho -I${FOAM_PATH}applications/solvers/compressible/rhoCentralFoam"
  libraries="${libraries} -fPIC -fuse-ld=bfd -Xlinker --add-needed -Xlinker --no-as-needed -L${FOAM_PATH}platforms/linux64GccDPInt32Opt/lib -L${FOAM_PATH}platforms/linux64GccDPInt32Opt/lib/dummy -lfiniteVolume -lfvOptions -lmeshTools -lcompressibleTransportModels -lfluidThermophysicalModels -lspecie -lrhoCentralFoam -lturbulenceModels -lcompressibleTurbulenceModels -ldynamicFvMesh -ltopoChangerFvMesh -lOpenFOAM -ldl -lm "
else
  FLAGS="${FLAGS} -pedantic"
fi
### Removed flags -Dlinux64 -DWM_ARCH_OPTION=64 -g -pg -Wold-style-cast
### Removed Include Paths
### Removed Libraries -lPstream -lOpenFOAM -lsurfMesh -lfileFormats
### From rhoCentralFoam allwmake g++ -std=c++11 -m64 -pthread -DOPENFOAM=2112 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -Wno-attributes -Wno-unknown-pragmas  -O3  -DNoRepository -ftemplate-depth-100
### -IBCs/lnInclude -I/home/helene/openfoam/src/finiteVolume/lnInclude -I/home/helene/openfoam/src/meshTools/lnInclude -I/home/helene/openfoam/src/transportModels/compressible/lnInclude -I/home/helene/openfoam/src/thermophysicalModels/basic/lnInclude -I/home/helene/openfoam/src/thermophysicalModels/specie/lnInclude -I/home/helene/openfoam/src/TurbulenceModels/turbulenceModels/lnInclude -I/home/helene/openfoam/src/TurbulenceModels/compressible/lnInclude -I/home/helene/openfoam/src/dynamicFvMesh/lnInclude -iquote. -IlnInclude -I/home/helene/openfoam/src/OpenFOAM/lnInclude -I/home/helene/openfoam/src/OSspecific/POSIX/lnInclude   -fPIC -c rhoCentralFoam.C -o /home/helene/openfoam/build/linux64GccDPInt32Opt/applications/solvers/compressible/rhoCentralFoam/rhoCentralFoam.o
### g++ -std=c++11 -m64 -pthread -DOPENFOAM=2112 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -Wno-attributes -Wno-unknown-pragmas  -O3  -DNoRepository -ftemplate-depth-100 -IBCs/lnInclude -I/home/helene/openfoam/src/finiteVolume/lnInclude -I/home/helene/openfoam/src/meshTools/lnInclude -I/home/helene/openfoam/src/transportModels/compressible/lnInclude -I/home/helene/openfoam/src/thermophysicalModels/basic/lnInclude -I/home/helene/openfoam/src/thermophysicalModels/specie/lnInclude -I/home/helene/openfoam/src/TurbulenceModels/turbulenceModels/lnInclude -I/home/helene/openfoam/src/TurbulenceModels/compressible/lnInclude -I/home/helene/openfoam/src/dynamicFvMesh/lnInclude -iquote. -IlnInclude -I/home/helene/openfoam/src/OpenFOAM/lnInclude -I/home/helene/openfoam/src/OSspecific/POSIX/lnInclude   -fPIC -Xlinker --add-needed -Xlinker --no-as-needed  /home/helene/openfoam/build/linux64GccDPInt32Opt/applications/solvers/compressible/rhoCentralFoam/rhoCentralFoam.o -L/home/helene/openfoam/platforms/linux64GccDPInt32Opt/lib \
    ### -lfiniteVolume -lfvOptions -lmeshTools -lcompressibleTransportModels -lfluidThermophysicalModels -lspecie -lrhoCentralFoam -lturbulenceModels -lcompressibleTurbulenceModels -ldynamicFvMesh -ltopoChangerFvMesh -lOpenFOAM -ldl  \
    ### -lm -o /home/helene/openfoam/platforms/linux64GccDPInt32Opt/bin/rhoCentralFoam


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

cd ${BUILD_PATH}
${compiler} -c ${MAMICO_PATH}/coupling/solvers/CoupledMolecularDynamicsSimulation.cpp ${FLAGS} ${includes} -o ${BUILD_PATH}/CoupledMolecularDynamicsSimulation.o
objects="${objects} ${BUILD_PATH}/CoupledMolecularDynamicsSimulation.o"

### builds, linking, objects for coupled simulation with MaMiCo
cd ${BUILD_PATH}
${compiler} -c ${MAMICO_PATH}/coupling/configurations/ParticleInsertionConfiguration.cpp ${FLAGS} ${includes} -o ${BUILD_PATH}/ParticleInsertionConfiguration.o
${compiler} -c ${MAMICO_PATH}/coupling/tests/main_evaporation.cpp ${FLAGS} ${includes} -o ${BUILD_PATH}/main_evaporation.o
objects="${objects} ${BUILD_PATH}/ParticleInsertionConfiguration.o ${BUILD_PATH}/main_evaporation.o"

${compiler} ${objects} ${libraries} -o ${BUILD_PATH}/test
