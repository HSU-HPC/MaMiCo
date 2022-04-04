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

compiler=""
libraries=""
#libraries="-ldmalloc -ldmallocxx" uncomment this if you use dmalloc for debugging purposes
objects=""
includes="-I${MAMICO_PATH}"

### specify flags, includes, libraries,compiler for parallel or sequential build
if [ "${parallel}" == "parallel" ]
then
    if [ -v MARDYN_PATH ]
    then
      FLAGS="-DLS1_MARDYN -DMAMICO_COUPLING -DMARDYN_AUTOPAS -DMDDim3 -std=c++17 -Wall -Wfatal-errors -Wno-unknown-pragmas -O3 -DMARDYN_DPDP -DENABLE_MPI -DMDCoupledParallel -DTarchParallel -DMPICH_IGNORE_CXX_SEEK -Wno-sign-conversion -Wno-gnu-array-member-paren-init" # todo put -Werror, O0 for debug
      includes="${includes} -I${LIB_EIGEN_PATH} -I${MARDYN_PATH}/src -I${MARDYN_PATH}/libs/rapidxml -I${MARDYN_PATH}/build/_deps/autopasfetch-src/src -I${MARDYN_PATH}/build/_deps/spdlog-src/include"
      libraries="-L${MPI_LIB_PATH} -l${LIB_MPI}"
      compiler="mpicxx"
      #FLAGS="-DLS1_MARDYN -DMAMICO_COUPLING -DMDDim3 -std=c++17 -Wall -Wfatal-errors -Wno-unknown-pragmas -O3 -DMARDYN_DPDP -DENABLE_MPI -DMDCoupledParallel -DTarchParallel -DMPICH_IGNORE_CXX_SEEK -Wno-sign-conversion -Wno-gnu-array-member-paren-init" # todo put -Werror, O0 for debug
      #includes="${includes} -I${LIB_EIGEN_PATH} -I${MARDYN_PATH}/src -I${MARDYN_PATH}/libs/rapidxml -I${MARDYN_PATH}/build/_deps/spdlog-src/include"
      #libraries="-L${MPI_LIB_PATH} -l${LIB_MPI}"
      #compiler="mpicxx"
    else
      # note: we need to set MDDim3 for ALL Simulations since we use the configuration classes from SimpleMD
      FLAGS="-DSIMPLE_MD -DMDDim3 -std=c++1z -Werror -Wno-unknown-pragmas -Wno-int-in-bool-context -Wno-maybe-uninitialized -Wall -Wfatal-errors -DMDCoupledParallel -DTarchParallel -DMPICH_IGNORE_CXX_SEEK -DENABLE_POST_MULTI_INSTANCE_FILTERING -O3"
      # -DMDCoupledDebug"
      includes="${includes} -I${MPI_INCLUDE_PATH} -I${LIB_EIGEN_PATH}"
      libraries="-L${MPI_LIB_PATH} -l${LIB_MPI}"
      compiler="mpicxx"
    fi
else
    if [ -v MARDYN_PATH ]
    then
      FLAGS="-DLS1_MARDYN -DMAMICO_COUPLING -DMARDYN_AUTOPAS -DMDDim3 -std=c++17 -Wall -Wfatal-errors -Wno-unknown-pragmas -O3 -DMARDYN_DPDP -Wno-gnu-array-member-paren-init" # todo put -Werror
      includes="${includes} -I${LIB_EIGEN_PATH} -I${MARDYN_PATH}/src -I${MARDYN_PATH}/libs/rapidxml -I${MARDYN_PATH}/build/_deps/autopasfetch-src/src -I${MARDYN_PATH}/build/_deps/spdlog-src/include"
      compiler="g++"
    else
      FLAGS="-DSIMPLE_MD -DMDDim3 -std=c++1z -Wall -Wfatal-errors -Wno-unknown-pragmas -O3" 
      # -Werror
      includes="${includes} -I${LIB_EIGEN_PATH} "
      compiler="g++"
    fi
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
  if [ -v MARDYN_PATH ]
  then
    scons compiler=clang target=libsimplemd dim=3 build=release parallel=yes -j4
    libraries="${libraries} -L${SIMPLEMD_PARALLEL_PATH} -l${LIBSIMPLEMD}"
    FLAGS="${FLAGS} -DMDParallel"

     ### ls1 mardyn
      libraries="${libraries} -L${MARDYN_PATH}/build/src -l:libMarDyn.a"
    ### autopas
      libraries="${libraries} -L${MARDYN_PATH}/build/_deps/autopasfetch-build/src/autopas -l:libautopas.a"
    ### spdlog (for ls1)
      libraries="${libraries} -L${MARDYN_PATH}/build/_deps/spdlog-build -l:libspdlog.a -lpthread" #in debug, libspdlogd, else libspdlog
    ### harmony (for autopas)
      libraries="${libraries} -L${MARDYN_PATH}/build/_deps/autopasfetch-build/libs/harmony/include/lib -l:libharmony.a"
  else
    scons target=libsimplemd dim=3 build=release parallel=yes -j4
    libraries="${libraries} -L${SIMPLEMD_PARALLEL_PATH} -l${LIBSIMPLEMD}"
    FLAGS="${FLAGS} -DMDParallel"
  fi
else
  if [ -v MARDYN_PATH ]
  then
    scons compiler=clang target=libsimplemd dim=3 build=release parallel=no -j4
    libraries="${libraries} -L${SIMPLEMD_SEQUENTIAL_PATH} -l${LIBSIMPLEMD}"
    ### tinyxml2

    ### ls1 mardyn
      libraries="${libraries} -L${MARDYN_PATH}/build/src -l:libMarDyn.a"
    ### autopas
      libraries="${libraries} -L${MARDYN_PATH}/build/_deps/autopasfetch-build/src/autopas -l:libautopas.a"
    ### spdlog (for ls1)
      libraries="${libraries} -L${MARDYN_PATH}/build/_deps/spdlog-build -l:libspdlog.a -lpthread" #in debug, libspdlogd, else libspdlog
    ### harmony (for autopas)
      libraries="${libraries} -L${MARDYN_PATH}/build/_deps/autopasfetch-build/libs/harmony/include/lib -l:libharmony.a"
  else
    scons target=libsimplemd dim=3 build=release parallel=no -j4
    libraries="${libraries} -L${SIMPLEMD_SEQUENTIAL_PATH} -l${LIBSIMPLEMD}"
  fi
fi

#FLAGS="${FLAGS} -pthread"

#FLAGS="${FLAGS} -g"
#FLAGS="${FLAGS} -fsanitize=address"
#libraries="${libraries} -fsanitize=address"
libraries="${libraries} -lxerces-c"

cd ${BUILD_PATH}
if ! [ -v MARDYN_PATH ]; then
  ${compiler} ${MAMICO_PATH}/coupling/solvers/CoupledMolecularDynamicsSimulation.cpp ${FLAGS} ${includes} -c -o ${BUILD_PATH}/CoupledMolecularDynamicsSimulation.o
  objects="${objects} ${BUILD_PATH}/CoupledMolecularDynamicsSimulation.o"
fi

### builds, linking, objects for coupled simulation with MaMiCo
cd ${BUILD_PATH}
${compiler} -c ${MAMICO_PATH}/coupling/configurations/ParticleInsertionConfiguration.cpp ${FLAGS} ${includes} -o ${BUILD_PATH}/ParticleInsertionConfiguration.o
${compiler} -c ${MAMICO_PATH}/coupling/tests/main_couette.cpp ${FLAGS} ${includes} -o ${BUILD_PATH}/main_couette.o
objects="${objects} ${BUILD_PATH}/ParticleInsertionConfiguration.o ${BUILD_PATH}/main_couette.o"

${compiler} ${objects} ${libraries} -o ${BUILD_PATH}/test
