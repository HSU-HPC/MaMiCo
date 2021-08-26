#!/bin/bash
set -v
 
### local settings like path variables
SETTINGS=../../../personal_settings

if test -f "$SETTINGS"; then
	source ../../../personal_settings
else
	echo "ERROR! No personal settings file found at $SETTINGS ."
	exit -1
fi

# path variables automatically set using the ones from above
BUILD_PATH=${MAMICO_PATH}/coupling/tests/build_espresso;
ESPRESSO_PATH_SRC=${ESPRESSO_PATH}/src;
buildscriptsdir=${ESPRESSO_PATH}/scripts;
DEFINES="-D ESPRESSO_SCRIPTS_DEFAULT=\"${buildscriptsdir}\" -DMDCoupledParallel -std=c++11";
ESPRESSOTCL_PATH=${ESPRESSO_PATH_SRC}/tcl;
ESPRESSO_LIB_PATH=${ESPRESSO_PATH_SRC}/libEspresso.la;
ESPRESSOTCL_LIB_PATH=${ESPRESSOTCL_PATH}/libEspressoTcl.la;
COMPILER=mpic++;
LIBTOOL=libtool;

rm ${BUILD_PATH}/test_espresso;
rm ${BUILD_PATH}/*.o;

${COMPILER} ${MAMICO_PATH}/coupling/configurations/ParticleInsertionConfiguration.cpp ${DEFINES} -I${MPI_INLCUDE_PATH} -I${MAMICO_PATH} -c -o ${BUILD_PATH}/ParticleInsertionConfiguration.o
${COMPILER} ${MAMICO_PATH}/tarch/tinyxml2/tinyxml2.cpp ${DEFINES} -I${MPI_INLCUDE_PATH} -I${MAMICO_PATH} -c -o ${BUILD_PATH}/tinyxml2.o
${COMPILER} ${MAMICO_PATH}/tarch/utils/RandomNumberService.cpp ${DEFINES} -I${MPI_INLCUDE_PATH} -I${MAMICO_PATH} -c -o ${BUILD_PATH}/RandomNumberService.o
${COMPILER} ${MAMICO_PATH}/coupling/solvers/DummySolverInterfaceService.cpp ${DEFINES} -I${MPI_INLCUDE_PATH} -I${MAMICO_PATH}  -I${ESPRESSO_PATH_SRC} -c -o ${BUILD_PATH}/DummySolverInterfaceService.o
${COMPILER} ${ESPRESSOTCL_PATH}/scriptsdir.cpp ${DEFINES} -I${MPI_INLCUDE_PATH} -I${MAMICO_PATH} -I${ESPRESSO_PATH_SRC} -I${ESPRESSOTCL_PATH} -c -o ${BUILD_PATH}/scriptsdir.o
${COMPILER} -g ${MAMICO_PATH}/coupling/tests/main_espresso.cpp ${DEFINES} -I${TCL_INLCUDE_PATH} -I${MPI_INLCUDE_PATH} -I${MAMICO_PATH} -I${ESPRESSOTCL_PATH} -I${ESPRESSO_PATH_SRC} -c -o ${BUILD_PATH}/main_espresso.o
${LIBTOOL} --tag=CXX --mode=link ${COMPILER} -g ${BUILD_PATH}/tinyxml2.o ${BUILD_PATH}/RandomNumberService.o ${BUILD_PATH}/ParticleInsertionConfiguration.o ${BUILD_PATH}/DummySolverInterfaceService.o ${BUILD_PATH}/scriptsdir.o  ${BUILD_PATH}/main_espresso.o -L${TCL_LIB_PATH} -l${LIB_TCL} -L${MPI_LIB_PATH} -l${LIB_MPI} ${ESPRESSOTCL_LIB_PATH} ${ESPRESSO_LIB_PATH} -o ${BUILD_PATH}/test_espresso
