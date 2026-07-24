if(NOT LAMMPS_DIR)
    message(FATAL_ERROR "Could not find the lammps directory. Please set the CMake variable LAMMPS_DIR to the correct location.")
endif()

find_path(LAMMPS_INCLUDE_DIR lammps.h HINTS ${LAMMPS_DIR} ${LAMMPS_DIR}/src)

find_library(LAMMPS_LIBRARY lammps HINTS ${LAMMPS_DIR} ${LAMMPS_DIR}/build)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LAMMPS
    DEFAULT_MSG
    LAMMPS_LIBRARY
    LAMMPS_INCLUDE_DIR)
