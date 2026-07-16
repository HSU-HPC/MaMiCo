find_path(LAMMPS_INCLUDE_DIR lammps.h HINTS ${LAMMPS_DIR} ${LAMMPS_DIR}/build ${LAMMPS_DIR}/build/includes)

find_library(LAMMPS_LIBRARY lammps HINTS ${LAMMPS_DIR} ${LAMMPS_DIR}/build)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LAMMPS
    DEFAULT_MSG
    LAMMPS_LIBRARY
    LAMMPS_INCLUDE_DIR)
