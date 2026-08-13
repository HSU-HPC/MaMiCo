
set(ALLRepoPath https://gitlab.version.fz-juelich.de/SLMS/loadbalancing.git)

FetchContent_Declare(
        allfetch
        GIT_REPOSITORY ${ALLRepoPath}
        GIT_TAG v0.9.4
)

FetchContent_GetProperties(allfetch)
if (NOT allfetch_POPULATED)
    FetchContent_Populate(allfetch)
endif()

find_path(ALL_INCLUDE_DIR ALL.hpp HINTS ${allfetch_SOURCE_DIR}/include/)
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ALL
    DEFAULT_MSG
    ALL_INCLUDE_DIR)
