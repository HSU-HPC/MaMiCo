if(NOT DEFINED OUTPUT_FILE)
    message(FATAL_ERROR "Usage: -DOUTPUT_FILE=<path>")
endif()

set(HASH "unknown")

find_package(Git QUIET)

if(GIT_FOUND)
    execute_process(
        COMMAND ${GIT_EXECUTABLE} rev-parse --is-inside-work-tree
        RESULT_VARIABLE INSIDE_REPO
        OUTPUT_QUIET
        ERROR_QUIET
    )

    if(INSIDE_REPO EQUAL 0)
        execute_process(
            COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
            OUTPUT_VARIABLE HASH
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )

        execute_process(
            COMMAND ${GIT_EXECUTABLE} status --porcelain
            OUTPUT_VARIABLE STATUS
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )

        if(NOT "${STATUS}" STREQUAL "")
            set(HASH "${HASH}-dirty")
        endif()
    endif()
endif()
message(STATUS ">> MAMICO_COMMIT_HASH: " ${HASH})

# Avoid re-build if hash is unchanged
file(WRITE "${OUTPUT_FILE}.tmp"
"#pragma once\n#define MAMICO_COMMIT_HASH \"${HASH}\"\n")
execute_process(
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${OUTPUT_FILE}.tmp"
            "${OUTPUT_FILE}"
)

file(REMOVE "${OUTPUT_FILE}.tmp")
