#!/usr/bin/env bash

#if any compilation fails, exit
set -e

if (( $# != 4 )); then
    echo Usage: $0 '"${COMPILETEST_COMPILERS}" "${COMPILETEST_MODES}" "${COMPILETEST_MPI_MODES}" "${COMPILETEST_KOKKOS_TARGETS}"' >&2
    exit 1
fi

for build in $2; do
	for compiler in $1; do
		for kokkos_target in $4; do
			# cmake must be called twice if the compiler changes (the other variables are deleted, so we need to set them again later)
			cmake . -D CMAKE_CXX_COMPILER=${compiler} -D KOKKOS_TARGET=${kokkos_target} || { #try/catch
			>&2 echo "CompileTest: Makefile generation of CouetteTest failed for compiler/kokkos target: ${compiler}/${kokkos_target} with MPI ${mpi}"
			exit 1 
			}

			make clean

			for mpi in $3; do
				test_type="build type=${build}, kokkos target=${kokkos_target}, compiler=${compiler}, MPI=${mpi}"
				echo "CompileTest: Running test with ${test_type}"
				start=$SECONDS

				cmake \
					-D CMAKE_CXX_COMPILER=${compiler} \
					-D BUILD_WITH_MPI=${mpi} \
					-D BUILD_WITH_PYBIND11=${mpi} \
					-D CMAKE_BUILD_TYPE=${build} \
					-D KOKKOS_TARGET=${kokkos_target} \
					. || { #try/catch
				>&2 echo "CompileTest: Makefile generation of CouetteTest failed for ${test_type}"
				exit 1 
				}

				make -j4 || { #try/catch
					>&2 echo "CompileTest: Compilation of CouetteTest failed for  ${test_type}" 
					exit 1 
				}

				elapsed=$((SECONDS - start))
				make clean # Clean up after timing
				echo "CompileTest: Test with ${test_type} took ${elapsed}s"
				echo
			done
		done
	done
done

echo "CompileTest finished succesfully!"