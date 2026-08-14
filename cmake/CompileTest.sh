#!/usr/bin/env bash

#if any compilation fails, exit
set -e

for build in $2; do
	for compiler in $1; do
		# cmake must be called twice if the compiler changes (the other variables are deleted, so we need to set them again later)
		cmake . -D CMAKE_CXX_COMPILER=${compiler} || { #try/catch
		>&2 echo "CompileTest: Makefile generation of CouetteTest failed for compiler: ${compiler} with MPI ${mpi}"
		exit 1 
		}

		for mpi in ON OFF; do
			echo "CompileTest: Running test with build type=${build}, compiler=${compiler}, MPI=${mpi}"
			start=$SECONDS
			
			cmake . -D CMAKE_CXX_COMPILER=${compiler} -D BUILD_WITH_MPI=${mpi} -D BUILD_WITH_PYBIND11=${mpi} -D CMAKE_BUILD_TYPE=${build} || { #try/catch
			>&2 echo "CompileTest: Makefile generation of CouetteTest failed for compiler: ${compiler} with MPI ${mpi}"
			exit 1 
			}

			make -j4 || { #try/catch
				>&2 echo "CompileTest: Compilation of CouetteTest failed for compiler: ${compiler} with MPI ${mpi}" 
				exit 1 
			}

			elapsed=$((SECONDS - start))
			make clean # Clean up after timing
    		echo "CompileTest: Test with build type=${build}, compiler=${compiler}, MPI=${mpi} took ${elapsed}s"
		done
	done
done

make clean

echo "CompileTest finished succesfully!"