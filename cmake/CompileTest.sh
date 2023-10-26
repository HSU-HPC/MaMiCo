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
			cmake . -D CMAKE_CXX_COMPILER=${compiler} -D BUILD_WITH_MPI=${mpi} -D CMAKE_BUILD_TYPE=${build} || { #try/catch
			>&2 echo "CompileTest: Makefile generation of CouetteTest failed for compiler: ${compiler} with MPI ${mpi}"
			exit 1 
			}

			make clean

			make -j4 || { #try/catch
				>&2 echo "CompileTest: Compilation of CouetteTest failed for compiler: ${compiler} with MPI ${mpi}" 
				exit 1 
			}
		done
	done
done

make clean

echo "CompileTest finished succesfully!"