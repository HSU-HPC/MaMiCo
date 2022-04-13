#!/usr/bin/env bash

#if any compilation fails, exit
set -e

for compiler in $1; do

	cmake . -D CMAKE_CXX_COMPILER=${compiler} -D BUILD_WITH_MPI=OFF || { #try/catch
		>&2 echo "CompileTest: Makefile generation of sequential CouetteTest failed for compiler: ${compiler}"
		exit 1 
	}

	make clean

	make || { #try/catch
		>&2 echo "CompileTest: Compilation of sequential CouetteTest failed for compiler: ${compiler}" 
		exit 1 
	}

	cmake . -D CMAKE_CXX_COMPILER=${compiler} -D BUILD_WITH_MPI=ON || { #try/catch
		>&2 echo "CompileTest: Makefile generation of parallel CouetteTest failed for compiler: ${compiler}" 
		exit 1 
	}

	make clean

	make || { #try/catch
	 	>&2 echo "CompileTest: Compilation of parallel CouetteTest failed for compiler: ${compiler}"
	 	exit 1 
	}

done

make clean

echo "CompileTest finished succesfully!"
