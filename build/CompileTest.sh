#!/usr/bin/env bash

#if any compilation fails, exit
set -e

#TODO: 
#-Support more compilers
#-In case of failure, act accordingly. e.g. write error output to file
#-Support parallel/sequential for CouetteTest

#TODO: Add more compilers
for compiler in g++ clang++ mpicxx; do

	#COUETTE
	cmake .. -D CMAKE_CXX_COMPILER=${compiler} || #try/catch
	echo "MaMiCo: CompileTest: Makefile generation of CouetteTest failed for compiler: ${compiler}"

	make clean

	make || #try/catch
	echo "MaMiCo: CompileTest: Compilation of CouetteTest failed for compiler: ${compiler}"

	#PYBIND11
	cmake .. -D CMAKE_CXX_COMPILER=${compiler} -D USE_PYBIND11=True || #try/catch
	echo "MaMiCo: CompileTest: Makefile generation for pybind11 library failed for compiler: ${compiler}"

	make clean

	make || #try/catch
	echo "MaMiCo: CompileTest: Compilation of pybind11 library failed for compiler: ${compiler}"

done

make clean
