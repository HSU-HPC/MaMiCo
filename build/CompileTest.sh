#!/usr/bin/env bash

#if any compilation fails, exit
set -e

#TODO: 
#-Support more compilers
#-In case of failure, act accordingly. e.g. write error output to file

#TODO: Add more compilers
for compiler in g++ clang++ mpicxx icc; do

	cmake .. -D CMAKE_CXX_COMPILER=${compiler} || #try/catch
	echo "MaMiCo: CompileTest: Makefile generation failed for compiler: ${compiler}"

	make clean

	#TODO: build target parametrization?
	make || #try/catch
	echo "MaMiCo: CompileTest: Compilation failed for compiler: ${compiler}"

done

make clean


#Fragen:
# -eigene build dir für jeden compiler oder make clean nach jedem build?
# -alternative: nicht immer neue makefiles generieren, sondern separat speichern und erneut make clean && make
# 	-problem: wenn sich ordner-/filestruktur ändert
