
The macro-micro-coupling tool for multiscale coupled molecular-continuum flow simulation.

## Documentation
* The source code documentation can be found [here](https://hsu-hpc.github.io/MaMiCo/)

## Requirements and optional dependencies
To compile and execute MaMiCo on a linux system, you need at least:
* a C++17 compiler installed (e.g. GCC 7 or newer)
* CMake version >= 3.0

Optional dependencies are:
* [MPI](https://www.open-mpi.org/), highly recommended for parallel runs, e.g. on clusters.
* [Eigen 3](http://eigen.tuxfamily.org/), a linear algebra library used for the POD noise filter.
* [pybind11](https://pybind11.readthedocs.io/en/stable/), for the MaMiCo Python bindings.
* [OpenFOAM](https://develop.openfoam.com/Development/openfoam/blob/develop/doc/Build.md), for coupling to CFD simulations with OpenFOAM.
* [preCICE](https://precice.org/), for coupling with other solvers using the preCICE coupling library.
* [ls1-MarDyn](https://www.ls1-mardyn.de/home.html), for coupling to MD simulations with ls1.
* [LCOV](https://github.com/linux-test-project/lcov), for test code coverage analysis.
* [genhtml](https://linux.die.net/man/1/genhtml), to generate HTML view from LCOV coverage data files.

## Build instructions
* First clone this repository and create a new build directory:

        git clone git@github.com:HSU-HPC/MaMiCo.git
        cd MaMiCo
        mkdir build
        cd build
* If you don't have a github account or don't want to use it, you can alternatively use `git clone https://github.com/HSU-HPC/MaMiCo.git` (read-only). 

* Then initialise CMake in your build directory, pointing it to the root directory of the repository. (Note that on some older systems, the command may be named `cmake3` instead of `cmake`.)

        cmake ..

* This will configure the build system with default settings and generate makefiles in the build directory. Optionally, you might want to change the build configuration now:

        ccmake .

* Here you can enable optional dependencies, e.g. activate BUILD_WITH_MPI (default OFF), or modifiy which compiler and flags you want to use. When you are done, press the configure [c] and generate [g] keys. Now you can start the compilation:

        make

### Additional instructions to build with ls1
* After cloning the repository, initialize the ls1 submodule with

        git submodule init
        git submodule update

* Follow the instructions on the [ls1 repository](https://github.com/ls1mardyn/ls1-mardyn) to build with cmake, however remember to enable the MAMICO_COUPLING flag, and provide the MaMiCo base directory in the MAMICO_SRC_DIR variable.

* Make MaMiCo as normal, choosing LS1_MARDYN as your MD library.

## Input file construction and first runs
* The build instructions have created an executable of the standard Couette flow test case, the file is called 'couette'. It expects an XML configuration file named 'couette.xml' in the current working directory. 
* An example simulation configuration file can be found [here](https://github.com/HSU-HPC/MaMiCo/blob/master/examples/couette.xml.template), other template input files are in the [examples](https://github.com/HSU-HPC/MaMiCo/blob/master/examples) folder as well. 
* Copy this file to your build folder, rename it 'couette.xml'. 
* The available options and features are listed directly in the template file via XML comments, so that you can modify the configuration to suit your needs. 
* Start the simulation by executing (sequential case) `./couette` or e.g. (MPI-parallel) `mpirun -n 8 ./couette`.
* If you get the error message 'ERROR MoleculeService::MoleculeService: Could not open file CheckpointSimpleMD_10000_reflecting_0.checkpoint!', copy the file of the same name from the 'examples' folder into your build folder.
* Depending on the configuration, you will obtain various output files in CSV, VTK or other formats. 

## Papers to cite
* P. Jarmatz, P. Neumann: [MaMiCo: Parallel Noise Reduction for Multi-instance Molecular-Continuum Flow Simulation](https://link.springer.com/chapter/10.1007/978-3-030-22747-0_34), International Conference on Computational Science. Springer, Cham, 2019
* P. Neumann, X. Bian: [MaMiCo: Transient multi-instance molecular-continuum flow simulation on supercomputers](https://doi.org/10.1016/j.cpc.2017.06.026), Computer Physics Communications 220 (2017): 390-402
* P. Neumann, H. Flohr, R. Arora, P. Jarmatz: [MaMiCo: Software design for parallel molecular-continuum flow simulations](https://doi.org/10.1016/j.cpc.2015.10.029) Computer Physics Communications 200 (2016): 324-335
