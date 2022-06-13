
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
* [OpenFOAM](https://openfoam.org/), for coupling to CFD simulations with OpenFOAM.
* [preCICE](https://precice.org/), for coupling with other solvers using the preCICE coupling library.

## Build instructions and first steps
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

* The previous step has created an executable of the standard Couette flow test case, the file is called 'couette'. It expects an XML configuration file named 'couette.xml' in the current working directory. An example simulation configuration can be found [here](https://github.com/HSU-HPC/MaMiCo/blob/master/examples/couette.xml.template).
* Copy this file to your build folder, rename it 'couette.xml' and modify the configuration to suit your needs.
* Start the simulation by executing (sequential case) `./couette` or e.g. (MPI-parallel) `mpirun -n 8 ./couette`. 
* Depending on the configuration, you will obtain various output files in CSV, VTK or other formats. 

## Papers to cite
* P. Jarmatz, P. Neumann: [MaMiCo: Parallel Noise Reduction for Multi-instance Molecular-Continuum Flow Simulation](https://link.springer.com/chapter/10.1007/978-3-030-22747-0_34), International Conference on Computational Science. Springer, Cham, 2019
* P. Neumann, X. Bian: [MaMiCo: Transient multi-instance molecular-continuum flow simulation on supercomputers](https://doi.org/10.1016/j.cpc.2017.06.026), Computer Physics Communications 220 (2017): 390-402
* P. Neumann, H. Flohr, R. Arora, P. Jarmatz: [MaMiCo: Software design for parallel molecular-continuum flow simulations](https://doi.org/10.1016/j.cpc.2015.10.029) Computer Physics Communications 200 (2016): 324-335
