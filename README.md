# About SensitivityRuns

The macro-micro-coupling tool for multiscale coupled molecular-continuum flow simulation.
This branch contains additional material for the research on coupling parameters. One may use it to deplicate the results, which were analysed. 

## Requirements and dependencies
* C++17 compiler

## First Steps First

1. Clone this repository
2. Install or compile OpenFOAM (best version 7) on your system
3. Set the correct path to your OpenFOAM library in the personal_settings file (located in the main folder)
4. Update the path in coupling/tests/build_couette/couette.xml in line 80 according to your system

## Simulation Setup

1. Create the OpenFOAM mesh
   - Open coupling/tests/build_couette/FoamSetup in the terminal
   - Execute 'blockMesh'
3. Compile the couette test: Go to coupling/tests/build_couette and execute ./buildtest.sh parallel
3. If you want to run the default setup, do it, and go to 5.
4. If you want to run a different simulation, change the parameters in the couette.xml
5. Run the simulation by executing: ./test or mpirun -np [numberOfProcesses] ./test
 
## Papers to cite
* P. Jarmatz, P. Neumann: [MaMiCo: Parallel Noise Reduction for Multi-instance Molecular-Continuum Flow Simulation](https://link.springer.com/chapter/10.1007/978-3-030-22747-0_34), International Conference on Computational Science. Springer, Cham, 2019
* P. Neumann, X. Bian: [MaMiCo: Transient multi-instance molecular-continuum flow simulation on supercomputers](https://doi.org/10.1016/j.cpc.2017.06.026), Computer Physics Communications 220 (2017): 390-402
* P. Neumann, H. Flohr, R. Arora, P. Jarmatz: [MaMiCo: Software design for parallel molecular-continuum flow simulations](https://doi.org/10.1016/j.cpc.2015.10.029) Computer Physics Communications 200 (2016): 324-335
