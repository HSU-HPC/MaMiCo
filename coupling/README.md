# About SensitivityRuns

This branch contains additional material for the research on coupling parameters. One may use it to deplicate the results, which were analysed. 

## First Steps First

1. Clone this repository
2. Install or compile OpenFOAM (best version 7) on your system
3. Set the correct path to your OpenFOAM library in the personal_settings file (located in the main folder)
4. Update the path in MaMiCo/coupling/tests/build_couette/couette.xml in line 80 according to your system

## Simulation Setup

1. Create the OpenFOAM mesh
   - Open MaMiCo/coupling/tests/build_couette/FoamSetup in the terminal
   - Execute 'blockMesh'
3. Compile the couette test: Go to MaMiCo/coupling/tests/build_couette and execute ./buildtest.sh parallel
3. If you want to run the default setup, do it, and go to 5.
4. If you want to run a different simulation, change the parameters in the couette.xml
5. Run the simulation by executing: ./test or mpirun -np [numberOfProcesses] ./test
