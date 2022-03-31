// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_TESTS_TESTESPRESSO2DUMMYCOUPLING_H_
#define _MOLECULARDYNAMICS_COUPLING_TESTS_TESTESPRESSO2DUMMYCOUPLING_H_
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <sstream>
#include <string>
#include <tcl.h>

#include "tarch/la/Vector.h"

// Espresso Files
#include "communication.hpp"
#include "global_tcl.hpp"
#include "initialize.hpp"
#include "integrate.hpp"
#include "integrate_tcl.hpp"
#include "interaction_data_tcl.hpp"
#include "particle_data.hpp"
#include "thermostat_tcl.hpp"
#include "utils.hpp"

// MaMiCo Files
#include "MaMiCo.h"
#include "coupling/interface/impl/Espresso/EspressoMDMolecule.h"
#include "coupling/interface/impl/Espresso/EspressoMDMoleculeIterator.h"
#include "coupling/interface/impl/Espresso/EspressoMDSolverInterface.h"
#include "coupling/solvers/DummySolverInterfaceService.h"
#include "coupling/tests/Test.h"

/** This Test Class provides an implementation of Espresso to Dummy Solver
 * coupling.
 *  Note that this test class is not inherited from Test Espresso and is
 * implemented separately
 *  @author Rahul Arora
 */

// Forward declaration of tcl commands, from src/tcl folder
void register_tcl_commands(Tcl_Interp *interp);
void register_global_variables(Tcl_Interp *interp);

class TestEspresso2DummyCoupling : public Test {
public:
  TestEspresso2DummyCoupling(std::string name, int argc, char **argv) : Test(name), _argc(argc), _argv(argv) {}

  virtual ~TestEspresso2DummyCoupling(){};

  virtual void run() {
    // Load Coupled Test Configuration
    std::cout << "Loading Coupled Configuration " << std::endl;
    loadCoupledTestConfiguration();
  }

private:
  int _argc;
  char **_argv;

protected:
  void loadCoupledTestConfiguration() {
    // Installation of tcl commands
    Tcl_Interp *interp;
    interp = Tcl_CreateInterp();
    if (TCL_OK != Tcl_Init(interp)) {
      std::cout << "Tcl_Init error: " << Tcl_GetStringResult(interp) << std::endl;
      exit(EXIT_FAILURE);
    }
    register_tcl_commands(interp);
    register_global_variables(interp);

    // Espresso Initialization
    mpi_init(&_argc, &_argv);
    on_program_start();

    /* Flag that is used to break the mpi_loop() in Espresso communications.cpp,
   in order to run the coupling test, I require a start stop
        * functionality in Espresso, however in espresso the slave processors are
   always in a loop and wait for instructions from the master 		* processor,
   thus we have to break this loop in order!!*/
    break_flag = 1;

    // I do initilaization of Espresso here and then stop Espresso such that all
    // the data remain on the particular processor
    if (this_node == 0) {

      // Initialize the cell system
      int cellsystemflag = Tcl_Eval(interp, "cellsystem domain_decomposition -no_verlet_list;");
      if (cellsystemflag == 1) {
        std::cout << "ERROR: cellsystem definition " << std::endl;
      }

      // Set global variables using setmd
      int globalvarsflag = Tcl_Eval(interp, "setmd time_step 0.002;") || Tcl_Eval(interp, "setmd skin      0.0;") ||
                           Tcl_Eval(interp, "setmd box_l 30 30 30;") || Tcl_Eval(interp, "setmd max_num_cells 3456;") ||
                           Tcl_Eval(interp, "setmd periodic 1 1 1;");
      if (globalvarsflag == 1) {
        std::cout << "ERROR: global variables definition " << std::endl;
      }

      // Set the thermostat
      // Espresso's thermostat is switched of during this calculation, however
      // there is a possibilty to use the thermostat provided.
      // int thermostatflag = Tcl_Eval(interp, "thermostat langevin 1.8 0.5");
      int thermostatflag = Tcl_Eval(interp, "thermostat off");
      if (thermostatflag == 1) {
        std::cout << "ERROR: Thermostat definition " << std::endl;
      }

      // Set up LJ interactions
      int interflag = Tcl_Eval(interp, "inter 0 0 lennard-jones 1.0 1.0 1.12246 auto;");
      if (interflag == 1) {
        std::cout << "ERROR: interactions definition " << std::endl;
      }

      // Read particle positions and velocity from input file
      std::stringstream ss;
      std::string line;
      std::string inputFile;
      inputFile = "espresso_input_positions.config";
      std::ifstream in(inputFile.c_str());
      if (!in.is_open()) {
        std::cout << "ERROR init(): Could not open " << inputFile << "!" << std::endl;
        exit(EXIT_FAILURE);
      }
      while (in.good() && (getline(in, line))) {
        std::istringstream iss(line);
        int id;
        double x, y, z, vx, vy, vz, fx, fy, fz;
        iss >> id >> x >> y >> z >> vx >> vy >> vz;
        ss.str("");
        ss << "part " << id << " pos " << x << " " << y << " " << z << " v " << vx << " " << vy << " " << vz << ";";
        int partflag = Tcl_Eval(interp, ss.str().c_str());
        if (partflag == 1) {
          std::cout << "ERROR: reading particle data from file " << std::endl;
        }
      }
      in.close();
      std::cout << "The cell grid is " << dd.cell_grid[0] << " " << dd.cell_grid[1] << " " << dd.cell_grid[2] << std::endl;
      mpi_stop();
    } else {
      mpi_loop();
    }

    /****************************************************************************************************************************/
    /***************************************** Coupled Integration
     * **************************************************************/
    /****************************************************************************************************************************/

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

#ifdef MAMICO_DEBUG
    int count = 0;
    for (unsigned int i = 0; i < local_cells.n; i++) {
      ParticleList *list = local_cells.cell[i];
      count += list->n;
    }

    std::cout << "Rank " << rank << " Total number of particles is " << n_part << " and no of particles on this node are " << count
              << " and the number of cells on this node are " << local_cells.n << std::endl;

    std::cout << "Node Grid " << node_grid[0] << " " << node_grid[1] << " " << node_grid[2] << " My node is " << this_node << " Node position " << node_pos[0]
              << " " << node_pos[1] << " " << node_pos[2] << std::endl;
#endif

    // Initialize MaMiCo the coupling tool, called on all processors
    initialize();

    // We need a barrier here, so that all processors are synchronized
    MPI_Barrier(MPI_COMM_WORLD);

    // Add information to buffers for send operation from Dummy Solver to MD
    if (this_node == 0) {
      std::cout << "Velocity and density profile before coupling " << std::endl;
      tarch::la::Vector<3, unsigned int> loop0(0);
      loop0[0] = 9;
      loop0[2] = 9;
      for (loop0[1] = 0; loop0[1] < 18; loop0[1]++) {
        tarch::la::Vector<3, double> tempVel = _dummySolver.getVelocity(loop0[0], loop0[1], loop0[2]);
        std::cout << loop0[1] << " " << tempVel[0] << " " << _dummySolver.getDensity(loop0[0], loop0[1], loop0[2]) << std::endl;
      }
      storeDummySolverDataInSendBuffer();
    }

    // Send Dummy Solver information to MD
    _macroscopicCellService->sendFromMacro2MD(DummySolverInterfaceService::getInstance().getSendBuffer(),
                                              DummySolverInterfaceService::getInstance().getGlobalCellIndices4SendBuffer());

    int outflag;

    // Integration loop, restarting espresso
    // Inside this loop, I also calculate and write the KE, TE and Temp to
    // espresso_output.obs for analysis of simulation
    for (unsigned int i = 0; i < 20; i++) {
      // Call to macroscopic cell service to sample, distribute momentum and
      // apply temperature
      _macroscopicCellService->processInnerMacroscopicCellAfterMDTimestep();
      // TODO: in future (e.g., for tests with momentum transfer other than
      // velocity grad. relaxation): move distributeMomentum after force
      // computation
      _macroscopicCellService->distributeMass(i);
      _macroscopicCellService->distributeMomentum(i);
      _macroscopicCellService->applyTemperatureToMolecules(i);

      break_flag = 1;
      integrate_vv(1);

      MPI_Barrier(MPI_COMM_WORLD);

      if (i % 2 == 0) {
        if (this_node == 0) {
          std::cout << "Iteration Number " << i << std::endl;
          if (total_energy.init_status == 0) {
            init_energies(&total_energy);
            energy_calc((double *)total_energy.data.e);
          }
          double KE = total_energy.data.e[0];
          std::cout << "Temperature is " << KE / (1.5 * n_part) << std::endl;
          /*
                                        double value = total_energy.data.e[0];
                                        for (int j = 1; j < total_energy.data.n; j++){
                                        value += total_energy.data.e[j];
                                        }

                                        for (int j = 0; j < n_external_potentials; j++) {
                                        value+=external_potentials[j].energy;
                                        }
          */
          std::ofstream outfile;
          if (i == 0) {
            outfile.open("espresso_output.obs", std::ios::out);
            // outfile << "time\ttemp\ttotal energy\tkinetic energy\n";
            outfile << "time\ttemp\tkinetic energy\n";
            outfile.close();
          }
          outfile.open("espresso_output.obs", std::ios::out | std::ios::app);
          // outfile << i*time_step << "\t" << KE/(1.5*n_part) << "\t" << value
          //<< "\t" << KE << "\n";
          outfile << i * time_step << "\t" << KE / (1.5 * n_part) << "\t" << KE << "\n";
          outfile.close();
        } else {
          energy_calc(NULL);
        }
      }
    }

    // Send MD (after coupled iterations) information to Dummy Solver
    _macroscopicCellService->sendFromMD2Macro(DummySolverInterfaceService::getInstance().getReceiveBuffer(),
                                              DummySolverInterfaceService::getInstance().getGlobalCellIndices4ReceiveBuffer());

    break_flag = 0;

    // Write particles positions and velocities to espresso_output_$rank.config,
    // for Testing purposes
    std::stringstream filename;
    filename.str("");
    filename << "espresso_output_" << rank << ".config";
    std::ofstream outfile;
    outfile.open(filename.str().c_str(), std::ios::out);

    for (unsigned int i = 0; i < local_cells.n; i++) {
      ParticleList *list = local_cells.cell[i];
      Particle *part;
      part = list->part;
      for (unsigned int j = 0; j < list->n; j++) {
        outfile << part[j].p.identity << "  " << part[j].r.p[0] << " " << part[j].r.p[1] << " " << part[j].r.p[2] << "  " << part[j].m.v[0] << " "
                << part[j].m.v[1] << " " << part[j].m.v[2] << "\n";
      }
    }
    outfile.close();

    MPI_Barrier(MPI_COMM_WORLD);

    // Print out the velocity values in the dummy solver after coupling to check
    // whether coupling was correct or not
    if (this_node == 0) {
      writeReceiveBufferDataToDummySolver();

      std::cout << std::endl;
      std::cout << "Velocity and density profile after coupling" << std::endl;
      tarch::la::Vector<3, unsigned int> loop(0);
      loop[0] = 9;
      loop[2] = 9;
      for (loop[1] = 0; loop[1] < 18; loop[1]++) {
        tarch::la::Vector<3, double> tempVel = _dummySolver.getVelocity(loop[0], loop[1], loop[2]);
        std::cout << loop[1] << " " << tempVel[0] << " " << _dummySolver.getDensity(loop[0], loop[1], loop[2]) << std::endl;
      }
    }

    // Shutdown MaMiCo
    shutdown();

    break_flag = 0;
    if (this_node == 0) {
    } else {
      mpi_loop();
    }
  }
};

#endif // _MOLECULARDYNAMICS_COUPLING_TESTS_TESTESPRESSO2DUMMYCOUPLING_H_
