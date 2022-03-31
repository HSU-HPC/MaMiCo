// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_TESTS_TESTESPRESSOMOLECULEITERATOR_H_
#define _MOLECULARDYNAMICS_COUPLING_TESTS_TESTESPRESSOMOLECULEITERATOR_H_
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>

#include "utils.hpp"
#include <tcl.h>
#include "initialize.hpp"
#include "communication.hpp"
#include "integrate.hpp"
#include "global_tcl.hpp"
#include "interaction_data_tcl.hpp"
#include "integrate_tcl.hpp"
#include "thermostat_tcl.hpp"
#include "particle_data.hpp"

#include "coupling/interface/impl/Espresso/EspressoMDMoleculeIterator.h"
#include "coupling/tests/Test.h"

/** This class is used to test the iterator over molecules in Espresso.
 *  I used the same input file as for LAMMPS
 * inputpositions3D_moleculeiterator.xyz
 *  @author Rahul Arora
 */

void register_tcl_commands(Tcl_Interp *interp);
void register_global_variables(Tcl_Interp *interp);
int tclcommand_cellsystem(ClientData data, Tcl_Interp *interp, int argc,
                          char **argv);
int tclcommand_blockfile(ClientData data, Tcl_Interp *interp, int argc,
                         char *argv[]);
int tclcommand_part(ClientData data, Tcl_Interp *interp, int argc, char **argv);

class TestEspressoMoleculeIterator : public Test {
private:
  int _argc;
  char **_argv;

public:
  TestEspressoMoleculeIterator(std::string name, int argc, char **argv)
      : Test(name), _argc(argc), _argv(argv) {}
  virtual ~TestEspressoMoleculeIterator() {}

  virtual void run() {
    // Initialize the Tcl interpreter to read tcl commands in C++
    Tcl_Interp *interp;
    interp = Tcl_CreateInterp();
    if (TCL_OK != Tcl_Init(interp)) {
      std::cout << "Tcl_Init error: " << Tcl_GetStringResult(interp)
                << std::endl;
      exit(EXIT_FAILURE);
    }
    register_tcl_commands(interp);
    register_global_variables(interp);

    // Espresso Initialization
    mpi_init(&_argc, &_argv);
    on_program_start();

    break_flag = 0;

    if (this_node == 0) {
      // Initialize the cell system
      int cellsystemflag =
          Tcl_Eval(interp, "cellsystem domain_decomposition -no_verlet_list;");
      if (cellsystemflag == 1) {
        std::cout << "ERROR: cellsystem definition " << std::endl;
      }

      // Set global variables using setmd
      int globalvarsflag = Tcl_Eval(interp, "setmd time_step 0.002;") ||
                           Tcl_Eval(interp, "setmd skin      0.0;") ||
                           Tcl_Eval(interp, "setmd box_l 10 10 10;") ||
                           Tcl_Eval(interp, "setmd max_num_cells 64;");
      if (globalvarsflag == 1) {
        std::cout << "ERROR: global variables definition " << std::endl;
      }

      // Set the thermostat
      int thermostatflag = Tcl_Eval(interp, "thermostat off ;");
      if (thermostatflag == 1) {
        std::cout << "ERROR: Thermostat definition " << std::endl;
      }

      // Set up LJ interactions
      int interflag =
          Tcl_Eval(interp, "inter 0 0 lennard-jones 1.0 1.0 1.12246 auto;");
      int interflag2 = Tcl_Eval(interp, "inter forcecap 0;");
      if (interflag == 1 || interflag2 == 1) {
        std::cout << "ERROR: interactions definition " << std::endl;
      }

      // Read Particle data from file
      std::stringstream ss;
      std::string line;
      std::string inputFile;
      inputFile = "testpositions_espresso_iterator.config";
      std::ifstream in(inputFile.c_str());
      if (!in.is_open()) {
        std::cout << "ERROR init(): Could not open " << inputFile << "!"
                  << std::endl;
        exit(EXIT_FAILURE);
      }
      while (in.good() && (getline(in, line))) {
        std::istringstream iss(line);
        int id;
        double x, y, z, vx, vy, vz;
        iss >> id >> x >> y >> z >> vx >> vy >> vz;
        ss.str("");
        ss << "part " << id << " pos " << x << " " << y << " " << z << " v "
           << vx << " " << vy << " " << vz << ";";
        int partflag = Tcl_Eval(interp, ss.str().c_str());
        if (partflag == 1) {
          std::cout << "ERROR: reading particle data from file " << std::endl;
        }
      }
      in.close();

      // Test EspressoMDMoleculeInterator
      testMolecules();
    } else {
      mpi_loop();
    }
  }

private:
  /** loops over all macroscopic cells on this process and checks for the number
   * of molecules in each cell */
  void testMolecules() {
    unsigned int flag = 1;
    for (unsigned int i = 0; i < local_cells.n; i++) {
      unsigned int numberOfMoleculesFound = 0;
      ParticleList *list = local_cells.cell[i];
      coupling::interface::EspressoMDMoleculeIterator it(*list);
      Particle *part;
      part = list->part;
      unsigned int numberOfMoleculesInCurrentCell = list->n;
      for (it.begin(); it.continueIteration(); it.next()) {
        numberOfMoleculesFound++;
      }
      if (numberOfMoleculesFound != numberOfMoleculesInCurrentCell) {
        flag = 0;
        std::cout << "ERROR: The number of molecules in the current cell " << i
                  << " are expected to be " << numberOfMoleculesInCurrentCell
                  << " while our implementation of molecule iterator returns "
                  << numberOfMoleculesFound << std::endl;
      }
    }
    if (flag == 1) {
      std::cout << "The test was successful " << std::endl;
    } else {
      std::cout << "The test was not successful " << std::endl;
    }
  }
};

#endif // _MOLECULARDYNAMICS_COUPLING_TESTS_TESTESPRESSOMOLECULEITERATOR_H_
