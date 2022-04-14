// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_TESTS_TESTESPRESSO_H_
#define _MOLECULARDYNAMICS_COUPLING_TESTS_TESTESPRESSO_H_
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "communication.hpp"
#include "global_tcl.hpp"
#include "initialize.hpp"
#include "integrate.hpp"
#include "integrate_tcl.hpp"
#include "interaction_data_tcl.hpp"
#include "thermostat_tcl.hpp"
#include "utils.hpp"
#include <tcl.h>

#include "coupling/tests/Test.h"

/** general test class for Espresso. Here, I provide method to initialisation routines which can be re-used by other inherited test classes.
 *  In this class, I provide an example of how to start Espresso calculations by parsing Tcl commands, which is later used in
 *  Dummy solver to Espresso test
 *  @author Rahul Arora
 */

// Forward declaration of tcl commands, from src/tcl folder
void register_tcl_commands(Tcl_Interp* interp);
void register_global_variables(Tcl_Interp* interp);
int tclcommand_cellsystem(ClientData data, Tcl_Interp* interp, int argc, char** argv);
int tclcommand_blockfile(ClientData data, Tcl_Interp* interp, int argc, char* argv[]);
int tclcommand_part(ClientData data, Tcl_Interp* interp, int argc, char** argv);

class TestEspresso : public Test {
public:
  TestEspresso(std::string name, int argc, char** argv) : Test(name), _argc(argc), _argv(argv) {}

  virtual ~TestEspresso(){};

  virtual void run() {
    // Load Espresso Test Configuration
    std::cout << "Loading Espresso Configuration " << std::endl;
    loadEspressoTestConfiguration();
  }

private:
  int _argc;
  char** _argv;

protected:
  void loadEspressoTestConfiguration() {
    // Installation of tcl commands
    Tcl_Interp* interp;
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

    // break_flag enables us to have a start-stop functionality in Espresso (for further details see espresso-mamico.pdf)
    break_flag = 0;

    // All commands will be executed on the master processor, while the other processors wait for instructions from master processor
    if (this_node == 0) {
      // Initialize the cell system
      int cellsystemflag = Tcl_Eval(interp, "cellsystem domain_decomposition -no_verlet_list;");
      if (cellsystemflag == 1) {
        std::cout << "ERROR: cellsystem definition " << std::endl;
      }

      // Set global variables using setmd
      int globalvarsflag = Tcl_Eval(interp, "setmd time_step 0.002;") || Tcl_Eval(interp, "setmd skin      0.0;") || Tcl_Eval(interp, "setmd box_l 2 2 2;") ||
                           Tcl_Eval(interp, "setmd max_num_cells 8;");
      if (globalvarsflag == 1) {
        std::cout << "ERROR: global variables definition " << std::endl;
      }

      // Set the thermostat
      int thermostatflag = Tcl_Eval(interp, "thermostat off ;");
      if (thermostatflag == 1) {
        std::cout << "ERROR: Thermostat definition " << std::endl;
      }

      // Set up LJ interactions
      int interflag = Tcl_Eval(interp, "inter 0 0 lennard-jones 1.0 1.0 1.0 auto;");
      int interflag2 = Tcl_Eval(interp, "inter forcecap 0;");
      if (interflag == 1 || interflag2 == 1) {
        std::cout << "ERROR: interactions definition " << std::endl;
      }

      // Particle data
      int partflag = Tcl_Eval(interp, "part 0 pos 0.5 0.5 0.5 type 0") || Tcl_Eval(interp, "part 1 pos 0.5 0.5 1.5 type 0") ||
                     Tcl_Eval(interp, "part 2 pos 0.5 1.5 1.5 type 0") || Tcl_Eval(interp, "part 3 pos 0.5 1.5 0.5 type 0") ||
                     Tcl_Eval(interp, "part 4 pos 1.5 0.5 0.5 type 0") || Tcl_Eval(interp, "part 5 pos 1.5 0.5 1.5 type 0") ||
                     Tcl_Eval(interp, "part 6 pos 1.5 1.5 1.5 type 0") || Tcl_Eval(interp, "part 7 pos 1.5 1.5 0.5 type 0");
      if (partflag == 1) {
        std::cout << " ERROR: setting up particle properties" << std::endl;
      }

      // Integrate command
      int integrateflag = Tcl_Eval(interp, "integrate 100");
      if (integrateflag == 1) {
        std::cout << "ERROR: integrate command " << std::endl;
      }

      // Writing back the positions and velocities to test.config
      std::stringstream ss;
      int outflag = Tcl_Eval(interp, "set trajectory [open \"test.config\" \"w\"]");
      outflag = Tcl_Eval(interp, "flush $trajectory");
      for (unsigned int i = 0; i < 8; i++) {
        ss.str("");
        ss << "puts $trajectory \"" << i << " [part " << i << " print pos] [part " << i << " print v]\"";
        outflag = Tcl_Eval(interp, ss.str().c_str());
      }
      outflag = Tcl_Eval(interp, "close $trajectory");
      if (outflag == 1) {
        std::cout << "ERROR: integrate command " << std::endl;
      }
    } else {
      mpi_loop();
    }
  }
};

#endif // _MOLECULARDYNAMICS_COUPLING_TESTS_TESTESPRESSO_H_
