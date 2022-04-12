// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_TESTS_TESTLAMMPS_H_
#define _MOLECULARDYNAMICS_COUPLING_TESTS_TESTLAMMPS_H_

#include "coupling/CouplingMDDefinitions.h"
#include "coupling/configurations/MaMiCoConfiguration.h"
#include "coupling/interface/MamicoInterfaceProvider.h"
#include "coupling/interface/impl/macroscopictestsolvers/VoidMacroscopicSolverInterface.h"
#include "coupling/services/MacroscopicCellService.h"
#include "coupling/tests/Test.h"
#include "tarch/configuration/ParseConfiguration.h"
#include <sstream>

#include "atom.h"
#include "input.h"
#include "lammps.h"
#include "mamico_cell.h"
#include "mpi.h"

/** general test class for lammps. Here, we provide tests for most get..(..)
 * methods of the MDSolverInterface as well as initialisation routines which can
 * be re-used by other inherited test cases.
 *  @author Philipp Neumann
 */
template <unsigned int dim> class TestLammps : public Test {
public:
  TestLammps(int argc, char** argv, std::string name)
      : Test(name), _lammps(NULL), _macroSolverInterface(NULL), _macroscopicCellService(NULL), _argc(argc), _argv(argv) {}

  /** if existent, deletes test configuration */
  virtual ~TestLammps() {
    if (_macroscopicCellService != NULL) {
      delete _macroscopicCellService;
      _macroscopicCellService = NULL;
    }
    if (_macroSolverInterface != NULL) {
      delete _macroSolverInterface;
      _macroSolverInterface = NULL;
    }
    coupling::interface::MamicoInterfaceProvider<LAMMPS_NS::MamicoCell, dim>::getInstance().setMacroscopicSolverInterface(NULL);
    coupling::interface::MamicoInterfaceProvider<LAMMPS_NS::MamicoCell, dim>::getInstance().setMDSolverInterface(NULL);
    coupling::interface::MamicoInterfaceProvider<LAMMPS_NS::MamicoCell, dim>::getInstance().setMacroscopicCellService(NULL);
    if (_lammps != NULL) {
      delete _lammps;
      _lammps = NULL;
    }
  }

  /** runs various tests on get-methods of the MD interface. Explicitly tested
   * are: getGlobalMDDomainSize() getGlobalMDDomainOffset() getMoleculeMass() ->
   * here we only allow molecules/atoms to be of type "1" getKB() -> is
   * always 1.0 (LJ case) getMoleculeSigma() -> is always 1.0 (LJ case)
   *  getMoleculeEpsilon() -> is always 1.0 (LJ case)
   *  getDt()
   *  Implicitly, we also test the method getLinkedCell() since this one is
   * automatically invoked in the setup phase of the macroscopic cell service;
   *  We can hence check the output from the couts in a first place.
   */
  virtual void run() {
    double tmp = -1.0;                          // test variable
    tarch::la::Vector<dim, double> tmpVec(0.0); // test vector
    bool valid = true;
    std::cout << "Load lammps test configuration..." << std::endl;
    if (dim == 2) {
      loadLammpsTestConfiguration("inputpositionsonly2D.xyz", 4);
    } else {
      loadLammpsTestConfiguration("inputpositionsonly3D.xyz", 8);
    }
    std::cout << "Load void macroscopic solver configuration..." << std::endl;
    loadMacroscopicSolverConfiguration();
    std::cout << "Load MaMiCo configuration and init macroscopic cell service..." << std::endl;
    loadMamicoTestConfiguration();
    // test molecule mass
    tmp = coupling::interface::MamicoInterfaceProvider<LAMMPS_NS::MamicoCell, dim>::getInstance().getMDSolverInterface()->getMoleculeMass();
    std::cout << "Molecule mass: " << tmp << std::endl;
    if (tmp != 2.0) {
      std::cout << "Wrong mass value, mass=" << tmp << ", instead of 2.0!" << std::endl;
      exit(EXIT_FAILURE);
    }
    // test Boltzmann constant kB
    tmp = coupling::interface::MamicoInterfaceProvider<LAMMPS_NS::MamicoCell, dim>::getInstance().getMDSolverInterface()->getKB();
    std::cout << "Boltzmann constant: " << tmp << std::endl;
    if (tmp != 1.0) {
      std::cout << "Wrong boltzmann constant, kB=" << tmp << ", should be 1.0!" << std::endl;
      exit(EXIT_FAILURE);
    }
    // TODO put tests for sigma+epsilon here
    // test time step dt
    tmp = coupling::interface::MamicoInterfaceProvider<LAMMPS_NS::MamicoCell, dim>::getInstance().getMDSolverInterface()->getDt();
    std::cout << "Time step: " << tmp << std::endl;
    if (tmp != 1.0e-14) {
      std::cout << "Wrong time step, time step=" << tmp << ", should be 1.0e-14!" << std::endl;
      exit(EXIT_FAILURE);
    }
    // test domain size
    tmpVec = coupling::interface::MamicoInterfaceProvider<LAMMPS_NS::MamicoCell, dim>::getInstance().getMDSolverInterface()->getGlobalMDDomainSize();
    valid = (tmpVec[0] == 10.0);
    for (int d = 1; d < dim; d++) {
      valid = valid && (tmpVec[d] == 10.0);
    }
    if (!valid) {
      std::cout << "Wrong domain size " << tmpVec << ", should be 10.0 in each component!" << std::endl;
      exit(EXIT_FAILURE);
    }
    // test domain offset
    tmpVec = coupling::interface::MamicoInterfaceProvider<LAMMPS_NS::MamicoCell, dim>::getInstance().getMDSolverInterface()->getGlobalMDDomainOffset();
    valid = (tmpVec[0] == 1.0);
    for (int d = 1; d < dim; d++) {
      valid = valid && (tmpVec[d] == 1.0);
    }
    if (!valid) {
      std::cout << "Wrong domain offset " << tmpVec << ", should be 1.0 in each component!" << std::endl;
      exit(EXIT_FAILURE);
    }
  }

protected:
  /** loads the LAMMPS test configuration. We provide a respective script in
   * input->one()-statements. Besides, we parse input positions for molecules
   *  from a file inputpositionsonly2D.xyz or inputpositions3D.xyz,
   * respectively. This method supports dim=2 and dim=3. The arguments
   * correspond to the file containing the molecule setup and the number of
   * atoms contained.
   */
  void loadLammpsTestConfiguration(std::string inputmoleculeconfiguration, int numberAtoms) {
    int size = -1;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if ((dim != 2) && (dim != 3)) {
      std::cout << "ERROR TestLammps::loadLammpsTestConfiguration: only dim=2 "
                   "or dim=3 supported!"
                << std::endl;
      exit(EXIT_FAILURE);
    }
    // allocate lammps
    if (_lammps != NULL) {
      delete _lammps;
      _lammps = NULL;
    }
    _lammps = new LAMMPS_NS::LAMMPS(_argc, _argv, MPI_COMM_WORLD);
    if (_lammps == NULL) {
      std::cout << "ERROR TestLammps::loadLammpsTestConfiguration: Could not "
                   "allocate _lammps!"
                << std::endl;
      exit(EXIT_FAILURE);
    }
    // read input -> this should also initialise the md solver interface
    //_lammps->input->file();
    _lammps->input->one("units lj");
    _lammps->input->one("atom_style atomic");
    _lammps->input->one("processors * * * map xyz");
    if (dim == 2) {
      _lammps->input->one("dimension 2");
    }
    if (dim == 3) {
      _lammps->input->one("region box block 1.0 11.0 1.0 11.0 1.0 11.0");
    } else {
      _lammps->input->one("region box block 1.0 11.0 1.0 11.0 0.0 10.0");
    }
    _lammps->input->one("create_box 1 box");
    if (dim == 2) {
      std::stringstream ss;
      ss << "create_atoms 1 random " << numberAtoms << " 123 box";
      _lammps->input->one(ss.str().c_str());
      ss.str("");
      ss.clear();
      ss << "read_dump " << inputmoleculeconfiguration << " 0 x y replace yes box no format xyz";
      _lammps->input->one(ss.str().c_str());
    } else {
      std::stringstream ss;
      ss << "create_atoms 1 random " << numberAtoms << " 123 box";
      _lammps->input->one(ss.str().c_str());
      ss.str("");
      ss.clear();
      ss << "read_dump " << inputmoleculeconfiguration << " 0 x y z replace yes box no format xyz";
      _lammps->input->one(ss.str().c_str());
    }
    _lammps->input->one("mass 1 2.0");
    _lammps->input->one("velocity all create 1.8 12345678 loop geom");
    _lammps->input->one("pair_style lj/cut 2.5");
    _lammps->input->one("pair_modify shift yes");
    _lammps->input->one("pair_coeff 1 1 1.0 1.0 2.5");
    // use a very small time step to yield nearly no position changes after 1 or
    // 2 time steps
    _lammps->input->one("timestep 1.0e-14");
    _lammps->input->one("neighbor 0.5 bin");
    _lammps->input->one("neigh_modify delay 0 every 20 check no");
    //_lammps->input->one("dump myDump all custom 100 atom*.dat id type x y z vx
    // vy vz"); _lammps->input->one("restart 100 checkpoint1 checkpoint2");
    _lammps->input->one("fix 1 all nve");
    // always use a minimal number of local mamico cells in fix mamico. This
    // implies that we have to choose this depending on the number of procs
    if (dim == 2) {
      switch (size) {
      case 1:
        _lammps->input->one("fix 2 all mamico 36 123456 2.5");
        break;
      case 4:
        _lammps->input->one("fix 2 all mamico 16 123456 2.5");
        break;
      case 16:
        _lammps->input->one("fix 2 all mamico  9 123456 2.5");
        break;
      default:
        std::cout << "ERROR TestLammps: only 1,4,16 procs supported in 2D!" << std::endl;
        exit(EXIT_FAILURE);
        break;
      }
    } else if (dim == 3) {
      switch (size) {
      case 1:
        _lammps->input->one("fix 2 all mamico 216 123456 2.5");
        break;
      case 8:
        _lammps->input->one("fix 2 all mamico  64 123456 2.5");
        break;
      case 64:
        _lammps->input->one("fix 2 all mamico  27 123456 2.5");
        break;
      default:
        std::cout << "ERROR TestLammps: only 1,8,64 procs supported in 3D!" << std::endl;
        exit(EXIT_FAILURE);
        break;
      }
    }
  }

  /** load mamico configuration and initialse macroscopic cell service; should
   * be called after loadLammpsTestConfiguration() and
   * loadMacroscopicSolverConfiguration(). We assume a valid configuration file
   *  for this test setup in the file mamico_lammps_test_configuration.xml.
   */
  void loadMamicoTestConfiguration() {
    std::string mamicoLammpsTestConfiguration = "";
    if (dim == 2) {
      mamicoLammpsTestConfiguration = "mamico_lammps_test_configuration2D.xml";
    } else if (dim == 3) {
      mamicoLammpsTestConfiguration = "mamico_lammps_test_configuration3D.xml";
    }
    tarch::la::Vector<dim, unsigned int> numberProcesses(0);
    int rank = 0;
    int size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (dim == 2) {
      switch (size) {
      case 1:
        for (int d = 0; d < dim; d++) {
          numberProcesses[d] = 1;
        }
        break;
      case 4:
        for (int d = 0; d < dim; d++) {
          numberProcesses[d] = 2;
        }
        break;
      case 16:
        for (int d = 0; d < dim; d++) {
          numberProcesses[d] = 4;
        }
        break;
      default:
        std::cout << "ERROR TestLammps: only 1,4,16 procs supported in 2D!" << std::endl;
        exit(EXIT_FAILURE);
        break;
      }
    } else if (dim == 3) {
      switch (size) {
      case 1:
        for (int d = 0; d < dim; d++) {
          numberProcesses[d] = 1;
        }
        break;
      case 8:
        for (int d = 0; d < dim; d++) {
          numberProcesses[d] = 2;
        }
        break;
      case 64:
        for (int d = 0; d < dim; d++) {
          numberProcesses[d] = 4;
        }
        break;
      default:
        std::cout << "ERROR TestLammps: only 1,8,64 procs supported in 3D!" << std::endl;
        exit(EXIT_FAILURE);
        break;
      }
    } else {
      std::cout << "ERROR TestLammps::loadMamicoTestConfiguration: only 2D and "
                   "3D supported!"
                << std::endl;
      exit(EXIT_FAILURE);
    }
    const unsigned int numberMDTimestepsPerCouplingCycle = 100000;
    coupling::configurations::MaMiCoConfiguration<dim> config;
    std::cout << "Parse config " << mamicoLammpsTestConfiguration << std::endl;
    tarch::configuration::ParseConfiguration::parseConfiguration<coupling::configurations::MaMiCoConfiguration<dim>>(mamicoLammpsTestConfiguration, "mamico",
                                                                                                                     config);
    std::cout << "Init macroscopic cell service..." << std::endl;
    if (_macroscopicCellService != NULL) {
      delete _macroscopicCellService;
      _macroscopicCellService = NULL;
    }
    _macroscopicCellService = new coupling::services::MacroscopicCellServiceImpl<LAMMPS_NS::MamicoCell, dim>(
        0, coupling::interface::MamicoInterfaceProvider<LAMMPS_NS::MamicoCell, dim>::getInstance().getMDSolverInterface(),
        coupling::interface::MamicoInterfaceProvider<LAMMPS_NS::MamicoCell, dim>::getInstance().getMacroscopicSolverInterface(), numberProcesses, rank,
        config.getParticleInsertionConfiguration(), config.getMomentumInsertionConfiguration(), config.getTransferStrategyConfiguration(),
        config.getParallelTopologyConfiguration(), numberMDTimestepsPerCouplingCycle, config.getMacroscopicCellConfiguration());
    if (_macroscopicCellService == NULL) {
      std::cout << "ERROR TestLammps: _macroscopicCellService==NULL!" << std::endl;
      exit(EXIT_FAILURE);
    }
    std::cout << "Set macroscopic cell service in MamicoInterfaceProvider..." << std::endl;
    coupling::interface::MamicoInterfaceProvider<LAMMPS_NS::MamicoCell, dim>::getInstance().setMacroscopicCellService(_macroscopicCellService);
  }

  /** this loads some dummy solver interface which basically does nothing at
   * all. */
  void loadMacroscopicSolverConfiguration() {
    // set up dummy macroscopic solver interface
    if (_macroSolverInterface != NULL) {
      delete _macroSolverInterface;
      _macroSolverInterface = NULL;
    }
    _macroSolverInterface = new coupling::interface::VoidMacroscopicSolverInterface<dim>();
    if (_macroSolverInterface == NULL) {
      std::cout << "ERROR TestLammps: macroSolverInterface==NULL!" << std::endl;
      exit(EXIT_FAILURE);
    }
    coupling::interface::MamicoInterfaceProvider<LAMMPS_NS::MamicoCell, dim>::getInstance().setMacroscopicSolverInterface(_macroSolverInterface);
  }

  /** prints the molecule information to the screen. If "printGhostMolecules" is
   * true, also ghost atom information is printed. */
  void printMolecules(bool printGhostMolecules = false) {
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Barrier(MPI_COMM_WORLD);

    int numberAtoms = _lammps->atom->nlocal;
    if (printGhostMolecules)
      numberAtoms += _lammps->atom->nghost;

    if (rank == 0)
      std::cout << "Print molecule information: " << std::endl;
    for (int r = 0; r < size; r++) {
      if (r == rank) {
        for (int i = 0; i < numberAtoms; i++) {
          std::cout << "Rank " << rank << ": molecule " << i << ": " << _lammps->atom->x[i][0] << "," << _lammps->atom->x[i][1] << "," << _lammps->atom->x[i][2]
                    << std::endl;
        }
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }

  /** test configuration for LAMMPS runs */
  LAMMPS_NS::LAMMPS* _lammps;
  /** test solver interface */
  coupling::interface::VoidMacroscopicSolverInterface<dim>* _macroSolverInterface;
  coupling::services::MacroscopicCellService<dim>* _macroscopicCellService;
  /** arguments from commandline */
  int _argc;
  char** _argv;
};

#endif // _MOLECULARDYNMAMICS_COUPLING_TESTS_TESTLAMMPS_H_
