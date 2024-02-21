#pragma once

#include "USER-MAMICO/mamico_cell.h"
#include "USER-MAMICO/mamico_lammps_md_solver_interface.h"
#include "coupling/interface/MDSimulation.h"
#include "lammps/input.h"
#include "lammps/lammps.h"
#include <mpi.h>

namespace coupling {
namespace interface {

/** md simulation based on LAMMPS, encapsulated in this simulation interface */
class LammpsMDSimulation : public coupling::interface::MDSimulation {
private:
  LAMMPS_NS::LAMMPS* _lmp;
  const simplemd::configurations::MolecularDynamicsConfiguration& _configuration;
  const coupling::configurations::MaMiCoConfiguration<MDSIMULATIONFACTORY_DIMENSION>& _mamicoConfiguration;
  const double _tolerance; // tolerance for checks of configs

public:
  LammpsMDSimulation(const simplemd::configurations::MolecularDynamicsConfiguration& configuration,
                     const coupling::configurations::MaMiCoConfiguration<MDSIMULATIONFACTORY_DIMENSION>& mamicoConfiguration
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
                     ,
                     MPI_Comm localComm
#endif
                     )
      : coupling::interface::MDSimulation(), _lmp(new LAMMPS_NS::LAMMPS(0, NULL, localComm)), _configuration(configuration),
        _mamicoConfiguration(mamicoConfiguration), _tolerance(1.0e-8) {
  }

  virtual void writeCheckpoint(const std::string& filestem, const unsigned int& t) {
    std::stringstream command;
    command << "write_dump all atom restart_checkpoint.dump";
    _lmp->input->one(command.str().c_str());
  } // TODO

  // switch off coupling in simulation -> we assume the mamico fix to have ID=2
  virtual void switchOffCoupling() {
    // nop -> this only works if we need this for initial equilibration
  }

  // switch on coupling -> add mamico fix with ID=2; we further check
  // correctness of number of cells
  // -> afterwards, the MDSolverInterface of mamico will be initialised and
  // provided via MamicoInterfaceProvider.
  virtual void switchOnCoupling() {
    std::stringstream ss;
    unsigned int localNumberCells = 1;
    const tarch::la::Vector<MDSIMULATIONFACTORY_DIMENSION, unsigned int> numberProcs = _configuration.getMPIConfiguration().getNumberOfProcesses();
    const tarch::la::Vector<MDSIMULATIONFACTORY_DIMENSION, double> domainSize = _configuration.getDomainConfiguration().getGlobalDomainSize();
    const tarch::la::Vector<MDSIMULATIONFACTORY_DIMENSION, double> meshsize = _mamicoConfiguration.getMacroscopicCellConfiguration().getCouplingCellSize();
    // determine max. number of cells required on each process
    for (unsigned int d = 0; d < MDSIMULATIONFACTORY_DIMENSION; d++) {
      const unsigned int cells = (unsigned int)ceil(domainSize[d] / meshsize[d] / numberProcs[d]);
      if (!tarch::la::equals(cells * numberProcs[d] * meshsize[d], domainSize[d], _tolerance)) {
        std::cout << "ERROR switchOnCoupling(): cells and "
                     "domainsize/meshsize/numberProcs do not match!"
                  << std::endl;
        exit(EXIT_FAILURE);
      }
      // take cells+2 as we also need to allocate ghost layer
      localNumberCells = localNumberCells * (cells + 2);
    }

    ss << "fix 2 all mamico " << localNumberCells << " 123456 " << _configuration.getDomainConfiguration().getCutoffRadius();
    _lmp->input->one(ss.str().c_str());
  }

  // execute LAMMPS time steps
  virtual void simulateTimesteps(const unsigned int& numberTimesteps, const unsigned int& firstTimestep) {
    std::stringstream ss;
    ss << "run " << numberTimesteps << " post no";
    _lmp->input->one(ss.str().c_str());
  }

  virtual void sortMoleculesIntoCells() {
    LAMMPS_NS::MamicoLammpsMDSolverInterface<MDSIMULATIONFACTORY_DIMENSION>* interface =
        (LAMMPS_NS::MamicoLammpsMDSolverInterface<MDSIMULATIONFACTORY_DIMENSION>*)
            coupling::interface::MamicoInterfaceProvider<LAMMPS_NS::MamicoCell, MDSIMULATIONFACTORY_DIMENSION>::getInstance()
                .getMDSolverInterface();
    if (interface == NULL) {
      std::cout << "ERROR sortMoleculesIntoCells(): interface==NULL!" << std::endl;
      exit(EXIT_FAILURE);
    }
    coupling::services::MacroscopicCellService<MDSIMULATIONFACTORY_DIMENSION>* macroscopicCellService =
        coupling::interface::MamicoInterfaceProvider<LAMMPS_NS::MamicoCell, MDSIMULATIONFACTORY_DIMENSION>::getInstance().getMacroscopicCellService();
    if (macroscopicCellService == NULL) {
      std::cout << "ERROR sortMoleculesIntoCells(): macroscopicCellService==NULL!" << std::endl;
      exit(EXIT_FAILURE);
    }
    interface->updateAllCells(macroscopicCellService->getIndexConversion());
  }

  // set macroscopic cell service to MamicoInterfaceProvider
  virtual void setCouplingCellService(coupling::services::MacroscopicCellService<MDSIMULATIONFACTORY_DIMENSION>* macroscopicCellService) {
    coupling::interface::MamicoInterfaceProvider<LAMMPS_NS::MamicoCell, MDSIMULATIONFACTORY_DIMENSION>::getInstance().setCouplingCellService(
        macroscopicCellService);
  }

  // init LAMMPS simulation for LJ as provided in xml-config
  virtual void init() { initSingleSimulation(MPI_COMM_WORLD, 0); }

  virtual void init(const tarch::utils::MultiMDService<MDSIMULATIONFACTORY_DIMENSION>& multiMDService, unsigned int localMDSimulation) {
    initSingleSimulation(multiMDService.getLocalCommunicator(), localMDSimulation);
  }

  // nop required
  virtual void shutdown() {
    if (_lmp != NULL) {
      delete _lmp;
      _lmp = NULL;
    }
  }

private:
  void initSingleSimulation(MPI_Comm comm, unsigned int localMDSimulation) {
    int rank;
    MPI_Comm_rank(comm, &rank);
    std::stringstream ss;
    std::string inputFile;

    // use logging separately for each LAMMPS instance
    ss.str("");
    ss << "log log.lammps" << localMDSimulation;
    _lmp->input->one(ss.str().c_str());
    ss.str("");
    ss << "log log.lammps" << localMDSimulation << " append";
    _lmp->input->one(ss.str().c_str());

    // define lj atomic input

    _lmp->input->one("units lj");
    _lmp->input->one("atom_style atomic");
    // define number of processes
    ss.str("");
    tarch::la::Vector<3, unsigned int> procs(1);
    for (int d = 0; d < MDSIMULATIONFACTORY_DIMENSION; d++) {
      procs[d] = _configuration.getMPIConfiguration().getNumberOfProcesses()[d];
    }
    ss << "processors " << procs[0] << " " << procs[1] << " " << procs[2] << " map xyz";
    _lmp->input->one(ss.str().c_str());

    // define 2D if necessary
    if (MDSIMULATIONFACTORY_DIMENSION == 2) {
      _lmp->input->one("dimension 2");
    } else {
      _lmp->input->one("dimension 3");
    }

    // incorporate boundary conditions: currently, only reflecting and periodic
    // boundaries are supported; to be extended to 2D
    ss.str("");
    ss << "boundary";
    if (MDSIMULATIONFACTORY_DIMENSION == 3) {
      const tarch::la::Vector<MD_LINKED_CELL_NEIGHBOURS, simplemd::BoundaryType> bc = _configuration.getDomainConfiguration().getBoundary();
      int indices[3] = {12, 10, 4};
      for (int i = 0; i < MDSIMULATIONFACTORY_DIMENSION; i++) {
        if (bc[indices[i]] == simplemd::PERIODIC_BOUNDARY) {
          ss << " p";
        } else if (bc[indices[i]] == simplemd::REFLECTING_BOUNDARY) {
          ss << " f";
        } else {
          std::cout << "ERROR LammpsMDSimulation::initSingleSimulation(): only "
                       "reflecting and periodic boundaries supported!"
                    << std::endl;
          exit(EXIT_FAILURE);
        }
      }
    } else {
      std::cout << "ERROR LammpsMDSimulation::initSingleSimulation(): only 3D "
                   "supported at the moment!"
                << std::endl;
    }
    _lmp->input->one(ss.str().c_str());

    // create atoms based on input configuration; currently, only creation of
    // molecules on a lattice is supported -------------------
    tarch::la::Vector<3, double> domainOffset(0.0);
    for (int d = 0; d < MDSIMULATIONFACTORY_DIMENSION; d++) {
      domainOffset[d] = _configuration.getDomainConfiguration().getGlobalDomainOffset()[d];
    }
    tarch::la::Vector<3, double> domainSize(1.0);
    for (int d = 0; d < MDSIMULATIONFACTORY_DIMENSION; d++) {
      domainSize[d] = _configuration.getDomainConfiguration().getGlobalDomainSize()[d];
    }
    if (!_configuration.getDomainConfiguration().initFromCheckpoint() && !_configuration.getDomainConfiguration().initFromSequentialCheckpoint()) {
      tarch::la::Vector<3, unsigned int> moleculesPerDirection(0);
      for (int d = 0; d < MDSIMULATIONFACTORY_DIMENSION; d++) {
        moleculesPerDirection[d] = _configuration.getDomainConfiguration().getMoleculesPerDirection()[d];
      }
      double particleDensity = 1.0;
      for (int d = 0; d < MDSIMULATIONFACTORY_DIMENSION; d++) {
        particleDensity = particleDensity * _configuration.getDomainConfiguration().getMoleculesPerDirection()[d] / domainSize[d];
      }
    }
    // --------------------------------------------------------------------------------------------------------------------------------------------

    // define molecule and LJ parameters
    ss.str("");
    ss << "pair_style lj/cut " << _configuration.getDomainConfiguration().getCutoffRadius();
    _lmp->input->one(ss.str().c_str());
    _lmp->input->one("pair_modify shift yes");
    ss.str("");
    ss << "pair_coeff 1 1 " << _configuration.getMoleculeConfiguration().getEpsilon() << " ";
    ss << _configuration.getMoleculeConfiguration().getSigma() << " " << _configuration.getDomainConfiguration().getCutoffRadius();
    _lmp->input->one(ss.str().c_str());

    // set time step
    ss.str("");
    ss << "timestep " << _configuration.getSimulationConfiguration().getDt();
    _lmp->input->one(ss.str().c_str());

    // set some verlet-specific values; currently hard-coded
    _lmp->input->one("neighbor 0.5 bin");
    _lmp->input->one("neigh_modify delay 0 every 20 check no");

    // write output; though this is not vtk format, we use this config for that purpose
    if (_configuration.getVTKConfiguration().getWriteEveryTimestep() != 0) {
      ss.str("");
      ss << "dump myDump all custom " << _configuration.getVTKConfiguration().getWriteEveryTimestep() << " ";
      ss << _configuration.getVTKConfiguration().getFilename() << "_" << localMDSimulation << "_"
         << "*.dat id type x y z vx vy vz";
      _lmp->input->one(ss.str().c_str());
    }

    // define simulation region
    ss.str("");
    ss << "region box block " << domainOffset[0] << " " << (domainOffset[0] + domainSize[0]) << " " << domainOffset[1] << " "
       << (domainOffset[1] + domainSize[1]) << " ";
    ss << domainOffset[2] << " " << (domainOffset[2] + domainSize[2]) << " units box";
    _lmp->input->one(ss.str().c_str());

    if (_configuration.getDomainConfiguration().initFromCheckpoint() || _configuration.getDomainConfiguration().initFromSequentialCheckpoint()) {
      std::stringstream command;
      command << "read_restart " << _configuration.getDomainConfiguration().getCheckpointFilestem() << ".restart";
      _lmp->input->one(command.str().c_str());
    } else {
      _lmp->input->one("create_box 1 box");
      _lmp->input->one("create_atoms 1 box");
      // set mass and correct temperature
      ss.str("");
      ss << "mass 1 " << _configuration.getMoleculeConfiguration().getMass();
      _lmp->input->one(ss.str().c_str());
      // initialise velocities on molecules; for parallel simulations, we do
      // this according to global MD simulation number (we add 1 as 0 is not a
      // valid seed)
      ss.str("");
      ss << "velocity all create " << _configuration.getMoleculeConfiguration().getTemperature() << " " << localMDSimulation + 1 << " loop geom";
      _lmp->input->one(ss.str().c_str());
    }

    // incorporate reflecting boundaries, pt.2 (this needs to be done after the
    // box was created); to be extended to 2D
    if (MDSIMULATIONFACTORY_DIMENSION == 3) {
      const tarch::la::Vector<MD_LINKED_CELL_NEIGHBOURS, simplemd::BoundaryType> bc = _configuration.getDomainConfiguration().getBoundary();
      int indices[3] = {12, 10, 4};
      std::string faces[3] = {"x", "y", "z"};
      for (int i = 0; i < MDSIMULATIONFACTORY_DIMENSION; i++) {
        if (bc[indices[i]] == simplemd::REFLECTING_BOUNDARY) {
          ss.str("");
          ss << "fix reflection" << faces[i] << " all wall/reflect " << faces[i] << "lo EDGE " << faces[i] << "hi EDGE units box";
          _lmp->input->one(ss.str().c_str());
        }
      }
    }
    // --------------------------------------------------------------------------------------------------------------------------------------------

    // define molecule and LJ parameters
    ss.str("");
    ss << "pair_style lj/cut " << _configuration.getDomainConfiguration().getCutoffRadius();
    _lmp->input->one(ss.str().c_str());
    _lmp->input->one("pair_modify shift yes");
    ss.str("");
    ss << "pair_coeff 1 1 " << _configuration.getMoleculeConfiguration().getEpsilon() << " ";
    ss << _configuration.getMoleculeConfiguration().getSigma() << " " << _configuration.getDomainConfiguration().getCutoffRadius();
    _lmp->input->one(ss.str().c_str());

    // set time step
    ss.str("");
    ss << "timestep " << _configuration.getSimulationConfiguration().getDt();
    _lmp->input->one(ss.str().c_str());

    // set some verlet-specific values; currently hard-coded
    _lmp->input->one("neighbor 0.5 bin");
    _lmp->input->one("neigh_modify delay 0 every 1 check no");

    // write output; though this is not vtk format, we use this config for that
    // purpose
    if (_configuration.getVTKConfiguration().getWriteEveryTimestep() != 0) {
      ss.str("");
      ss << "dump myDump all custom " << _configuration.getVTKConfiguration().getWriteEveryTimestep() << " ";
      ss << _configuration.getVTKConfiguration().getFilename() << "_" << localMDSimulation << "_"
         << "*.dat id type x y z vx vy vz";
      _lmp->input->one(ss.str().c_str());
    }
    // write checkpointing
    if (_configuration.getCheckpointConfiguration().getWriteEveryTimestep() != 0) {
      ss.str("");
      ss << "restart " << _configuration.getCheckpointConfiguration().getWriteEveryTimestep() << " ";
      ss << _configuration.getCheckpointConfiguration().getFilename() << "_" << localMDSimulation << "_"
         << "1.restart " << _configuration.getCheckpointConfiguration().getFilename() << "_" << localMDSimulation << "_"
         << "2.restart";
      _lmp->input->one(ss.str().c_str());
    }
    // set nve ensemble
    _lmp->input->one("fix 1 all nve");
  }
};
} // namespace interface
} // namespace coupling