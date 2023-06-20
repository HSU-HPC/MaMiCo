// Copyright (C) 2016 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MDSIMULATIONFACTORY_H_
#define _MDSIMULATIONFACTORY_H_

#include <ctime>
#include <math.h>
#include <string>

#include "coupling/CouplingMDDefinitions.h"
#include "coupling/configurations/MaMiCoConfiguration.h"
#include "coupling/interface/MamicoInterfaceProvider.h"
#include "coupling/services/MacroscopicCellService.h"
#include "simplemd/configurations/MolecularDynamicsConfiguration.h"
#include "tarch/la/ScalarOperations.h"
#include "tarch/utils/MultiMDService.h"

// hacked: "currently, only 3D is supported!"
#define MDSIMULATIONFACTORY_DIMENSION 3

#if defined(SIMPLE_MD)
#include "coupling/interface/impl/SimpleMD/SimpleMDSolverInterface.h"
#include "coupling/solvers/CoupledMolecularDynamicsSimulation.h"
#include "simplemd/LinkedCell.h"
#define MY_LINKEDCELL simplemd::LinkedCell
#elif defined(LAMMPS_MD) || defined(LAMMPS_DPD)
#include "impl/LAMMPS/USER-MAMICO/mamico_cell.h"
#include "impl/LAMMPS/USER-MAMICO/mamico_lammps_md_solver_interface.h"
#include "lammps/input.h"
#include "lammps/lammps.h"
#include <mpi.h>
#define MY_LINKEDCELL LAMMPS_NS::MamicoCell
#elif defined(LS1_MARDYN)
#include "coupling/interface/impl/ls1/LS1MDSolverInterface.h"
#include "coupling/interface/impl/ls1/LS1RegionWrapper.h"
#include "coupling/interface/impl/ls1/LS1StaticCommData.h"
#include "ls1/src/Simulation.h"
#include "ls1/src/plugins/MamicoCoupling.h"
#define MY_LINKEDCELL ls1::LS1RegionWrapper
#endif

/** interface for different MD solvers.
 *  @author Philipp Neumann
 */
namespace coupling {
namespace interface {

/**
 *	@brief generic interface class for different microscopic (MD) solvers.
 *  @author Philipp Neumann
 */
class MDSimulation {
public:
  /** Destructor */
  virtual ~MDSimulation() {}

  /** switches coupling off*/
  virtual void switchOffCoupling() = 0;

  /** switches coupling on*/
  virtual void switchOnCoupling() = 0;

  /** simulates numberTimesteps time steps and starts at time step no.
   *firstTimestep
   *	@param numberTimesteps
   *	@param firstTimestep
   */
  virtual void simulateTimesteps(const unsigned int& numberTimesteps, const unsigned int& firstTimestep) = 0;

  // simulates a single time step
  // virtual void simulateTimestep(const unsigned int &thisTimestep ){const
  // unsigned int steps=1; simulateTimesteps(thisTimestep,steps);} TODO BUG

  /** sortMoleculesIntoCells*/
  virtual void sortMoleculesIntoCells() = 0;

  /** setMacroscopicCellService
   *	@param macroscopicCellService
   */
  virtual void setMacroscopicCellService(coupling::services::MacroscopicCellService<MDSIMULATIONFACTORY_DIMENSION>* macroscopicCellService) = 0;

  /** initialises the _molecularDynamicsSimulation solver
   *	@sa simplemd::MolecularDynamicsSimulation::initServices()
   *	@todo Philipp ??
   */
  virtual void init() = 0;

  /** initialises the _molecularDynamicsSimulation solver
   *	@param multiMDService
   *	@param localMDSimulation
   *	@sa simplemd::MolecularDynamicsSimulation::initServices(const
   *tarch::utils::MultiMDService<MD_DIM>& multiMDService,unsigned int
   *localMDSimulation)
   *	@todo Philipp ??
   */
  virtual void init(const tarch::utils::MultiMDService<MDSIMULATIONFACTORY_DIMENSION>& multiMDService, unsigned int localMDSimulation) = 0;

  /** shuts down the MD simulation*/
  virtual void shutdown() = 0;

  /** Saves the simulation result as check point in the file filestem
   *	@param filestem
   *	@param t
   */
  virtual void writeCheckpoint(const std::string& filestem, const unsigned int& t) = 0;
};

/** define MD simulation from default MD code */
#if defined(SIMPLE_MD)
/** defines MD simulation from default MD code. Derived from the class
 *MDSimulation.
 *	@brief defines MD simulation from default MD code.
 *  @author Philipp Neumann
 */
class SimpleMDSimulation : public coupling::interface::MDSimulation {
public:
  SimpleMDSimulation(const simplemd::configurations::MolecularDynamicsConfiguration& configuration)
      : coupling::interface::MDSimulation(), _molecularDynamicsSimulation(configuration) {}
  virtual ~SimpleMDSimulation() {}

  virtual void switchOffCoupling() { _molecularDynamicsSimulation.switchOffCoupling(); }
  virtual void switchOnCoupling() { _molecularDynamicsSimulation.switchOnCoupling(); }
  virtual void simulateTimesteps(const unsigned int& numberTimesteps, const unsigned int& firstTimestep) {
    for (unsigned int t = firstTimestep; t < firstTimestep + numberTimesteps; t++) {
      _molecularDynamicsSimulation.simulateOneCouplingTimestep(t);
    }
  }
  virtual void sortMoleculesIntoCells() {
    // nop required, since the linked cells are very tightly linked to mamico
  }

  virtual void setMacroscopicCellService(coupling::services::MacroscopicCellService<MDSIMULATIONFACTORY_DIMENSION>* macroscopicCellService) {
    _molecularDynamicsSimulation.setMacroscopicCellService(macroscopicCellService);
    // set the cell service also in singleton of mamico interface provider ->
    // typically not required in coupling, but makes the simulation state more
    // consistent compared to using LAMMPS
    coupling::interface::MamicoInterfaceProvider<simplemd::LinkedCell, MDSIMULATIONFACTORY_DIMENSION>::getInstance().setMacroscopicCellService(
        macroscopicCellService);
  }
  virtual void init() { _molecularDynamicsSimulation.initServices(); }
  virtual void init(const tarch::utils::MultiMDService<MDSIMULATIONFACTORY_DIMENSION>& multiMDService, unsigned int localMDSimulation) {
    _molecularDynamicsSimulation.initServices(multiMDService, localMDSimulation);
  }
  virtual void shutdown() { _molecularDynamicsSimulation.shutdownServices(); }

  virtual void writeCheckpoint(const std::string& filestem, const unsigned int& t) {
    getMoleculeService().writeCheckPoint(getParallelTopologyService(), filestem, t);
  }

  // function particularly needed to init MD solver interface -> should only be
  // called from factory
  simplemd::BoundaryTreatment& getBoundaryTreatment() { return _molecularDynamicsSimulation.getBoundaryTreatment(); }
  simplemd::services::ParallelTopologyService& getParallelTopologyService() { return _molecularDynamicsSimulation.getParallelTopologyService(); }
  simplemd::services::MoleculeService& getMoleculeService() { return _molecularDynamicsSimulation.getMoleculeService(); }
  simplemd::services::LinkedCellService& getLinkedCellService() { return _molecularDynamicsSimulation.getLinkedCellService(); }
  const simplemd::services::MolecularPropertiesService& getMolecularPropertiesService() { return _molecularDynamicsSimulation.getMolecularPropertiesService(); }

private:
  coupling::solvers::CoupledMolecularDynamicsSimulation _molecularDynamicsSimulation;
};
#endif

/** md simulation based on LAMMPS, encapsulated in this simulation interface */
#if defined(LAMMPS_MD)
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
    const tarch::la::Vector<MDSIMULATIONFACTORY_DIMENSION, double> meshsize = _mamicoConfiguration.getMacroscopicCellConfiguration().getMacroscopicCellSize();
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
  virtual void setMacroscopicCellService(coupling::services::MacroscopicCellService<MDSIMULATIONFACTORY_DIMENSION>* macroscopicCellService) {
    coupling::interface::MamicoInterfaceProvider<LAMMPS_NS::MamicoCell, MDSIMULATIONFACTORY_DIMENSION>::getInstance().setMacroscopicCellService(
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
#endif

#if defined(LAMMPS_DPD)
class LammpsDPDSimulation : public coupling::interface::MDSimulation {
private:
  LAMMPS_NS::LAMMPS* _lmp;
  const simplemd::configurations::MolecularDynamicsConfiguration& _configuration;
  const coupling::configurations::MaMiCoConfiguration<MDSIMULATIONFACTORY_DIMENSION>& _mamicoConfiguration;
  const double _tolerance; // tolerance for checks of configs

public:
  LammpsDPDSimulation(const simplemd::configurations::MolecularDynamicsConfiguration& configuration,
                      const coupling::configurations::MaMiCoConfiguration<MDSIMULATIONFACTORY_DIMENSION>& mamicoConfiguration
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
                      ,
                      MPI_Comm localComm
#endif
                      )
      : coupling::interface::MDSimulation(), _lmp(new LAMMPS_NS::LAMMPS(0, NULL, localComm)), _configuration(configuration),
        _mamicoConfiguration(mamicoConfiguration), _tolerance(1.0e-8) {
  }

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
    const tarch::la::Vector<MDSIMULATIONFACTORY_DIMENSION, double> meshsize = _mamicoConfiguration.getMacroscopicCellConfiguration().getMacroscopicCellSize();
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
  virtual void setMacroscopicCellService(coupling::services::MacroscopicCellService<MDSIMULATIONFACTORY_DIMENSION>* macroscopicCellService) {
    coupling::interface::MamicoInterfaceProvider<LAMMPS_NS::MamicoCell, MDSIMULATIONFACTORY_DIMENSION>::getInstance().setMacroscopicCellService(
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
    // particle insertion is currently only supported for single-centered LJ
    // particles
    if (_mamicoConfiguration.getParticleInsertionConfiguration().getParticleInsertionType() !=
        coupling::configurations::ParticleInsertionConfiguration::NO_INSERTION) {
      std::cout << "ERROR coupling::LammpsDPDSImulation::initSingleSimulation(): "
                   "Particle insertion currently not supported!"
                << std::endl;
      exit(EXIT_FAILURE);
    }

    // set logs for each LAMMPS instance separately
    ss.str("");
    ss << "log log.lammps" << localMDSimulation;
    _lmp->input->one(ss.str().c_str());
    ss.str("");
    ss << "log log.lammps" << localMDSimulation << " append";
    _lmp->input->one(ss.str().c_str());
    _lmp->input->one("variable ndim equal 3");
    ss.str("");
    ss << "variable kB equal " << _configuration.getDomainConfiguration().getKB();
    _lmp->input->one(ss.str().c_str());
    ss.str("");
    ss << "variable kBT equal " << _configuration.getMoleculeConfiguration().getTemperature() * _configuration.getDomainConfiguration().getKB();
    _lmp->input->one(ss.str().c_str());
    ss.str("");
    ss << "variable T equal " << _configuration.getMoleculeConfiguration().getTemperature();
    _lmp->input->one(ss.str().c_str());
    ss.str("");
    ss << "variable rc equal " << _configuration.getDomainConfiguration().getCutoffRadius();
    _lmp->input->one(ss.str().c_str());
    _lmp->input->one("units lj");
    _lmp->input->one("dimension ${ndim}");

    // set up boundary treatment, Pt1: for "periodic", we just set the
    // "p"-option and "f" otherwise; for the latter, we further have to
    // introduce the respective walls, see boundary treatment, Pt2, further down
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
          std::cout << "ERROR LammpsDPDSimulation::initSingleSimulation(): "
                       "only reflecting and periodic boundaries supported!"
                    << std::endl;
          exit(EXIT_FAILURE);
        }
      }
    } else {
      std::cout << "ERROR LammpsDPDSimulation::initSingleSimulation(): only 3D "
                   "supported at the moment!"
                << std::endl;
    }
    _lmp->input->one(ss.str().c_str());

    _lmp->input->one("atom_style  atomic");
    // define number of processes with xyz-topology
    ss.str("");
    ss << "processors " << _configuration.getMPIConfiguration().getNumberOfProcesses()[0] << " "
       << _configuration.getMPIConfiguration().getNumberOfProcesses()[1] << " ";
    ss << _configuration.getMPIConfiguration().getNumberOfProcesses()[2] << " map xyz";
    _lmp->input->one(ss.str().c_str());
    ss.str("");
    ss << "variable    seed1 equal " << localMDSimulation + rank + 1;
    _lmp->input->one(ss.str().c_str()); // currently fixed random seed
    ss.str("");
    ss << "variable    seed2 equal " << localMDSimulation + rank + 10;
    _lmp->input->one(ss.str().c_str()); // new seed for velocity generation;
                                        // currently deduced from _randSeed+10
    _lmp->input->one("print \"seed1: ${seed1}\"");
    _lmp->input->one("print \"seed2: ${seed2}\"");
    _lmp->input->one("pair_style dpd/general ${kBT} ${rc} ${seed1}"); // DPD scheme with
                                                                      // forcing by Xin

    const tarch::la::Vector<MDSIMULATIONFACTORY_DIMENSION, double> domainSize(_configuration.getDomainConfiguration().getGlobalDomainSize());
    const tarch::la::Vector<MDSIMULATIONFACTORY_DIMENSION, double> domainOffset(_configuration.getDomainConfiguration().getGlobalDomainOffset());
    ss.str("");
    ss << "variable  Lx equal " << domainSize[0];
    _lmp->input->one(ss.str().c_str());
    ss.str("");
    ss << "variable  Ly equal " << domainSize[1];
    _lmp->input->one(ss.str().c_str());
    ss.str("");
    ss << "variable  Lz equal " << domainSize[2];
    _lmp->input->one(ss.str().c_str());
    // The following command computes the average density based on
    // molecules-per-direction. However, particles will be randomly distributed
    // - they are not arranged on a grid
    const tarch::la::Vector<MDSIMULATIONFACTORY_DIMENSION, unsigned int> molecules(_configuration.getDomainConfiguration().getMoleculesPerDirection());
    unsigned int numberMolecules = 1;
    for (unsigned int d = 0; d < MDSIMULATIONFACTORY_DIMENSION; d++) {
      numberMolecules = numberMolecules * molecules[d];
    }
    ss.str("");
    ss << "variable Npart equal " << numberMolecules;
    _lmp->input->one(ss.str().c_str());
    _lmp->input->one("variable  Ntype equal 1");

    ss.str("");
    ss << "region region_dpd block " << domainOffset[0] << " " << domainOffset[0] + domainSize[0] << " " << domainOffset[1] << " "
       << domainOffset[1] + domainSize[1] << " ";
    ss << domainOffset[2] << " " << domainOffset[2] + domainSize[2] << " units box";
    _lmp->input->one(ss.str().c_str());
    _lmp->input->one("create_box    ${Ntype} region_dpd");
    _lmp->input->one("create_atoms 1 random ${Npart} 12345 NULL"); // new randomized
                                                                   // initialisation

    // boundary treatment, Pt2: introduce reflecting walls if required (fixes
    // need to be defined after box definition)
    if (MDSIMULATIONFACTORY_DIMENSION == 3) {
      const tarch::la::Vector<MD_LINKED_CELL_NEIGHBOURS, simplemd::BoundaryType> bc = _configuration.getDomainConfiguration().getBoundary();
      int indices[3] = {12, 10, 4};
      if (bc[indices[0]] == simplemd::REFLECTING_BOUNDARY) {
        _lmp->input->one("fix reflectionX all wall/reflect xlo EDGE xhi EDGE");
      }
      if (bc[indices[1]] == simplemd::REFLECTING_BOUNDARY) {
        _lmp->input->one("fix reflectionY all wall/reflect ylo EDGE yhi EDGE");
      }
      if (bc[indices[2]] == simplemd::REFLECTING_BOUNDARY) {
        _lmp->input->one("fix reflectionZ all wall/reflect zlo EDGE zhi EDGE");
      }
    }

    ss.str("");
    ss << "variable m equal " << _configuration.getMoleculeConfiguration().getMass();
    _lmp->input->one(ss.str().c_str());

    ss.str("");
    ss << "variable alpha equal " << _configuration.getMoleculeConfiguration().getEpsilon(); // configurable via epsilon parameter in xml config
    _lmp->input->one(ss.str().c_str());
    ss.str("");
    ss << "variable sigma equal " << _configuration.getMoleculeConfiguration().getSigma();
    _lmp->input->one(ss.str().c_str());
    _lmp->input->one("variable  gamma  equal ${sigma}^2/(2*${kBT})");

    _lmp->input->one("variable power equal 0.25");

    ss.str("");
    ss << "variable Ntime equal " << _configuration.getSimulationConfiguration().getNumberOfTimesteps();
    _lmp->input->one(ss.str().c_str());
    ss.str("");
    ss << "variable Nrepeat equal " << (int)(0.8 * _configuration.getSimulationConfiguration().getNumberOfTimesteps());
    _lmp->input->one(ss.str().c_str());

    _lmp->input->one("pair_coeff      1 1 ${alpha} ${gamma} ${power} ${rc}");

    _lmp->input->one("neighbor    0.5 bin");
    _lmp->input->one("neigh_modify    delay 0 every 1 check no");

    _lmp->input->one("communicate single vel yes");
    _lmp->input->one("group               extra_atoms id > ${Npart}"); // remove particles to
                                                                       // reach correct
                                                                       // density
    _lmp->input->one("delete_atoms        group extra_atoms");

    ss.str("");
    ss << "timestep " << _configuration.getSimulationConfiguration().getDt();
    _lmp->input->one(ss.str().c_str());

    _lmp->input->one("mass                1 ${m}");
    _lmp->input->one("velocity    all create ${kBT} ${seed2} loop geom");

    _lmp->input->one("thermo 1");
    _lmp->input->one("thermo_style custom step time ke temp press");
    _lmp->input->one("thermo_modify       flush yes");

    // wall fix bottom
    if (_mamicoConfiguration.getBoundaryForceConfiguration().getBoundary()[4]) {
      if (_configuration.getDomainConfiguration().getBoundary()[4] != simplemd::REFLECTING_BOUNDARY) {
        std::cout << "ERROR "
                     "coupling::MDSimulationFactory::LammpsDPDSimulation::"
                     "initSingleSimulation(): bottom boundary is not reflecting!"
                  << std::endl;
        exit(EXIT_FAILURE);
      }
      ss.str("");
      ss << "fix lower_wall all dpd/boundary/force ${rc} 100 dpd_plane_fc.dat "
            "dpd_plane_fd1.dat dpd_plane_fd2.dat 2 0 "
         << domainOffset[2] << " 0.0";
      _lmp->input->one(ss.str().c_str());
    }
    // wall fix top
    if (_mamicoConfiguration.getBoundaryForceConfiguration().getBoundary()[5]) {
      if (_configuration.getDomainConfiguration().getBoundary()[21] != simplemd::REFLECTING_BOUNDARY) {
        std::cout << "ERROR "
                     "coupling::MDSimulationFactory::LammpsDPDSimulation::"
                     "initSingleSimulation(): top boundary is not reflecting!"
                  << std::endl;
        exit(EXIT_FAILURE);
      }
      ss.str("");
      ss << "fix upper_wall all dpd/interface/force ${rc} 100 dpd_plane_fc.dat "
            "dpd_plane_fd1.dat dpd_plane_fd2.dat 2 1 "
         << domainOffset[2] + domainSize[2] << " 0.0";
      _lmp->input->one(ss.str().c_str());
    }

    _lmp->input->one("fix 1 all nve");
  }
};
#endif

#if defined(LS1_MARDYN)
class LS1MDSimulation : public coupling::interface::MDSimulation {
private:
  const simplemd::configurations::MolecularDynamicsConfiguration& _configuration;
  Simulation* simulation; //cannot name this _simulation, a global preprocessor marco with the name _simulation expands to *global_simulation
  MamicoCoupling* ls1MamicoPlugin; //the plugin is only initialized after the simulation object reads xml, so cannot use it before that point
  bool internalCouplingState;
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
  MPI_Comm comm;
#endif

public:
  LS1MDSimulation(const simplemd::configurations::MolecularDynamicsConfiguration& configuration
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
                    , MPI_Comm localComm
#endif
  ) : coupling::interface::MDSimulation(), _configuration(configuration) {

#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
    comm = localComm;
    LS1StaticCommData::getInstance().setLocalCommunicator(comm); //needs to be done before creating the simulation object
    const tarch::la::Vector<3, unsigned int> numberProcs = _configuration.getMPIConfiguration().getNumberOfProcesses();
    LS1StaticCommData::getInstance().setDomainGridDecompAtDim(0, numberProcs[0]);
    LS1StaticCommData::getInstance().setDomainGridDecompAtDim(1, numberProcs[1]);
    LS1StaticCommData::getInstance().setDomainGridDecompAtDim(2, numberProcs[2]);
#endif

    simulation = new Simulation();
    global_simulation = simulation;
    simulation->disableFinalCheckpoint();
    internalCouplingState = false;
    ls1MamicoPlugin = nullptr;
  }
  virtual ~LS1MDSimulation() {
    if(simulation != nullptr) {
      delete simulation;
      simulation = nullptr;
    }
  }
  /** switches coupling on/off*/
  virtual void switchOffCoupling() override { 
    //coupling::interface::LS1MamicoCouplingSwitch::getInstance().setCouplingStateOff();
    internalCouplingState = false;
    if(ls1MamicoPlugin != nullptr)
      ls1MamicoPlugin->switchOffCoupling();
  }
  virtual void switchOnCoupling() override { 
    //coupling::interface::LS1MamicoCouplingSwitch::getInstance().setCouplingStateOn();
    internalCouplingState = true;
    if(ls1MamicoPlugin != nullptr)
      ls1MamicoPlugin->switchOnCoupling();
  }

  /** simulates numberTimesteps time steps and starts at time step no.
   * firstTimestep*/
  virtual void simulateTimesteps(const unsigned int& numberTimesteps, const unsigned int& firstTimestep) override {
    global_simulation = simulation;
    for (unsigned int i = 0; i < numberTimesteps; i++) {
      simulation->simulateOneTimestep();
    }
  }
  /** simulates a single time step*/
  // virtual void simulateTimestep(const unsigned int &thisTimestep ){const
  // unsigned int steps=1; simulateTimesteps(thisTimestep,steps);} TODO BUG
  virtual void sortMoleculesIntoCells() override {}

  virtual void setMacroscopicCellService(coupling::services::MacroscopicCellService<MDSIMULATIONFACTORY_DIMENSION>* macroscopicCellService) override {
    //coupling::interface::MamicoInterfaceProvider<ls1::LS1RegionWrapper, MDSIMULATIONFACTORY_DIMENSION>::getInstance().setMacroscopicCellService(
    //    macroscopicCellService);
    global_simulation = simulation;
    PluginBase* searchedPlugin = simulation->getPlugin("MamicoCoupling");
    if(searchedPlugin == nullptr)
    {
      std::cout << "ERROR: MaMiCo plugin not found!" << std::endl;
      exit(EXIT_FAILURE);
    }
    ls1MamicoPlugin = dynamic_cast<MamicoCoupling*>(searchedPlugin);
    if(ls1MamicoPlugin != nullptr) {
      ls1MamicoPlugin->setMamicoMacroscopicCellService(macroscopicCellService);
    }
    else {
      std::cout << "ERROR: Cast to Mamico plugin unsuccessful!" << std::endl;
      exit(EXIT_FAILURE);
    }
    //since this is the first time the plugin is accessed, set whatever preexisting coupling variable we had here for the first time
    if(internalCouplingState)
      ls1MamicoPlugin->switchOnCoupling();
    else
      ls1MamicoPlugin->switchOffCoupling();
  }
  virtual void init() override {
    global_simulation = simulation;
    // parse file
    const std::string filename = coupling::interface::LS1StaticCommData::getInstance().getConfigFilename();
    simulation->readConfigFile(filename);
    // after this point the mamico plugin exists and is accessible
    simulation->prepare_start();
    simulation->preSimLoopSteps();
  }
  virtual void init(const tarch::utils::MultiMDService<MDSIMULATIONFACTORY_DIMENSION>& multiMDService, unsigned int localMDSimulation) override { init(); }
  virtual void shutdown() override {
    global_simulation = simulation;
    simulation->markSimAsDone();
    simulation->postSimLoopSteps();
    simulation->finalize();
  }
  virtual void writeCheckpoint(const std::string& filestem, const unsigned int& t) override {
    // configure through ls1 config file, using plugins
  }
};
#endif

/** factory to produced md simulation, md solver interface (for mamico) and the
 *macroscopic cell service using singleton pattern.
 *	@brief factory to produced md simulation, md solver interface (for
 *mamico) and the macroscopic cell service
 *  @author Philipp Neumann
 */
class SimulationAndInterfaceFactory {
public:
  /** @returns the SimulationAndInterfaceFactory object
   *	@note singleton pattern
   */
  static SimulationAndInterfaceFactory& getInstance() {
    static SimulationAndInterfaceFactory singleton;
    return singleton;
  }

  /** returns a pointer to the md simulation.
   *  @param configuration
   *  @param mamicoConfiguration
   *  @param localComm
   *	@remark This will create a new simulation which needs to be deleted at
   *the end
   */
  coupling::interface::MDSimulation* getMDSimulation(const simplemd::configurations::MolecularDynamicsConfiguration& configuration,
                                                     const coupling::configurations::MaMiCoConfiguration<MDSIMULATIONFACTORY_DIMENSION>& mamicoConfiguration
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
                                                     ,
                                                     MPI_Comm localComm
#endif
  ) {
#if defined(SIMPLE_MD)
    return new coupling::interface::SimpleMDSimulation(configuration);
#elif defined(LAMMPS_MD)
    return new coupling::interface::LammpsMDSimulation(configuration, mamicoConfiguration
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
                                                       ,
                                                       localComm
#endif
    );
#elif defined(LAMMPS_DPD)
    return new coupling::interface::LammpsDPDSimulation(configuration, mamicoConfiguration
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
                                                        ,
                                                        localComm
#endif
    );
#elif defined(LS1_MARDYN)
    return new coupling::interface::LS1MDSimulation(configuration
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
                                                    ,
                                                    localComm
#endif
    );
#else
    std::cout << "ERROR MDSimulationFactory::getMDSimulation(): Unknown MD "
                 "simulation!"
              << std::endl;
    exit(EXIT_FAILURE);
    return NULL;
#endif
  }

  /** returns the MD solver interface. This method should be called AFTER
   * initialising the MD simulation AND AFTER the equilibration of MD. We thus
   * expect that getMDSimulationInterface() and switchOnCoupling() of the
   * MDSimulationInterface() have been called before.
   *  @param configuration
   *  @param mamicoConfiguration
   *  @param mdSimulation
   */
  coupling::interface::MDSolverInterface<MY_LINKEDCELL, MDSIMULATIONFACTORY_DIMENSION>*
  getMDSolverInterface(const simplemd::configurations::MolecularDynamicsConfiguration& configuration,
                       const coupling::configurations::MaMiCoConfiguration<MDSIMULATIONFACTORY_DIMENSION>& mamicoConfiguration,
                       coupling::interface::MDSimulation* mdSimulation) {
    coupling::interface::MDSolverInterface<MY_LINKEDCELL, MDSIMULATIONFACTORY_DIMENSION>* mdSolverInterface = NULL;
#if defined(SIMPLE_MD)
    // for the simple MD code, we create a new md solver interface and also add
    // it to the mamico interface provider (the latter is not really required,
    // but makes the simulation state more consistent in the overall simulation)
    coupling::interface::SimpleMDSimulation* simpleMDSimulation = (coupling::interface::SimpleMDSimulation*)mdSimulation;
    if (simpleMDSimulation == NULL) {
      std::cout << "ERROR MDSimulationFactory::getMDSolverInterface(): Could "
                   "not cast to SimpleMDSimulation!"
                << std::endl;
      exit(EXIT_FAILURE);
    }
    mdSolverInterface = new coupling::interface::SimpleMDSolverInterface(
        mamicoConfiguration.getMacroscopicCellConfiguration().getNumberLinkedCellsPerMacroscopicCell(), simpleMDSimulation->getBoundaryTreatment(),
        simpleMDSimulation->getParallelTopologyService(), simpleMDSimulation->getMoleculeService(), simpleMDSimulation->getLinkedCellService(),
        simpleMDSimulation->getMolecularPropertiesService(), (simpleMDSimulation->getParallelTopologyService()).getLocalBoundaryInformation(),
        configuration.getSimulationConfiguration().getDt());
    coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL, MDSIMULATIONFACTORY_DIMENSION>::getInstance().setMDSolverInterface(mdSolverInterface);
#elif defined(LAMMPS_MD)
    // as switchOnCoupling() should have been called before this method, the
    // fix-mamico should be already initialised. hence the MDSolverInterface of
    // the mamico interface provider should be initialised at this stage
    mdSolverInterface = coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL, MDSIMULATIONFACTORY_DIMENSION>::getInstance().getMDSolverInterface();
#elif defined(LAMMPS_DPD)
    mdSolverInterface = coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL, MDSIMULATIONFACTORY_DIMENSION>::getInstance().getMDSolverInterface();
#elif defined(LS1_MARDYN)
    mdSolverInterface = new coupling::interface::LS1MDSolverInterface();
    coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL, MDSIMULATIONFACTORY_DIMENSION>::getInstance().setMDSolverInterface(mdSolverInterface);
#endif

    if (mdSolverInterface == NULL) {
      std::cout << "ERROR MDSimulationFactory::getMDSolverInterface(): "
                   "mdSolverInterface==NULL!"
                << std::endl;
      exit(EXIT_FAILURE);
    }
    return mdSolverInterface;
  }

  /** shuts down the MD solver interfaces and deletes the interface if required
     (depending on the respective MD simulation) */
  void
  shutdownMDSolverInterface() {
#if defined(SIMPLE_MD)
    coupling::interface::MDSolverInterface<MY_LINKEDCELL, MDSIMULATIONFACTORY_DIMENSION>* interface =
        coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL, MDSIMULATIONFACTORY_DIMENSION>::getInstance().getMDSolverInterface();
    if (interface != NULL) {
      delete interface;
    }
    coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL, MDSIMULATIONFACTORY_DIMENSION>::getInstance().setMDSolverInterface(NULL);
#elif defined(LAMMPS_MD)
// nop, since the fix-mamico takes care of deleting the MD solver interface
#elif defined(LAMMPS_DPD)
// nop
#elif defined(LS1_MARDYN)
    coupling::interface::MDSolverInterface<MY_LINKEDCELL, MDSIMULATIONFACTORY_DIMENSION>* interface =
        coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL, MDSIMULATIONFACTORY_DIMENSION>::getInstance().getMDSolverInterface();
    if (interface != NULL) {
      delete interface;
    }
    coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL, MDSIMULATIONFACTORY_DIMENSION>::getInstance().setMDSolverInterface(NULL);
#endif
  }

private:
  /** Private constructor
   *  @note singelton pattern
   */
  SimulationAndInterfaceFactory() {}
  /** Private destructor
   *  @note singelton pattern
   */
  ~SimulationAndInterfaceFactory() {}
};

} // namespace interface
} // namespace coupling
#endif // _MDSIMULATIONFACTORY_H_
