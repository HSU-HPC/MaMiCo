#pragma once

#include "USER-MAMICO/mamico_cell.h"
#include "USER-MAMICO/mamico_lammps_md_solver_interface.h"
#include "coupling/interface/MDSimulation.h"
#include "lammps/input.h"
#include "lammps/lammps.h"
#include <mpi.h>

namespace coupling {
namespace interface {
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
} // namespace interface
} // namespace coupling
