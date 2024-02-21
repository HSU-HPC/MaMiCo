// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef MAMICO_LMP_LAMMPS_MD_SOLVER_INTERFACE_H
#define MAMICO_LMP_LAMMPS_MD_SOLVER_INTERFACE_H

#include "coupling/CouplingMDDefinitions.h"
#include "coupling/interface/MDSolverInterface.h"
#include "coupling/interface/MamicoInterfaceProvider.h"
#include "tarch/la/ScalarOperations.h"

#include "mamico_cell.h"
#include "mamico_lammps_molecule_iterator.h"

#include "lammps/atom_vec.h"
#include "lammps/comm.h"
#include "lammps/compute_pe_atom.h"
#include "lammps/domain.h"
#include "lammps/error.h"
#include "lammps/fix.h"
#include "lammps/force.h"
#include "lammps/input.h"
#include "lammps/lammps.h"
#include "lammps/modify.h"
#include "lammps/neighbor.h"
#include "lammps/random_park.h"
#include "lammps/timer.h"
#include "lammps/update.h"
#include "lammps/verlet.h"
#include "sorting.h"
// the following three lines were included since we had issues when compiling
// with walberla -> compiler flag stcxx and INT64_MAX do not fit
// -> use LAMMPS_SMALLSMALL or so as compiler flag to build library
#define __STDC_LIMIT_MACROS
#include <stdint.h>
// #define MAXBIGINT INT_MAX

#include <iomanip>
#include <iostream>
#include <sstream>

namespace LAMMPS_NS {

/** current restrictions on the interface implementation:
 *  - domain in 2D must have z=0.0
 *  - time integration must follow the verlet-scheme
 *  - currently, only simple LJ-atoms are supported (at least in particle
 * insertion)
 *  - parallel domain decomposition must be set to xyz
 *  - only a single molecule type is supported; this type must be set to 1.
 *  @author Philipp Neumann
 */
template <unsigned int dim> class MamicoLammpsMDSolverInterface : public coupling::interface::MDSolverInterface<LAMMPS_NS::MamicoCell, dim> {
public:
  MamicoLammpsMDSolverInterface(LAMMPS* lmp, int seed, int numberCells, double cutoffRadius)
      : coupling::interface::MDSolverInterface<LAMMPS_NS::MamicoCell, dim>(), _lmp(lmp), _sorting(numberCells, lmp),
        _cutOffRadiusSquared(cutoffRadius * cutoffRadius), _cutOffRadius(cutoffRadius), _atomType(1), _precision(10), _tolerance(1.0e-8),
        _random(new RanPark(lmp, seed)) {
    if (_random == NULL)
      _lmp->error->all(FLERR, "Could not allocated new RanPark!");
  }
  virtual ~MamicoLammpsMDSolverInterface() {
    delete _random;
    _random = NULL;
    _lmp = NULL;
  }

  /** returns a particular linked cell inside a coupling cell.
   *  The coupling cells are currently located on the same process as the
   * respective linked cells. However, several linked cells may be part of a
   * coupling cell. The coupling cells also contain a ghost layer which
   * surrounds each local domain; the very first coupling cell inside the
   * global MD domain (or local MD domain) is thus given by coordinates (1,1,1)
   * (or (1,1) in 2D, respectively). The index linkedCellInCouplingCell
   * corresponds to the coordinates of the linked cell inside the given
   * coupling cell. These coordinates thus lie in a range
   * (0,linkedCellsPerCouplingCell-1).
   */
  virtual LAMMPS_NS::MamicoCell& getLinkedCell(const tarch::la::Vector<dim, unsigned int>& couplingCellIndex,
                                               const tarch::la::Vector<dim, unsigned int>& linkedCellInCouplingCell,
                                               const tarch::la::Vector<dim, unsigned int>& linkedCellsPerCouplingCell,
                                               const coupling::IndexConversion<dim>& indexConversion) {
    // we currently embed one cell per coupling cell. It is hence sufficient
    // to use the same enumeration as used for the coupling cells
    const unsigned int index = indexConversion.getLocalCellIndex(couplingCellIndex);
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
    std::cout << "Return linked cell at macro-cell " << couplingCellIndex << " no. " << linkedCellInCouplingCell << ", assuming " << linkedCellsPerCouplingCell
              << " linked cells per cell: " << index << std::endl;
#endif
    return _sorting.getMamicoCell(index);
  }

  /** returns the global size of the box-shaped MD domain */
  virtual tarch::la::Vector<dim, double> getGlobalMDDomainSize() const {
    if (_lmp->domain->triclinic == 1) {
      _lmp->error->all(FLERR, "Only orthogonal box shape domain is supported!");
    }
    tarch::la::Vector<dim, double> domainsize(0.0);
    for (unsigned int d = 0; d < dim; d++) {
      domainsize[d] = _lmp->domain->boxhi[d] - _lmp->domain->boxlo[d];
    }
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
    std::cout << "Domain size of MD domain: " << domainsize << std::endl;
#endif
    return domainsize;
  }

  /** returns the offset (i.e. lower,left corner) of MD domain */
  virtual tarch::la::Vector<dim, double> getGlobalMDDomainOffset() const {
    if (_lmp->domain->triclinic == 1) {
      _lmp->error->all(FLERR, "Only orthogonal box shape domain is supported!");
    }
    tarch::la::Vector<dim, double> domainoffset(0.0);
    for (unsigned int d = 0; d < dim; d++) {
      domainoffset[d] = _lmp->domain->boxlo[d];
    }
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
    std::cout << "Offset of MD domain: " << domainoffset << std::endl;
#endif
    return domainoffset;
  }

  /** returns the mass of a single fluid molecule */
  virtual double getMoleculeMass() const {
// currently, return mass of atom type no "1". We currently assume thus, that
// all our molecules are of type 1.
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
    std::cout << "Mass of an atom/molecule: " << _lmp->atom->mass[_atomType] << std::endl;
#endif
    return _lmp->atom->mass[_atomType];
  }

  /** returns Boltzmann's constant. For LJ fluids, this should always be 1.0,
   * see doc/units.html. */
  virtual double getKB() const {
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
    std::cout << "Boltzmann constant: " << _lmp->force->boltz << std::endl;
#endif
    return _lmp->force->boltz;
  }

  /** returns the sigma parameter of the LJ potential. We always return 1.0, see
   * also doc/units.html for LJ style. */
  virtual double getMoleculeSigma() const { return 1.0; }

  /** returns the epsilon parameter of the LJ potential. We always return 1.0,
   * see also doc/units.html for LJ style. */
  virtual double getMoleculeEpsilon() const { return 1.0; }

  /** sets a random velocity in the vector 'initialVelocity'. This velocity is
   * sampled from a Maxwellian assuming a mean flow velocity 'meanVelocity' and
   * a temperature 'temperature' of the fluid.
   */
  virtual void getInitialVelocity(const tarch::la::Vector<dim, double>& meanVelocity, const double& kB, const double& temperature,
                                  tarch::la::Vector<dim, double>& initialVelocity) const {
    const double stdDeviation = std::sqrt(dim * kB * temperature / getMoleculeMass()); // temperature-based std-deviation of
                                                                                       // gaussian distr.

    tarch::la::Vector<dim, double> ran(0.0); // buffer for one gaussian and dim-1 uniform random numbers
                                             // (latter ones are used for random rotation of velocity vector)
    ran[0] = _random->gaussian();
    for (unsigned int d = 1; d < dim; d++)
      ran[d] = _random->uniform();

    // init 2D/3D vector with avg meanVelocity and respective fluctuation
    if (dim == 2) {
      initialVelocity[0] = meanVelocity[0] + stdDeviation * (ran[0] * std::cos(ran[1]));
      initialVelocity[1] = meanVelocity[1] + stdDeviation * (ran[0] * std::sin(ran[1]));
    } else {
      initialVelocity[0] = meanVelocity[0] + stdDeviation * (ran[0] * std::sin(ran[1]) * std::cos(ran[2]));
      initialVelocity[1] = meanVelocity[1] + stdDeviation * (ran[0] * std::sin(ran[1]) * std::sin(ran[2]));
      initialVelocity[2] = meanVelocity[2] + stdDeviation * (ran[0] * std::cos(ran[1]));
    }
  }

  /** deletes the molecule from the MD simulation. This only handles the
   * process-local deletion. The synch. is carried out afterwards in
   *  synchronizeMoleculesAfterMassModification.
   */
  virtual void deleteMoleculeFromMDSimulation(const coupling::interface::Molecule<dim>& molecule, LAMMPS_NS::MamicoCell& cell) {
    const tarch::la::Vector<dim, double> pos = molecule.getPosition(); // position of molecule
    Atom* atom = _lmp->atom;

    // ------ the following code was copied and modified from delete_atoms.cpp
    // ------------------ delete atom from lammps store state before delete
    // bigint natoms_previous = atom->natoms;

    // delete local atom
    // reset nlocal
    AtomVec* avec = atom->avec;
    int nlocal = atom->nlocal;

    int i = 0;
    while (i < nlocal) {
      bool found = true;
      for (unsigned int d = 0; d < dim; d++) {
        found = found && tarch::la::equals(pos[d], atom->x[i][d], _tolerance);
      }

      if (found) {
        // do lammps stuff...
        avec->copy(nlocal - 1, i, 1);
        nlocal--;
      } else
        i++;
    }

    atom->nlocal = nlocal;
    // -------------------------------------------------------------------------------------------

    // delete molecule from list in mamico cells; as global enumeration may
    // change, we need to completely re-new the cell sorting in the NON-GHOST
    // cells
    // -> we cannot put this to synchronizeMoleculesAfterMassModification(),
    // since the local numbering of existing atoms might
    //    have already changed as well...
    const coupling::IndexConversion<dim>& indexConversion =
        coupling::interface::MamicoInterfaceProvider<LAMMPS_NS::MamicoCell, dim>::getInstance().getCouplingCellService()->getIndexConversion();
    _sorting.updateNonGhostCells(indexConversion);
  }

  /** adds the molecule to the MD simulation. */
  virtual void addMoleculeToMDSimulation(const coupling::interface::Molecule<dim>& molecule) {
    const coupling::IndexConversion<dim>& indexConversion =
        coupling::interface::MamicoInterfaceProvider<MamicoCell, dim>::getInstance().getCouplingCellService()->getIndexConversion();

    const int nlocal_previous = _lmp->atom->nlocal;              // no. of atoms before insertion
    tarch::la::Vector<dim, double> pos = molecule.getPosition(); // position of molecule
    double posVec[3] = {0.0, 0.0, 0.0};
    for (unsigned int d = 0; d < dim; d++)
      posVec[d] = pos[d];
    for (unsigned int d = dim; d < 3; d++)
      posVec[d] = 0.5 * (_lmp->domain->boxhi[d] + _lmp->domain->boxlo[d]);
    const tarch::la::Vector<dim, double> vel = molecule.getVelocity(); // velocity of molecule
    const tarch::la::Vector<dim, double> f = molecule.getForce();      // force of molecule

    // the following lines are extracted and adapted from create_atoms
    // ------------------------------------ this should yield a process-LOCAL
    // insertion of the atom
    _lmp->atom->avec->create_atom(_atomType, posVec);
    int nlocal = _lmp->atom->nlocal;
    for (int m = 0; m < _lmp->modify->nfix; m++) {
      Fix* fix = _lmp->modify->fix[m];
      if (fix->create_attribute)
        for (int i = nlocal_previous; i < nlocal; i++)
          fix->set_arrays(i);
    }
    // ----------------------------------------------------------------------------------------------------

    // find index of this molecule; since we expect that the molecule is
    // inserted at the last position, we start searching from the back
    int indexAtom = -1;
    for (int i = nlocal - 1; i > -1; i--) {
      bool found = true;
      for (unsigned int d = 0; d < dim; d++) {
        found = found && tarch::la::equals(_lmp->atom->x[i][d], pos[d], _tolerance);
      }
      if (found) {
        indexAtom = i;
        break;
      }
    }
    // if the molecules does not exist on this rank, we just jump out of this
    // method...
    if (indexAtom == -1) {
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
      std::cout << "Molecule not found on rank " << indexConversion.getThisRank() << std::endl;
#endif
      return;
    }

    // otherwise: write velocity and forces at respective positions
    for (unsigned int d = 0; d < dim; d++) {
      _lmp->atom->v[indexAtom][d] = vel[d];
      _lmp->atom->f[indexAtom][d] = f[d];
    }

    // add molecule to list in mamico cells; as global enumeration may change,
    // we need to completely re-new the cell sorting in all NON-GHOST cells
    // -> we cannot put this to synchronizeMoleculesAfterMassModification(),
    // since the local numbering of existing atoms might
    //    have already changed as well...
    _sorting.updateNonGhostCells(indexConversion);
  }

  /** sets up the potential energy landscape over the domain spanned by
   * indexOfFirstCouplingCell and rangeCoordinates. The first vector denotes
   * the position of the lower,left,front corner of the domain, rangeCoordinates
   * the number of coupling cells in each spatial direction in which the
   *  potential energy needs to be computed.
   */
  virtual void setupPotentialEnergyLandscape(const tarch::la::Vector<dim, unsigned int>& indexOfFirstCouplingCell,
                                             const tarch::la::Vector<dim, unsigned int>& rangeCouplingCells,
                                             const tarch::la::Vector<dim, unsigned int>& linkedCellsPerCouplingCell) {
    // nop, potential energy is evaluated for each atom separately when calling
    // getPotentialEnergy()
  }

  /** returns the local index vector (w.r.t. lexicographic ordering of the
   * linked cells in the MD simulation, see getLinkedCell()) for the linked cell
   * that the position 'position' belongs to. This method is crucial for the
   * USHER particle insertion scheme as we need to loop over all neighboured
   * linked cells and the cell itself to determine potential energy and forces
   * acting on the molecule at position 'position'.
   */
  virtual tarch::la::Vector<dim, unsigned int> getLinkedCellIndexForMoleculePosition(const tarch::la::Vector<dim, double>& position) {
    // currently, we have one mamico(linked) cell per coupling cell, so we
    // can just use the indexing of the coupling cells of MaMiCo note:
    // typically, one should not really access other coupling modules from an
    // interface; however, since this method is expected to be only called AFTER
    // all initialisation processes are finished, this is safe
    const coupling::IndexConversion<dim>& indexConversion =
        coupling::interface::MamicoInterfaceProvider<MamicoCell, dim>::getInstance().getCouplingCellService()->getIndexConversion();
    tarch::la::Vector<dim, unsigned int> vectorIndex = indexConversion.getGlobalVectorCellIndex(position);
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
    std::cout << "Linked cell index for position " << position << ": " << vectorIndex << "( global index), corresponds to local vector index of ";
    std::cout << indexConversion.convertGlobalToLocalVectorCellIndex(vectorIndex) << std::endl;
#endif
    return indexConversion.convertGlobalToLocalVectorCellIndex(vectorIndex);
  }

  /** assumes that a molecule is placed somewhere inside the linked cell at
   * index 'linkedCellIndex' and computes the force and potential energy
   * contributions from all molecules in the same linked cell and the
   * neighboured linked cells onto this molecule. The molecule "molecule" is
   * considered to NOT be part of the simulation domain and thus the linked
   * cells. Therefore, another molecule inside the linked cells may even
   * coincide with "molecule". The results are stored within the molecule.
   */
  virtual void calculateForceAndEnergy(coupling::interface::Molecule<dim>& molecule) {
    // helper variables and molecule position
    const coupling::IndexConversion<dim>& indexConversion =
        coupling::interface::MamicoInterfaceProvider<LAMMPS_NS::MamicoCell, dim>::getInstance().getCouplingCellService()->getIndexConversion();
    const tarch::la::Vector<dim, double> position1 = molecule.getPosition();
    const tarch::la::Vector<dim, unsigned int> cellIndex = getLinkedCellIndexForMoleculePosition(position1);

    // lennard jones parameters
    const double sigma6 = pow(getMoleculeSigma(), 6.0);
    const double epsilon = getMoleculeEpsilon();
    const double sigma6OverCutoff6 = sigma6 / (_cutOffRadiusSquared * _cutOffRadiusSquared * _cutOffRadiusSquared);
    const double cutOffEnergy = 4.0 * epsilon * sigma6OverCutoff6 * (sigma6OverCutoff6 - 1.0);

    // set force and energy to zero
    tarch::la::Vector<dim, double> force(0.0);
    double potEnergy = 0.0;

    // init start and end coordinates for neighbour search
    tarch::la::Vector<3, unsigned int> start(0);
    tarch::la::Vector<3, unsigned int> end(1);
    tarch::la::Vector<3, unsigned int> loop(0);
    for (unsigned int d = 0; d < dim; d++) {
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
      if (cellIndex[d] == 0) {
        std::cout << "ERROR "
                     "MamicoLammpsMDSolverInterface::calculateForceAndEnergy(): "
                     "Cannot compute force/energy for molecules in ghost layer!";
        std::cout << std::endl;
        exit(EXIT_FAILURE);
      }
#endif
      start[d] = cellIndex[d] - 1;
      end[d] = start[d] + 3;
    }

    for (loop[2] = start[2]; loop[2] < end[2]; loop[2]++) {
      for (loop[1] = start[1]; loop[1] < end[1]; loop[1]++) {
        for (loop[0] = start[0]; loop[0] < end[0]; loop[0]++) {
          LAMMPS_NS::MamicoCell& cell = _sorting.getMamicoCell(indexConversion.getLocalCellIndex(coupling::initDimVector<dim>(loop)));
          coupling::interface::MoleculeIterator<MamicoCell, dim>* it = getMoleculeIterator(cell);
          for (it->begin(); it->continueIteration(); it->next()) {
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
            std::cout << "MamicoLammpsMDSolverInterface::calculateForceAndEnergy(): "
                         "Compute force/energy with molecule at "
                      << it->getConst().getPosition() << std::endl;
#endif
            const tarch::la::Vector<dim, double> rij = it->getConst().getPosition() - position1;
            const double rij2 = tarch::la::dot(rij, rij);

            // compute LJ force/energy and add the contributions, only if center
            // is within cut-off radius AND if this particle is not contained in
            // this region already. We need the latter condition to re-use this
            // function for the evaluation of the potential energy in the
            // molecule-wrapper (mamico_lammps_molecule)
            if (rij2 > _cutOffRadiusSquared) {
            } else {
              const double rij6 = rij2 * rij2 * rij2;
              force += 24.0 * epsilon / rij2 * (sigma6 / rij6) * (1.0 - 2.0 * (sigma6 / rij6)) * rij;
              potEnergy += 0.5 * (4.0 * epsilon * (sigma6 / rij6) * ((sigma6 / rij6) - 1.0) - cutOffEnergy);
            }
          }
          delete it;
        }
      }
    }

    // set force and energy in molecule
    molecule.setForce(force);
    molecule.setPotentialEnergy(potEnergy);
  }

  /** is called each time when MaMiCo tried to insert/ delete molecules from the
   * MD simulation. As a consequence, a synchronization between processes or
   * with local boundary data might be necessary. Example: If in the lower left
   * cell of a 2D MD simulation a molecule was inserted and we use periodic
   *  boundary conditions, we need to provide the inserted molecule (amongst
   * others) to the upper right cell by adding this molecule to the ghost
   * boundary layer (when using the builtin MD simulation). For the builtin MD
   * simulation, the implementation of this method clears all molecules from the
   * ghost layers and re-fills the ghost layers again.
   */
  virtual void synchronizeMoleculesAfterMassModification() {
    // compute global number of atoms
    bigint nblocal = _lmp->atom->nlocal;
    MPI_Allreduce(&nblocal, &_lmp->atom->natoms, 1, MPI_LMP_BIGINT, MPI_SUM, _lmp->world);
    // if it's too many global atoms, throw error
    if (_lmp->atom->natoms < 0 || _lmp->atom->natoms > MAXBIGINT)
      _lmp->error->all(FLERR, "Too many total atoms");
    // if it's too many atoms to assign tags, switch this off...
    if (_lmp->atom->natoms > MAXSMALLINT) {
      if (_lmp->comm->me == 0)
        _lmp->error->warning(FLERR, "Total atom count exceeds ID limit, "
                                    "atoms will not have individual IDs");
      _lmp->atom->tag_enable = 0;
    } else {
      const bool compress_flag = true; // we always want to compress our data TODO CHECK
      if (_lmp->atom->molecular == 0 && compress_flag) {
        int* tag = _lmp->atom->tag;
        for (int i = 0; i < _lmp->atom->nlocal; i++)
          tag[i] = 0;
        _lmp->atom->tag_extend();
      }
    }

    // set map
    if (_lmp->atom->map_style) {
      _lmp->atom->nghost = 0;
      _lmp->atom->map_init();
      _lmp->atom->map_set();
    }

    // exchange information and re-build lists (extracted and modified from
    // verlet.cpp) ------------------------------ regular communication vs
    // neighbor list rebuild -> we definitely want to have a complete re-build
    // of lists, so we
    //  choose the "else"-branch from verlet.cpp
    {
      Verlet* verlet = (Verlet*)_lmp->update->integrate;
      if (verlet == NULL) {
        std::cout << "ERROR synchronizeMoleculesAfterMassModification(): Could not "
                     "cast Integrate* to Verlet*! Only verlet-integration supported!"
                  << std::endl;
        exit(EXIT_FAILURE);
      }
      if (_lmp->modify->n_pre_exchange)
        _lmp->modify->pre_exchange();
      if (_lmp->domain->triclinic)
        _lmp->domain->x2lamda(_lmp->atom->nlocal);
      _lmp->domain->pbc();
      if (_lmp->domain->box_change) {
        _lmp->domain->reset_box();
        _lmp->comm->setup();
        if (_lmp->neighbor->style)
          _lmp->neighbor->setup_bins();
      }
      _lmp->timer->stamp();
      _lmp->comm->exchange();
      // if (sortflag && ntimestep >= atom->nextsort) atom->sort();
      _lmp->comm->borders();
      if (_lmp->domain->triclinic)
        _lmp->domain->lamda2x(_lmp->atom->nlocal + _lmp->atom->nghost);
      _lmp->timer->stamp(Timer::COMM);
      if (_lmp->modify->n_pre_neighbor)
        _lmp->modify->pre_neighbor();
      _lmp->neighbor->build(1);
      _lmp->timer->stamp(Timer::NEB);
    }

    // one more sorting for all cells, since the molecules might have been
    // reordered in parallel exchange
    const coupling::IndexConversion<dim>& indexConversion =
        coupling::interface::MamicoInterfaceProvider<LAMMPS_NS::MamicoCell, dim>::getInstance().getCouplingCellService()->getIndexConversion();
    _sorting.updateAllCells(indexConversion);
  }

  /** is called each time when MaMiCo tried to insert momentum in the MD
   * simulation. For the builtin MD simulation, this method is empty as the
   * simulation does not need to synchronize changing momentum over the
   * processes within a timestep (only the positions and numbers of molecules in
   * possible ghost layers matter!). However, other solvers might need to
   * implement something in here.
   */
  virtual void synchronizeMoleculesAfterMomentumModification() {
    // nop here, this worked out for all of my settings so far...
  }

  /** returns the timestep of the MD simulation */
  virtual double getDt() {
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
    std::cout << "Timestep in MD: " << _lmp->update->dt << std::endl;
#endif
    return _lmp->update->dt;
  }

  /** returns a new molecule iterator for a certain linked cell */
  virtual coupling::interface::MoleculeIterator<LAMMPS_NS::MamicoCell, dim>* getMoleculeIterator(LAMMPS_NS::MamicoCell& cell) {
    LAMMPS_NS::MoleculeInformation info;
    if (cell.isGhostCell()) {
      info._x = _sorting.getGhostAtomPositions();
      info._v = NULL;
      info._f = NULL;
    } else {
      info._x = _lmp->atom->x;
      info._v = _lmp->atom->v;
      info._f = _lmp->atom->f;
    }
    return new MamicoLammpsMoleculeIterator<dim>(info, _cutOffRadius, cell);
  }

  /** update all cell information -> forward call to _sorting. This is required
     once per time step before the actual coupling */
  void updateAllCells(const coupling::IndexConversion<dim>& indexConversion) { _sorting.updateAllCells(indexConversion); }

  /** prints molecules in all cells/inner cells/only ghost cells. For debugging
   * purposes only. */
  void printMolecules(const coupling::IndexConversion<dim>& indexConversion, typename LAMMPS_NS::Sorting<dim>::PrintType printType) {
    _sorting.printMolecules(indexConversion, printType);
  }

private:
  LAMMPS* _lmp;
  Sorting<dim> _sorting;
  const double _cutOffRadiusSquared;
  const double _cutOffRadius;
  const unsigned int _atomType;  // this is the atom type used in the coupled simulation.
                                 // Currently, it is hard-coded and set to 1.
  const unsigned int _precision; // number of digits for number-to-string conversion
  const double _tolerance;       // tolerance in number-to-string conversion
  RanPark* _random;              // random number generator
};

} // namespace LAMMPS_NS
#endif // MAMICO_LMP_LAMMPS_MD_SOLVER_INTERFACE_H
