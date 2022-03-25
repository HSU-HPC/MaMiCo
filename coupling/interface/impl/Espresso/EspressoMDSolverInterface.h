// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_ESPRESSOMDSOLVERINTERFACE_H_
#define _MOLECULARDYNAMICS_COUPLING_ESPRESSOMDSOLVERINTERFACE_H_

#include "cells.hpp"
#include "coupling/interface/MDSolverInterface.h"
#include "coupling/interface/impl/Espresso/EspressoMDMoleculeIterator.h"
#include "domain_decomposition.hpp"
#include "forces.hpp"
#include "particle_data.hpp"
#include "tarch/la/Vector.h"
#include <cstdlib>

namespace coupling {
namespace interface {

/** Interface for Espresso Simulation
 *  @author Rahul Arora
 */
class EspressoMDSolverInterface : public MDSolverInterface<ParticleList, 3> {
public:
  ~EspressoMDSolverInterface() {}

  ParticleList &getLinkedCell(const tarch::la::Vector<3, unsigned int> &macroscopicCellIndex,
                              const tarch::la::Vector<3, unsigned int> &linkedCellInMacroscopicCell,
                              const tarch::la::Vector<3, unsigned int> &linkedCellsPerMacroscopicCell, const coupling::IndexConversion<3> &indexConversion) {
    for (unsigned int d = 0; d < 3; d++) {
      if (macroscopicCellIndex[d] == 0) {
        std::cout << "ERROR EspressoMDSolverInterface::getLinkedCell(): "
                     "macroscopic cell outside range for linked cells!"
                  << std::endl;
        exit(EXIT_FAILURE);
      }
    }
    tarch::la::Vector<3, unsigned int> index(0);

    for (unsigned int d = 0; d < 3; d++) {
      index[d] = index[d] + (macroscopicCellIndex[d] - 1) * linkedCellsPerMacroscopicCell[d] + linkedCellInMacroscopicCell[d];
    }
    int id = (index[0] + dd.cell_grid[0] * (index[1] + dd.cell_grid[1] * index[2]));

    ParticleList *cell;
    cell = local_cells.cell[id];
    return (*cell);
  }

  tarch::la::Vector<3, double> getGlobalMDDomainSize() const {
    tarch::la::Vector<3, double> domainsize(0.0);
    for (int i = 0; i < 3; i++) {
      domainsize[i] = box_l[i];
    }
    return domainsize;
  }

  tarch::la::Vector<3, double> getGlobalMDDomainOffset() const {
    tarch::la::Vector<3, double> domainoffset(0.0);
    /*        	for(int i=0; i<3; i++) {
                            domainoffset[i]=my_left[i];
                    }*/
    return domainoffset;
  }
  double getMoleculeMass() const { return 1.0; }

  double getKB() const { return 1.0; }

  double getMoleculeSigma() const { return 1.0; }

  double getMoleculeEpsilon() const { return 1.0; }

  void getInitialVelocity(const tarch::la::Vector<3, double> &meanVelocity, const double &kB, const double &temperature,
                          tarch::la::Vector<3, double> &initialVelocity) const {
    tarch::la::Vector<3, double> randomNumbers(0.0);
    double mass = 1.0; // use getMoleculeMass() when it will be defined

    double u, v, r, c, num;
    num = 0.0;
    while (r > 0 && r < 1) {
      u = ((double)rand() / (RAND_MAX)) * 2 - 1;
      v = ((double)rand() / (RAND_MAX)) * 2 - 1;
      r = u * u + v * v;
      c = sqrt(-2 * log(r) / r);
      num = u * c;
    }

    randomNumbers[0] = num;
    double stdDeviation = std::sqrt(3 * kB * temperature / mass);
    for (unsigned int d = 1; d < 3; d++) {
      randomNumbers[d] = 2.0 * PI * ((double)rand()) / ((double)RAND_MAX);
    }

    initialVelocity[0] = meanVelocity[0] + stdDeviation * (randomNumbers[0] * std::sin(randomNumbers[1]) * std::cos(randomNumbers[2]));
    initialVelocity[1] = meanVelocity[1] + stdDeviation * (randomNumbers[0] * std::sin(randomNumbers[1]) * std::sin(randomNumbers[2]));
    initialVelocity[2] = meanVelocity[2] + stdDeviation * (randomNumbers[0] * std::cos(randomNumbers[1]));
  }

  void deleteMoleculeFromMDSimulation(const coupling::interface::Molecule<3> &molecule, ParticleList &cell) {
    int cnt = 0;
    int number = cell.n;
    if (number == 0) {
      std::cout << "ERROR: Cell contains no particle, cannot delete the particle " << std::endl;
      exit(EXIT_FAILURE);
      return;
    }
    Particle part;

    tarch::la::Vector<3, double> moleculePosition = molecule.getPosition();
    while (cnt < number) {
      part = cell.part[cnt];
      tarch::la::Vector<3, double> itPosition(0.0);
      for (int i = 0; i < 3; i++) {
        itPosition[i] = part.r.p[i];
      }

      // Delete a particle at a particular location using remove_particle()
      if (moleculePosition == itPosition) { // tolerance level
        int id = part.p.identity;
        int flag = remove_particle(id);
        if (flag == 1) {
          std::cout << "ERROR: remove_particle could not delete molecule at position " << moleculePosition << "!" << std::endl;
          exit(EXIT_FAILURE);
        }
        return;
      }
      cnt++;
    }

    std::cout << "ERROR: Could not delete molecule at position " << moleculePosition << "!" << std::endl;
    exit(EXIT_FAILURE);
  }

  void addMoleculeToMDSimulation(const coupling::interface::Molecule<3> &molecule) {
    tarch::la::Vector<3, double> position = molecule.getPosition();
    tarch::la::Vector<3, double> velocity = molecule.getVelocity();
    tarch::la::Vector<3, double> force = molecule.getForce();

    int id = max_seen_particle + 1;
    double pos[3], vel[3], f[3];
    for (int i = 0; i < 3; i++) {
      pos[i] = position[i];
    }

    // Places a particle with id onto a given position
    int flag = place_particle(id, pos);

    if (flag == 1) {
      for (int i = 0; i < 3; i++) {
        vel[i] = velocity[i];
      }

      for (int i = 0; i < 3; i++) {
        f[i] = force[i];
      }

      int vflag = set_particle_v(id, vel);
      int fflag = set_particle_f(id, f);
      int tflag = set_particle_type(id, 0);
    } else {
      std::cout << "ERROR: particle addition was not successful " << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  void setupPotentialEnergyLandscape(const tarch::la::Vector<3, unsigned int> &indexOfFirstMacroscopicCell,
                                     const tarch::la::Vector<3, unsigned int> &rangeMacroscopicCells,
                                     const tarch::la::Vector<3, unsigned int> &linkedCellsPerMacroscopicCell) {
    // Empty
  }

  // Returns the index of the linked cell the particle is located in
  tarch::la::Vector<3, unsigned int> getLinkedCellIndexForMoleculePosition(const tarch::la::Vector<3, double> &position) {
    double lpos;
    tarch::la::Vector<3, unsigned int> cpos(0);
    for (unsigned int i = 0; i < 3; i++) {
      lpos = position[i] - my_left[i];
      // cpos[i] = (int)(lpos*dd.inv_cell_size[i])+1;
      cpos[i] = ceil(lpos * dd.inv_cell_size[i]);
    }
    std::cout << "Position " << position << " Index " << cpos << std::endl;
    return cpos;
  }

  // Compute the force and potential energy of the molecule using Espresso
  // methods calc_non_bonded_pair_energy and calc_non_bonded_pair_force
  void calculateForceAndEnergy(coupling::interface::Molecule<3> &molecule) {
    // Create a new temporary particle in Espresso and set its velocity and
    // position
    Particle temp_part;
    temp_part.p.identity = max_seen_particle + 1;

    tarch::la::Vector<3, double> pos(0.0);
    tarch::la::Vector<3, double> vel(0.0);

    pos = molecule.getPosition();
    vel = molecule.getVelocity();

    for (unsigned int i = 0; i < 3; i++) {
      temp_part.r.p[i] = pos[i];
      temp_part.m.v[i] = vel[i];
    }
    temp_part.p.type = 0; // Specific for LJ case, type 0

    double potentialEnergy = 0.0;
    double force[3];
    for (unsigned int i = 0; i < 3; i++) {
      force[i] = 0.0;
    }

    // Get the linked cell the particle is located in
    tarch::la::Vector<3, unsigned int> linkedCellIndex = getLinkedCellIndexForMoleculePosition(pos);
    unsigned int c = ((linkedCellIndex[0] - 1) + dd.cell_grid[0] * ((linkedCellIndex[1] - 1) + dd.cell_grid[1] * (linkedCellIndex[2] - 1)));

    if ((c < 0) || (c > local_cells.n)) {
      std::cout << "ERROR: could not find the cell the particle is located in" << std::endl;
      exit(EXIT_FAILURE);
    }

    // pl is the cell which contains the particle Molecule
    ParticleList *pl;
    pl = local_cells.cell[c];

    IA_Neighbor *neighbor;
    Particle *p2;
    int np2;

    double dist2, vec21[3], t1[3], t2[3];
    for (unsigned int i = 0; i < 3; i++) {
      t1[i] = 0.0;
      t2[i] = 0.0;
    }

    unsigned int ctr = 0;
    // loop over all the neighbor cells and calculate potential energy and force
    // contribution from each molecule in the cell
    for (unsigned int n = 0; n < dd.cell_inter[c].n_neighbors; n++) {
      neighbor = &dd.cell_inter[c].nList[n];
      p2 = neighbor->pList->part;
      np2 = neighbor->pList->n;

      for (unsigned int j = 0; j < np2; j++) {
        // The functions calc_non_bonded_pair_energy and
        // calc_non_bonded_pair_force have a check, if the two particles sent to
        // them have the same location or are the same, force or potential
        // energy is not calculated then.

        dist2 = distance2vec(temp_part.r.p, p2[j].r.p, vec21);
        IA_parameters *ia_params = get_ia_param(temp_part.p.type, p2->p.type);
        Particle *temp = &(p2[j]);
        Particle *p = &temp_part;
        potentialEnergy += calc_non_bonded_pair_energy(p, temp, ia_params, vec21, sqrt(dist2), dist2);
        calc_non_bonded_pair_force(p, temp, ia_params, vec21, sqrt(dist2), dist2, force, t1, t2);
      }
    }

    tarch::la::Vector<3, double> forcevec(0.0);
    for (unsigned int i = 0; i < 3; i++) {
      forcevec[i] = force[i];
    }

    molecule.setForce(forcevec);
    molecule.setPotentialEnergy(potentialEnergy);
  }

  // Synchronize Molecules after mass modification using cells_update_ghosts()
  // function in Espresso
  void synchronizeMoleculesAfterMassModification() { cells_update_ghosts(); }

  // Synchronize Molecules after Momentum modification is empty
  void synchronizeMoleculesAfterMomentumModification() {}

  double getDt() { return time_step; }

  // Returns a molecule iterator to the given cell
  coupling::interface::MoleculeIterator<ParticleList, 3> *getMoleculeIterator(ParticleList &cell) {
    return new coupling::interface::EspressoMDMoleculeIterator(cell);
  }
};
} // namespace interface
} // namespace coupling
#endif // _MOLECULARDYNAMICS_COUPLING_ESPRESSOMDSOLVERINTERFACE_H_
