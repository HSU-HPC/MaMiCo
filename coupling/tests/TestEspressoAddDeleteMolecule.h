// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_TESTS_TESTESPRESSOADDDELETEMOLECULE_H
#define _MOLECULARDYNAMICS_COUPLING_TESTS_TESTESPRESSOADDDELETEMOLECULE_H

#include "coupling/interface/impl/Espresso/EspressoMDMolecule.h"
#include "coupling/interface/impl/Espresso/EspressoMDMoleculeIterator.h"
#include "coupling/interface/impl/Espresso/EspressoMDSolverInterface.h"
#include "coupling/tests/TestEspresso.h"
#include "particle_data.hpp"
#include "tarch/la/Vector.h"

/** test class to check addition and deletion of molecule in espresso. It
 * derives from the scenario provided in TestEspresso.h, containing 8 particles
 * in 2*2*2 domain. We call the add particle and delete particle methods
 * provided in the interface classes and check the number of particles remaining
 * in our domain. Also this test class check the synchronization of ghosts cell
 *  due to addition and deletion of particles
 *  @author Rahul Arora
 */

class TestEspressoAddDeleteMolecule : public TestEspresso {
public:
  TestEspressoAddDeleteMolecule(std::string name, int argc, char **argv) : TestEspresso(name, argc, argv) {}
  ~TestEspressoAddDeleteMolecule() {}
  virtual void run() {

    int type = 1;

    // If type is one, we run a delete particle test, otherwise we run a add
    // particle test and check the ghost cells for synchronization

    if (type == 1) {
      loadEspressoTestConfiguration();
      Cell *cell;
      cell = local_cells.cell[0];
      deleteParticleTest(cell);
    }

    else {
      // Test particle addition
      loadEspressoTestConfiguration();
      addParticleTest();

      // Test Synchronize molecules after mass modification
      synchronizeMoleculesAfterMassModificationTest();
    }
  }

private:
  /* This test function delete the first molecule in cell 0 at time t = 0 and
   * then iterates over all the molecules, and counts the number of remaining
   * molecules */

  void deleteParticleTest(ParticleList *cell) {
    ParticleList &temp_cell = *cell;
    Particle *part;
    part = cell->part;

    coupling::interface::EspressoMDSolverInterface delete_molecule_test;
    coupling::interface::EspressoMDMolecule coupling_part(part);

    // Call the method deleteMoleculeFromMDSimulation provided in
    // EspressoMDSolverInterface
    delete_molecule_test.deleteMoleculeFromMDSimulation(coupling_part, temp_cell);

    // Iterate over all the particles in the local domain (excluding ghost
    // cells) in Espresso
    Cell *celltemp;
    int c, i, np, cnt = 0;
    for (c = 0; c < local_cells.n; c++) {
      celltemp = local_cells.cell[c];
      np = celltemp->n;
      for (i = 0; i < np; i++) {
        cnt++;
      }
    }
    if (cnt == 7) {
      std::cout << "The test was successful" << std::endl;
    } else {
      std::cout << "The test was not successful" << std::endl;
    }
  }

  /* This test function adds a molecule to the simulation domain at time t=0 and
   * at the location (0.05,0.05,0.05). Then it iterates over all the molecules
   * in the domain and calculates the number of molecules */

  void addParticleTest() {
    coupling::interface::EspressoMDSolverInterface _test;

    Particle newpart;
    newpart.p.identity = max_seen_particle + 1;
    tarch::la::Vector<3, double> pos(0.05);
    tarch::la::Vector<3, double> vel(0.0);
    for (unsigned int i = 0; i < 3; i++) {
      newpart.r.p[i] = pos[i];
      newpart.m.v[i] = vel[i];
      newpart.f.f[i] = vel[i];
    }
    newpart.p.type = 0;
    coupling::interface::EspressoMDMolecule new_part(&newpart);
    new_part.setPosition(pos);
    new_part.setVelocity(vel);
    new_part.setForce(vel);

    // Call the method addMoleculeToMDSimulation provided in
    // EspressoMDSolverInterface
    _test.addMoleculeToMDSimulation(new_part);

    // Iterate over all the particles in the local domain (excluding ghost
    // cells) in Espresso
    Cell *cell;
    int c, i, np, cnt = 0;
    for (c = 0; c < local_cells.n; c++) {
      cell = local_cells.cell[c];
      np = cell->n;
      for (i = 0; i < np; i++) {
        cnt++;
      }
    }
    if (cnt == 9) {
      std::cout << "The test was successful" << std::endl;
    } else {
      std::cout << "The test was not successful" << std::endl;
    }
  }

  /* This function, first call the addParticleTest() function and then updates
   * the ghost layers and prints them to see if the added particle has been
   * updated in the ghost layers */

  void synchronizeMoleculesAfterMassModificationTest() {
    coupling::interface::EspressoMDSolverInterface _test;
    // Call the method synchronizeMoleculesAfterMassModification provided in
    // EspressoMDSolverInterface, which call the Espresso synchronization
    // function
    _test.synchronizeMoleculesAfterMassModification();

    Cell *cell;
    int c, i, np, cnt = 0;
    Particle *part;

    std::cout << ghost_cells.n << " " << local_cells.n << std::endl;

    // Iterate over all the ghost cells in Espresso, and calculate the total no.
    // of particles in them
    for (c = 0; c < ghost_cells.n; c++) {
      cell = ghost_cells.cell[c];
      part = cell->part;
      np = cell->n;
      for (i = 0; i < np; i++) {
        cnt++;
      }
    }
    if (cnt > 56) {
      std::cout << "The test was successful with count " << cnt << std::endl;
    } else {
      std::cout << "The test was not successful with count " << cnt << std::endl;
    }
  }
};

#endif // _MOLECULARDYNAMICS_COUPLING_TESTS_TESTESPRESSOADDDELETEMOLECULE_H
