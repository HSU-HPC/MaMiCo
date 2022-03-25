// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_TESTS_TESTESPRESSOFORCEENERGYCALCULATION_H_
#define _MOLECULARDYNAMICS_COUPLING_TESTS_TESTESPRESSOFORCEENERGYCALCULATION_H_

#include "coupling/interface/impl/Espresso/EspressoMDMolecule.h"
#include "coupling/interface/impl/Espresso/EspressoMDMoleculeIterator.h"
#include "coupling/interface/impl/Espresso/EspressoMDSolverInterface.h"
#include "coupling/tests/TestEspresso.h"
#include "integrate_tcl.hpp"
#include "particle_data.hpp"
#include "tarch/la/Vector.h"

/** test class for force and potential energy calculation methods provided
 * inside the EspressoMDSolverInterface. In this class, I derive from
 * TestEspresso and have a very simple 2*2*2 grid, with 8 particles, all located
 * at the center of each cell I add a single particle at the center of the
 * domain. In this simple scenario, the net force acting on the new particle
 * should be zero, while its potential energy should be greater than zero. Also
 * I provide a very simple test to the method
 * getLinkedCellIndexForMoleculePositionTest
 *  @author Rahul Arora
 */

class TestEspressoForceEnergyCalculation : public TestEspresso {
public:
  TestEspressoForceEnergyCalculation(std::string name, int argc, char **argv) : TestEspresso(name, argc, argv) {}
  ~TestEspressoForceEnergyCalculation() {}
  virtual void run() {
    loadEspressoTestConfiguration();

    // Test getLinkedCellIndex from molecule position
    getLinkedCellIndexForMoleculePositionTest();

    // Test force and potential energy calculation
    calculateForceandEnergyTest();
  }

private:
  /* This test function is used on a test domain ,which consists of 8 cells and
   * each cell conatin a molecule in the center of the cell, in this function
   * first we add an molecule at the center of the domain , equidistant from all
   * the molecules at time t=0. Then we compute the force acting on that
   * molecule as well as its potential energy. The force should be zero ( the
   * forces are equal in magnitude and opposite in direction), while the
   * potential energy has to be non zero */

  void calculateForceandEnergyTest() {
    coupling::interface::EspressoMDSolverInterface test;
    Particle newpart;
    newpart.p.identity = max_seen_particle + 1;

    tarch::la::Vector<3, double> pos(1.0);
    tarch::la::Vector<3, double> vel(0.0);

    for (unsigned int i = 0; i < 3; i++) {
      newpart.r.p[i] = pos[i];
      newpart.m.v[i] = vel[i];
    }
    newpart.p.type = 0;

    coupling::interface::EspressoMDMolecule new_part(&newpart);
    new_part.setPosition(pos);

    for (int c = 0; c < local_cells.n; c++) {
      ParticleList *list = local_cells.cell[c];
      Particle *part;
      part = list->part;
      int number = list->n;
      for (int j = 0; j < number; j++) {
        std::cout << "Positon of particle " << j << " in cell " << c << " is " << part[j].r.p[0] << " " << part[j].r.p[1] << " " << part[j].r.p[2] << std::endl;
      }
    }

    // Call to the method calculateForceAndEnergy, provided in
    // EspressoMDSolverInterface
    test.calculateForceAndEnergy(new_part);

    tarch::la::Vector<3, double> force(0.0);
    force = new_part.getForce();
    std::cout << "the force acting on the particle " << force[0] << " " << force[1] << " " << force[2] << std::endl;
  }

  /* This function calcultes the linked cell index for each molecule in the
   * domain (for local cells only, no ghosts) and prints it out */
  void getLinkedCellIndexForMoleculePositionTest() {
    Cell *cell;
    Particle *part;
    cell = local_cells.cell[0];
    part = cell->part;
    coupling::interface::EspressoMDMolecule coupling_part(part);

    coupling::interface::EspressoMDSolverInterface test;
    tarch::la::Vector<3, unsigned int> index(1);
    tarch::la::Vector<3, double> position(0.0);

    position = coupling_part.getPosition();
    // Call to the method getLinkedCellIndexForMoleculePosition, provided in
    // EspressoMDSolverInterface
    index = test.getLinkedCellIndexForMoleculePosition(position);
    if ((index[0] == 1) && (index[1] == 1) && (index[2] == 1)) {
      std::cout << "getLinkedCellIndexForMoleculePositionTest was successful " << std::endl;
    } else {
      std::cout << "getLinkedCellIndexForMoleculePositionTest was not successful" << std::endl;
    }
  }
};

#endif //_MOLECULARDYNAMICS_COUPLING_TESTS_TESTESPRESSOFORCEENERGYCALCULATION_H_
