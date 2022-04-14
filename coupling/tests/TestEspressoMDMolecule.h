// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_TESTS_TESTESPRESSOMDMOLECULE_H_
#define _MOLECULARDYNAMICS_COUPLING_TESTS_TESTESPRESSOMDMOLECULE_H_

#include "coupling/interface/impl/Espresso/EspressoMDMolecule.h"
#include "coupling/interface/impl/Espresso/EspressoMDMoleculeIterator.h"
#include "coupling/interface/impl/Espresso/EspressoMDSolverInterface.h"
#include "coupling/tests/TestEspresso.h"
#include "particle_data.hpp"
#include "tarch/la/Vector.h"
#include <iostream>
#include <tcl.h>

/** Test class for EspressoMDMolecule class and functions, tries to check if the datastaructures defined
 *  for Espresso in the coupling tool are correct or not. To do that, I create two lists, one which gets data
 *  from Espresso data structures and other from coupling tool data structures and compare the two lists
 *  @author Rahul Arora
 */

class TestEspressoMDMolecule : public TestEspresso {
private:
public:
  TestEspressoMDMolecule(std::string name, int argc, char** argv) : TestEspresso(name, argc, argv) {}
  virtual ~TestEspressoMDMolecule() {}

  virtual void run() {
    std::cout << "Test EspressoMDMolecule" << std::endl;

    loadEspressoTestConfiguration();

    Cell* cell;
    int j, c, np, cnt = 0;
    Particle* part;
    int flag = 1;

    for (c = 0; c < local_cells.n; c++) {
      compareParticleLists(c, flag);
    }
    if (flag == 1) {
      std::cout << "EspressoMDMolecule Test was successful " << std::endl;
    } else {
      std::cout << "EspressoMDMolecule Test was not successful " << std::endl;
    }
  }

  /* compareParticleList() compares the list of particles in each cell, one directly read from Espresso Simulation and other from the list 	sent to the
   * coupling tool (EspressoMDMolecule) */

  void compareParticleLists(int c, int& flag) {
    int j = 0;
    int counter = 0;

    // Define a molecule iterator and iterate over all the molecules in that particular cell
    ParticleList* list = local_cells.cell[c];
    coupling::interface::EspressoMDMoleculeIterator it(*list);
    Particle* part;
    part = list->part;
    int number = list->n;

    // I define two lists, listEspresso and listEspressoCouplingInterface and add the respective data to them
    std::vector<tarch::la::Vector<3, double>> listEspresso;
    std::vector<tarch::la::Vector<3, double>> listEspressoCouplingInterface;
    tarch::la::Vector<3, double> velocity(0.0);
    tarch::la::Vector<3, double> position(0.0);

    for (it.begin(); it.continueIteration(); it.next()) {
      velocity = it.get().getVelocity();
      position = it.get().getPosition();

      listEspressoCouplingInterface.push_back(velocity);
      listEspressoCouplingInterface.push_back(position);

      for (int i = 0; i < 3; i++) {
        velocity[i] = part[j].m.v[i];
        position[i] = part[j].r.p[i];
      }
      listEspresso.push_back(velocity);
      listEspresso.push_back(position);

      j++;
      counter++;
    }

    if (listEspresso.size() != (unsigned int)2 * number) {
      std::cout << "ERROR: particle list has wrong size!" << std::endl;
      flag = 0;
    }
    if (listEspresso != listEspressoCouplingInterface) {
      std::cout << "ERROR: the particle positions and velocites obtianed by the Coupling interfaces and Espresso simulation do not match " << std::endl;
      flag = 0;
    }
  }
};

#endif // _MOLECULARDYNAMICS_COUPLING_TESTS_TESTESPRESSOMDMOLECULE_H_
