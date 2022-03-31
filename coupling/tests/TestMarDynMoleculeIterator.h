// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef TESTMARDYNMOLECULEITERATOR_H_
#define TESTMARDYNMOLECULEITERATOR_H_

#include "TestMarDyn.h"
#include "coupling/interface/impl/MarDyn/MarDynMoleculeIterator.h"

#include "particleContainer/LinkedCells.h"

/**
 * 	Test for the MarDynMoleculeIterator interface:
 * 	Loops through all ParticleCells of the simulation, sets an
 * MarDynMoleculeIterator interface for each one,
 * 	then compares the number of molecules of each ParticleCell with the number
 * of molecules obtained through the interface.
 * 	@author Hanno Flohr
 */
class TestMarDynMoleculeIterator : public TestMarDyn {
public:
  TestMarDynMoleculeIterator(int argc, char **argv, std::string name) : TestMarDyn(argc, argv, name) {}
  virtual ~TestMarDynMoleculeIterator() {}

  virtual void run() {
    std::cout << "Initializing MarDyn MD solver interface.." << std::endl;
    // init simulation and set md solver interface
    this->loadMarDynTestConfiguration("mardyn_dummy_coupling.cfg", 1);

    std::cout << "Testing MarDynMoleculeIterator interface:" << std::endl;
    testIterator();
  }

protected:
  void testIterator() {
    unsigned int errors = 0;
    // get particle container from simulation
    LinkedCells *cells = (LinkedCells *)this->_marDyn->getMolecules();

    ParticleCell pc;
    for (int i = 0; i < 2744; i++) {
      pc = cells->getCell(i);
      MarDynCell cell = MarDynCell(&pc, cells->getCutoff());

      // set up molecule iterator for testing
      MarDynMoleculeIterator mdmi(cell);

      // number of molecules in current particle cell
      unsigned int moleculeCount = pc.getMoleculeCount();

      // determine number of molecules through interface
      unsigned int moleculesInInterface = 0;
      for (mdmi.begin(); mdmi.continueIteration(); mdmi.next()) {
        moleculesInInterface++;

        // std::cout << "pos: " << mdmi.getConst().getPosition() << std::endl;
        // std::cout << "vel: " << mdmi.getConst().getVelocity() << std::endl;
        // std::cout << "------" << std::endl;
      }

      if (moleculeCount != 0)
        std::cout << "Test: " << moleculeCount << " =?= " << moleculesInInterface << std::endl;

      if (moleculeCount != moleculesInInterface) {
        errors++;
        std::cout << "ERROR: number of molecules in current cell " << i << "should be " << moleculeCount
                  << ". The interface instead has: " << moleculesInInterface << std::endl;
      }
    }

    if (errors == 0)
      std::cout << "MarDynMoleculeIterator test was successful!" << std::endl;
    else
      std::cout << "MarDynMoleculeIterator test was not successful! (Errors: " << errors << ")" << std::endl;
  }
};

#endif /* TESTMARDYNMOLECULEITERATOR_H_ */
