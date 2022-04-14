// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#include "coupling/tests/TestEspresso.h"
#include "coupling/tests/TestEspressoAddDeleteMolecule.h"
#include "coupling/tests/TestEspressoForceEnergyCalculation.h"
#include "coupling/tests/TestEspressoMDMolecule.h"
#include "coupling/tests/TestEspressoMoleculeIterator.h"

#include <mpi.h>

void runTest(Test* test) {
  if (test == NULL) {
    std::cout << "ERROR executeTest: test==NULL!" << std::endl;
    exit(EXIT_FAILURE);
  }
  test->run();
  delete test;
}

int main(int argc, char* argv[]) {
  int thisTest = atoi(argv[1]);
  switch (thisTest) {
  case 0:
    runTest(new TestEspresso("TestEspresso", argc, argv));
    break;
  case 1:
    runTest(new TestEspressoMDMolecule("TestEspressoMDMolecule", argc, argv));
    break;
  case 3:
    runTest(new TestEspressoAddDeleteMolecule("TestEspressoAddDeleteMolecule", argc, argv));
    break;
  case 4:
    runTest(new TestEspressoForceEnergyCalculation("TestEspressoForceEnergyCalculation", argc, argv));
    break;
  case 5:
    runTest(new TestEspressoMoleculeIterator("TestEspressoMoleculeIterator", argc, argv));
    break;
  default:
    std::cout << "Test number " << thisTest << " out of range..." << std::endl;
  }

  return 0;
}
