// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#include <iostream>
#include <stdlib.h>

#include "mpi.h"

#include "coupling/tests/Test.h"

#include "TestMarDyn.h"
#include "TestMarDynAddDeleteMolecule.h"
#include "TestMarDynMolecule.h"
#include "TestMarDynMoleculeIterator.h"

/* executes a newly created test and deletes it immediately again. */
void runTest(Test *test) {
  if (test == NULL) {
    std::cout << "ERROR: test==NULL !" << std::endl;
    exit(EXIT_FAILURE);
  }
  test->run();
  delete test;
}

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);

  std::cout << "Testing MarDyn coupling interfaces: " << std::endl;

  // Test for simple get methods of the MD solver interface
  runTest(new TestMarDyn(argc, argv, "MarDyn Test"));

  // Test for the Molecule iterator interface
  runTest(new TestMarDynMoleculeIterator(argc, argv, "MarDyn Iterator Test"));

  // Test for the Molecule interface
  runTest(new TestMarDynMolecule(argc, argv, "MarDyn Molecule Test"));

  // Test for molecule insertion/deletion
  runTest(new TestMarDynAddDeleteMolecule(argc, argv,
                                          "MarDyn Deletion/Insertion Test"));

  std::cout << "Finished!" << std::endl;

  MPI_Finalize();

  return 0;
}
