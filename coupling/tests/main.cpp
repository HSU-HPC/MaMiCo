// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MAIN_CPP_
#define _MAIN_CPP_

#include "coupling/CouplingMDDefinitions.h"
#include "coupling/tests/CompileTest.h"
#include "coupling/tests/FromMD2MacroTest.h"
#include "coupling/tests/FromMacro2MDTest.h"
#include "coupling/tests/IndexConversionTest.h"
#include "coupling/tests/MacroscopicCellsTest.h"
#include <cstdlib>
#include <iostream>
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
#include <mpi.h>
#endif

/** executes a newly created test and deletes it immediately again. */
void runTest(Test* test) {
  if (test == NULL) {
    std::cout << "ERROR executeTest: test==NULL!" << std::endl;
    exit(EXIT_FAILURE);
  }
  test->run();
  delete test;
}

int main(int argc, char* argv[]) {

#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
  MPI_Init(&argc, &argv);
#endif

  // run tests
  runTest(new IndexConversionTest());
  runTest(new FromMD2MacroTest());
  runTest(new FromMacro2MDTest());
  runTest(new CompileTest());
  runTest(new MacroscopicCellsTest());

#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
  MPI_Finalize();
#endif

  return 0;
}

#endif // _MAIN_CPP_