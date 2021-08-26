// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MAIN_CPP_
#define _MAIN_CPP_

#include <iostream>
#include <cstdlib>
#include "coupling/tests/TestLammps.h"
#include "coupling/tests/TestLammpsMoleculeIterator.h"
#include "coupling/tests/TestLammpsAddDeleteAtom.h"
#include "coupling/tests/TestLammpsInitialVelocity.h"
#include "coupling/tests/TestLammpsGhost.h"
#include "coupling/tests/TestLammpsCalculateForceEnergy.h"
#include <mpi.h>



/** executes a newly created test and deletes it immediately again. */
void runTest(Test *test){
  if (test==NULL){
    std::cout << "ERROR executeTest: test==NULL!" << std::endl; exit(EXIT_FAILURE);
  }
  test->run();
  delete test;
}



int main(int argc, char **argv){
  int size;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  const bool valid2D=(size==1) || (size==4) || (size==16);
  const bool valid3D=(size==1) || (size==8) || (size==64);

  // run tests
  if (valid3D) runTest(new TestLammps<3>(argc,argv,"TestLammps3D"));
  if (valid2D) runTest(new TestLammps<2>(argc,argv,"TestLammps2D"));

  if (valid3D) runTest(new TestLammpsMoleculeIterator<3>(argc,argv,"TestLammpsMoleculeIterator3D"));
  if (valid2D) runTest(new TestLammpsMoleculeIterator<2>(argc,argv,"TestLammpsMoleculeIterator2D"));

  if (valid3D) runTest(new TestLammpsAddDeleteAtom<3>(argc,argv,"TestLammpsAddDeleteAtom3D"));
  if (valid2D) runTest(new TestLammpsAddDeleteAtom<2>(argc,argv,"TestLammpsAddDeleteAtom2D"));

  if (valid3D) runTest(new TestLammpsInitialVelocity<3>(argc,argv,"TestLammpsInitialVelocity3D"));
  if (valid2D) runTest(new TestLammpsInitialVelocity<2>(argc,argv,"TestLammpsInitialVelocity2D"));

  if (valid3D) runTest(new TestLammpsGhost<3>(argc,argv,"TestLammpsGhost3D"));
  if (valid2D) runTest(new TestLammpsGhost<2>(argc,argv,"TestLammpsGhost2D"));

  if (valid3D) runTest(new TestLammpsCalculateForceEnergy<3>(argc,argv,"TestLammpsCalculateForceEnergy3D"));
  if (valid2D) runTest(new TestLammpsCalculateForceEnergy<2>(argc,argv,"TestLammpsCalculateForceEnergy2D"));

  MPI_Finalize();

  return 0;
}

#endif // _MAIN_CPP_
