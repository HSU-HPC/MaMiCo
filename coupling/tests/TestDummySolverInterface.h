// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _TESTS_TESTDUMMYSOLVERINTERFACE_H_
#define _TESTS_TESTDUMMYSOLVERINTERFACE_H_

#include "coupling/tests/Test.h"
#include "tarch/la/Vector.h"
#include "coupling/solvers/DummySolver.h"
#include "coupling/solvers/DummySolverInterface.h"
#include "coupling/solvers/DummySolverInterfaceService.h"
#include "coupling/CouplingMDDefinitions.h"
#include "coupling/datastructures/MacroscopicCell.h"

/** test class for dummy solver interface and dummy solver interface service.
 * Here I test all the functions defined in the dummy
 *	solver interface and interface service for a test scenario, consisting of
 * 18*18*18 mesh and allocate send and receive buffers
 *  and try to measure their size and check if they are performing as desired.
 *  @author Rahul Arora
 */

class TestDummySolverInterface : public Test {
private:
  int _argc;
  char **_argv;

public:
  TestDummySolverInterface(std::string name, int argc, char **argv)
      : Test(name), _argc(argc), _argv(argv) {}
  virtual ~TestDummySolverInterface() {}

  virtual void run() {
    MPI_Init(&_argc, &_argv);

    // Run dummy solver interface test
    std::cout << "Running test for dummy solver interface " << std::endl;
    testDummySolverInterface();

    // Run dummy solver interface service test
    std::cout << "Running test  for dummy solver interface service "
              << std::endl;
    testDummySolverInterfaceService();
    MPI_Finalize();
  }

  virtual void testDummySolverInterface() {
    // Simple test scenario, 12*12*12 mesh and check whether a cell sends to or
    // receivs from MD solver
    DummySolver _dummysolver(12, 12, 12, 0.5);
    tarch::la::Vector<3, unsigned int> num(12);
    DummySolverInterface _dummysolverinterface(num);
    tarch::la::Vector<3, unsigned int> id1(6);
    tarch::la::Vector<3, unsigned int> id2(1);

    bool flag1 =
        _dummysolverinterface.receiveMacroscopicQuantityFromMDSolver(id1);
    bool flag2 = _dummysolverinterface.sendMacroscopicQuantityToMDSolver(id2);

    if (flag1 == true && flag2 == true) {
      std::cout << "The Test was successful" << std::endl;
    } else {
      std::cout << "The Test was not successful" << std::endl;
    }
  }

  virtual void testDummySolverInterfaceService() {
    // Test scenario 18*18*18 mesh
    DummySolver _dummySolver(18, 18, 18, 0.5);

    tarch::la::Vector<3, unsigned int> numberProcesses(1);
    unsigned int rank = 0;
    tarch::la::Vector<3, double> globalMDDomainSize(0.0);
    tarch::la::Vector<3, double> globalMDDomainOffset(0.0);
    tarch::la::Vector<3, double> macroscopicCellSize(0.0);
    tarch::la::Vector<3, unsigned int> linkedCellsPerMacroscopicCell(2);
    for (unsigned int d = 0; d < 3; d++) {
      globalMDDomainSize[d] = 30.0;
      globalMDDomainOffset[d] = 0.0;
      macroscopicCellSize[d] = 1.25 * linkedCellsPerMacroscopicCell[d];
    }

    // Initialize the dummy solver interface
    DummySolverInterfaceService::getInstance().init(
        numberProcesses, rank, globalMDDomainSize, globalMDDomainOffset,
        macroscopicCellSize);

    tarch::la::Vector<3, unsigned int> loop(2);
    unsigned int sendCounter = 0;
    unsigned int recvCounter = 0;

    // Update the sendCounter or recvCounter based on whether the cell sends
    // data to or receives data from MD Solver
    for (loop[2] = 2; loop[2] < 16; loop[2]++) {
      for (loop[1] = 2; loop[1] < 16; loop[1]++) {
        for (loop[0] = 2; loop[0] < 16; loop[0]++) {
          const tarch::la::Vector<3, unsigned int> index =
              coupling::initDimVector<3>(loop);
          const double mass =
              _dummySolver.getDensity(index[0], index[1], index[2]);
          const tarch::la::Vector<3, double> momentum =
              _dummySolver.getVelocity(index[0], index[1], index[2]);
          bool flagsend =
              DummySolverInterfaceService::getInstance().addToSendBuffer(
                  mass, momentum, index);
          if (flagsend == true) {
            sendCounter++;
          }
          double massr = 1.0;
          tarch::la::Vector<3, double> momentumr(1.0);
          bool flagrecv =
              DummySolverInterfaceService::getInstance().getFromReceiveBuffer(
                  massr, momentumr, index);
          if (flagrecv == true) {
            recvCounter++;
          }
        }
      }
    }

    if (sendCounter == 2232 && recvCounter == 512) {
      std::cout << "The test was successful" << std::endl;
    } else {
      std::cout << "The test was not successful" << std::endl;
    }
    DummySolverInterfaceService::getInstance().shutdown();
  }
};
#endif
