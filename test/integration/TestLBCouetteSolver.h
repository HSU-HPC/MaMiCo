// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_TESTS_TESTLBCOUETTESOLVER_H_
#define _MOLECULARDYNAMICS_COUPLING_TESTS_TESTLBCOUETTESOLVER_H_

#include "Test.h"
#include "coupling/solvers/LBCouetteSolver.h"
#include "coupling/solvers/LBCouetteSolverInterface.h"

class TestLBCouetteSolver : public Test {
public:
  TestLBCouetteSolver() : Test("TestLBCouetteSolver") {}
  virtual ~TestLBCouetteSolver() {}

  virtual void run() {
    tarch::la::Vector<3, unsigned int> numberProcesses(1);
    int rank = 0;
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    numberProcesses = {2, 2, 1}; // tests run at -np=4, hence a decomp of 2x2x1
#endif
    // setup and parameterisation similar to NieTest
    const double channelheight = 50.0;
    const tarch::la::Vector<3, double> wallVelocity(1.5, 0.0, 0.0);
    const double density = 0.81;
    const double kinVisc = 2.14 / density; // dyn. viscosity in NieTest=2.14, density in NieTest=0.81
    const double dx = 2.5;
    const double dt = 0.005 * 100; // MD step in NieTest=0.005
    const int plotEveryTimestep = 1;
    const std::string filestem = "LBCouette";
    const tarch::la::Vector<3, unsigned int> processes(3, 2, 4);

    const int numberTimesteps = 10;

    coupling::solvers::LBCouetteSolver solver(channelheight, wallVelocity, kinVisc, dx, dt, plotEveryTimestep, filestem, processes, 4);

    for (int i = 0; i < numberTimesteps; i++) {
      solver.advance(dt);
      if (rank == 0) {
        std::cout << "Finish time step " << i << "..." << std::endl;
      }
    }
    const tarch::la::Vector<3, double> mdDomainOffset(10.0, 10.0, 2.5);
    const tarch::la::Vector<3, double> mdDomainSize(30.0, 30.0, 30.0);
    const tarch::la::Vector<3, double> couplingCellSize(dx, dx, dx);
    const tarch::la::Vector<3, unsigned int> globalNumberCouplingCells(floor(mdDomainSize[0] / dx + 0.5), floor(mdDomainSize[1] / dx + 0.5),
                                                                       floor(mdDomainSize[2] / dx + 0.5));
    const unsigned int overlapStrip = 2;
    const tarch::la::Vector<3, unsigned int> mdCellOffset(floor(mdDomainOffset[0] / dx + 0.5), floor(mdDomainOffset[1] / dx + 0.5),
                                                          floor(mdDomainOffset[2] / dx + 0.5));
    coupling::solvers::LBCouetteSolverInterface interface(solver.getAvgNumberLBCells(), processes, mdCellOffset, globalNumberCouplingCells, overlapStrip);
    std::vector<coupling::datastructures::CouplingCell<3>*> recvBuffer;
    IDXS.initWithMDSize(interface.getGlobalMDDomainSize(), interface.getGlobalMDDomainOffset(), numberProcesses, couplingCellSize,
                        coupling::paralleltopology::ZYX, 3, (unsigned int)rank);
    unsigned int* globalCellIndices = initRecvBuffer(recvBuffer, interface, density, dx);
    if (rank == 0) {
      std::cout << "Introduce MD domain at offset=" << mdDomainOffset << " and with size " << mdDomainSize << std::endl;
    }

    // test MD boundary and boundary values
    solver.setMDBoundary(mdDomainOffset, mdDomainSize, overlapStrip);
    solver.setMDBoundaryValues(recvBuffer, globalCellIndices);
    for (int i = numberTimesteps; i < 2 * numberTimesteps; i++) {
      solver.advance(dt);
      if (rank == 0) {
        std::cout << "Finish time step " << i << "..." << std::endl;
      }
    }
    outputRanks(interface, solver.getAvgNumberLBCells(), processes, globalNumberCouplingCells, mdCellOffset);
    // delete pointers
    if (globalCellIndices != NULL) {
      delete[] globalCellIndices;
      globalCellIndices = NULL;
    }
    for (unsigned int i = 0; i < recvBuffer.size(); i++) {
      if (recvBuffer[i] != NULL) {
        delete recvBuffer[i];
        recvBuffer[i] = NULL;
      }
    }
    recvBuffer.clear();
  }

private:
  void outputRanks(coupling::solvers::LBCouetteSolverInterface& interface, const tarch::la::Vector<3, unsigned int> avgNumberLBCells,
                   const tarch::la::Vector<3, unsigned int> processes, const tarch::la::Vector<3, unsigned int> globalNumberCouplingCells,
                   const tarch::la::Vector<3, unsigned int> mdCellOffset) {
// only carry out tests for interface on rank 0 (these tests do not depend on parallel parameters)
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank != 0) {
      return;
    }
#endif
    std::cout << "Test LBCouetteSolverInterface..." << std::endl;

    std::cout << "Avg. number LB cells=" << avgNumberLBCells << ", no. processes=" << processes << std::endl;
    std::cout << "MD offset=" << mdCellOffset << ", MD size=" << globalNumberCouplingCells << std::endl;
    for (auto idx : I01()) {
      tarch::la::Vector<3, unsigned int> coords(idx.get());
      std::vector<unsigned int> ranks = interface.getRanks(coords);
      std::cout << "Macro-cell " << coords << " corresponds to LB cell " << coords + mdCellOffset << " and is located on ranks ";
      for (unsigned int i = 0; i < ranks.size(); i++) {
        tarch::la::Vector<3, unsigned int> processCoords;
        processCoords[2] = ranks[i] / (processes[0] * processes[1]);
        processCoords[1] = (ranks[i] - processCoords[2] * processes[0] * processes[1]) / processes[0];
        processCoords[0] = ranks[i] - processCoords[2] * processes[0] * processes[1] - processCoords[1] * processes[0];
        std::cout << ranks[i] << " =" << processCoords << " ; ";
      }
      std::cout << std::endl;
    }

    unsigned int* initRecvBuffer(std::vector<coupling::datastructures::CouplingCell<3>*> & recvBuffer, coupling::solvers::LBCouetteSolverInterface & interface,
                                 const double density, const double dx) {
      // compute avg. mass in this cell
      const double mass = density * dx * dx * dx;
      tarch::la::Vector<3, unsigned int> coords(0);
      // determine local rank
      int rank = 0;
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

      // loop over all global coupling cells
      std::vector<I01> myIndex;
      recvBuffer.clear();
      for (auto idx : I01()) {
        tarch::la::Vector<3, unsigned int> coords(idx.get());
        // check if thi<s cell is located on the current (target) rank
        std::vector<unsigned int> ranks = interface.getTargetRanks(coords);
        bool contained = false;
        for (unsigned int i = 0; i < ranks.size(); i++) {
          contained = contained || (ranks[i] == (unsigned int)rank);
        }

        // if this cell shall be received by the current (target) rank, create a cell in the recvBuffer and store the index in the vector
        if (contained) {
          myIndex.push_back(I01{coords});
          recvBuffer.push_back(new coupling::datastructures::CouplingCell<3>());
          if (recvBuffer[recvBuffer.size() - 1] == NULL) {
            std::cout << "ERROR TestLBCouetteSolver::initRecvBuffer(): recBuffer[..]==NULL!" << std::endl;
            exit(EXIT_FAILURE);
          }
          recvBuffer[recvBuffer.size() - 1]->setMacroscopicMass(mass);
          // set a reference velocity between 0 and one in each component of the velocity vector for testing
          tarch::la::Vector<3, double> localVel(((double)coords[0]) / (globalNumberCouplingCells[0] + 2),
                                                ((double)coords[1]) / (globalNumberCouplingCells[1] + 2),
                                                ((double)coords[2]) / (globalNumberCouplingCells[2] + 2));
          std::cout << "Global cell=" << coords << ", mass=" << mass << ", vel=" << localVel << std::endl;
          recvBuffer[recvBuff > er.size() - 1]->setMacroscopicMomentum(mass * localVel);
        }
      }
      // allocate indices and copy entries from vector
      unsigned int* indices = new unsigned int[recvBuffer.size()];
      if (indices == NULL) {
        std::cout << "ERROR TestLBCouetteSolver::initRecvBuffer(): indices==NULL!" << std::endl;
        exit(EXIT_FAILURE);
      }
      for (unsigned int i = 0; i < recvBuffer.size(); i++) {
        indices[i] = myIndex[i];
      }
      return indices;
    }
  }
};

#endif
