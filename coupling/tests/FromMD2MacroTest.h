// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_TESTS_FROMMD2MACROTEST_H_
#define _MOLECULARDYNAMICS_COUPLING_TESTS_FROMMD2MACROTEST_H_

#include "coupling/CouplingMDDefinitions.h"
#include "coupling/IndexConversion.h"
#include "coupling/sendrecv/FromMD2Macro.h"
#include "coupling/tests/Test.h"
#include "coupling/tests/TestCell.h"
#include "coupling/tests/TestDataExchangeFromMD2Macro.h"
#include <cmath>
#include <unistd.h>
#include <vector>
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
#include <mpi.h>
#endif

/** tests the communication from a block of macroscopic cells (such as the block
 * of MaMiCo) to a selection of cells of a possible macroscopic solver.
 *  @author Philipp Neumann
 */
class FromMD2MacroTest : public Test {
private:
  template <unsigned int dim> void test() {
    // define current rank, total size (only allowed to be cubic)
    int rank = 0;
    int size = 1;
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    const unsigned int oneDir = (unsigned int)floor(pow(((double)size), (1.0 / ((double)dim))) + 0.5);
    if (pow(oneDir, dim) != size) {
      std::cout << "FromMD2MacroTest::test: this is not a cubic domain decomposition!" << std::endl;
      std::cout << "Size=" << size << ", in each direction: " << oneDir << std::endl;
      std::cout << dim << "D case is not carried out." << std::endl;
      return;
    }

    // define domain sizes for cells and MD
    const tarch::la::Vector<dim, unsigned int> globalNumberMacroscopicCells(5 * oneDir + 2);
    unsigned int numberCellsInclBoundary = globalNumberMacroscopicCells[0] + 2;
    for (unsigned int d = 1; d < dim; d++) {
      numberCellsInclBoundary = numberCellsInclBoundary * (globalNumberMacroscopicCells[d] + 2);
    }
    const tarch::la::Vector<dim, unsigned int> numberProcesses(oneDir);
    const tarch::la::Vector<dim, double> mdDomainSize(1.0);
    const tarch::la::Vector<dim, double> mdDomainOffset(0.0);
    // output information
    if (rank == 0) {
      std::cout << "Global number macroscopic cells: " << globalNumberMacroscopicCells << std::endl;
      std::cout << "Number processes: " << numberProcesses << std::endl;
    }

    // define functional objects
    coupling::IndexConversion<dim> indexConversion(globalNumberMacroscopicCells, numberProcesses, rank, mdDomainSize, mdDomainOffset,
                                                   coupling::paralleltopology::XYZ);
    coupling::sendrecv::FromMD2Macro<TestCell<dim>, dim> fromMD2Macro;
    TestDataExchangeFromMD2Macro<dim> testDataExchangeFromMD2Macro(10, indexConversion);

    // initialise the test cells of MaMiCo
    std::cout << "Init test cells of mamico on rank " << rank << std::endl;
    std::vector<TestCell<dim>*> cellsMamico;
    for (unsigned int i = 0; i < numberCellsInclBoundary; i++) {
      cellsMamico.push_back(new TestCell<dim>());
      if (cellsMamico[i] == NULL) {
        std::cout << "ERROR FromMD2Macro::test(): cellsMamico[i]==NULL!" << std::endl;
        exit(EXIT_FAILURE);
      }
      cellsMamico[i]->setBuffer1(tarch::la::Vector<dim, double>(i));
      cellsMamico[i]->setBuffer2(i);
    }

    // get receive buffers
    std::cout << "Get receive buffers on rank " << rank << std::endl;
    std::vector<TestCell<dim>*> receivedCells;
    unsigned int* receivedGlobalIndices = NULL;
    unsigned int numberReceivedCells = 0;
    testDataExchangeFromMD2Macro.getBuffer4MacroscopicSolverCells(numberReceivedCells, receivedCells, receivedGlobalIndices);
    // INCLUDE THE FOLLOWING LINES FOR BETTER OUTPUT
    sleep(indexConversion.getThisRank());
    std::cout << "Init receive cells on rank " << indexConversion.getThisRank() << ": " << std::endl;
    for (unsigned int i = 0; i < numberReceivedCells; i++) {
      std::cout << "Cell " << receivedGlobalIndices[i] << std::endl;
    }
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    // exchange quantities from MD to macro
    fromMD2Macro.sendFromMD2Macro(indexConversion, testDataExchangeFromMD2Macro, cellsMamico, receivedCells, receivedGlobalIndices);

    // INCLUDE THE FOLLOWING LINES FOR BETTER OUTPUT
    sleep(indexConversion.getThisRank());
    std::cout << "Received information on rank " << indexConversion.getThisRank() << ": " << std::endl;
    for (unsigned int i = 0; i < numberReceivedCells; i++) {
      std::cout << "Cell " << receivedGlobalIndices[i] << ": " << receivedCells[i]->getBuffer1() << " , " << receivedCells[i]->getBuffer2() << std::endl;
    }

    // free memory
    for (unsigned int i = 0; i < (unsigned int)receivedCells.size(); i++) {
      if (receivedCells[i] != NULL) {
        delete receivedCells[i];
        receivedCells[i] = NULL;
      }
    }
    if (receivedGlobalIndices != NULL) {
      delete[] receivedGlobalIndices;
      receivedGlobalIndices = NULL;
    }
    for (unsigned int i = 0; i < (unsigned int)cellsMamico.size(); i++) {
      if (cellsMamico[i] != NULL) {
        delete cellsMamico[i];
        cellsMamico[i] = NULL;
      }
    }
  }

public:
  FromMD2MacroTest() : Test("FromMD2MacroTest") {}
  virtual ~FromMD2MacroTest() {}

  virtual void run() {
    std::cout << "FromMD2MacroTest: Test 2D..." << std::endl;
    test<2>();
    std::cout << "FromMD2MacroTest: Test 3D..." << std::endl;
    test<3>();
  }
};
#endif // _MOLECULARDYNAMICS_COUPLING_TESTS_FROMMD2MACROTEST_H_
