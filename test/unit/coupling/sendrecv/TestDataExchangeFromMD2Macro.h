// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_TESTS_TESTDATAEXCHANGEFROMMD2MACRO_H_
#define _MOLECULARDYNAMICS_COUPLING_TESTS_TESTDATAEXCHANGEFROMMD2MACRO_H_

#include "coupling/CouplingMDDefinitions.h"
#include "coupling/IndexConversion.h"
#include "coupling/sendrecv/DataExchange.h"
#include "test/unit/coupling/sendrecv/TestCell.h"
#include <vector>
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
#include <mpi.h>
#endif

/** this class is used to handle test data exchange operations between dummy MaMiCo and a dummy macroscopic solver.
 *  For N processes, this class sends data from process n to n+1. The only data that are sent
 *  correspond to the points on the main "diagonal" of a field, that is cell data for cells at coordinates (n,n,n) (in global ordering).
 *  @author Philipp Neumann
 */
template <unsigned int dim> class TestDataExchangeFromMD2Macro : public coupling::sendrecv::DataExchange<TestCell<dim>, dim> {

public:
  TestDataExchangeFromMD2Macro(unsigned int tag, const coupling::IndexConversion<dim>& indexConversion)
      : coupling::sendrecv::DataExchange<TestCell<dim>, dim>(tag), _indexConversion(indexConversion) {}
  virtual ~TestDataExchangeFromMD2Macro() {}

  /** returns the ranks to which a particular cell (at index globalCellIndex) should be sent. */
  virtual std::vector<unsigned int> getTargetRanks(tarch::la::Vector<dim, unsigned int> globalCellIndex) {
    bool isSource = true;
    std::vector<unsigned int> vec;
    // empty target ranks if this is a global cell at the lower or upper end of the diagonal (inside ghost layer)
    if ((globalCellIndex == tarch::la::Vector<dim, unsigned int>(0)) ||
        (globalCellIndex == _indexConversion.getGlobalNumberMacroscopicCells() + tarch::la::Vector<dim, unsigned int>(1))) {
      return vec;
    }
    const unsigned int rank = _indexConversion.getUniqueRankForMacroscopicCell(globalCellIndex);
    int numberProcesses = 1;
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
    MPI_Comm_size(MPI_COMM_WORLD, &numberProcesses);
#endif

    for (unsigned int d = 1; d < dim; d++) {
      isSource = isSource && (globalCellIndex[d] == globalCellIndex[0]);
    }

    // if diagonal AND this is not the last rank, add rank+1 as target
    if (isSource) {
      vec.push_back((rank + 1) % numberProcesses);
    }
    return vec;
  }

  /** we obtain this value, if the entry is located on the "diagonal", i.e. the global coordinates correspond to
   *  (n,n,n).
   */
  virtual std::vector<unsigned int> getSourceRanks(tarch::la::Vector<dim, unsigned int> globalCellIndex) {
    bool isSource = true;
    std::vector<unsigned int> vec;
    // empty source ranks if this is a global cell at the lower or upper end of the diagonal (inside ghost layer)
    if ((globalCellIndex == tarch::la::Vector<dim, unsigned int>(0)) ||
        (globalCellIndex == _indexConversion.getGlobalNumberMacroscopicCells() + tarch::la::Vector<dim, unsigned int>(1))) {
      return vec;
    }

    const unsigned int rank = _indexConversion.getUniqueRankForMacroscopicCell(globalCellIndex);

    for (unsigned int d = 1; d < dim; d++) {
      isSource = isSource && (globalCellIndex[d] == globalCellIndex[0]);
    }

    // source is rank
    if (isSource) {
      vec.push_back(rank);
    }
    return vec;
  }

  /** local rule to read from macroscopic cell and write data to (e.g. send) buffer */
  virtual void readFromCell(double* const buffer, const TestCell<dim>& cell) {
    tarch::la::Vector<dim, double> b1 = cell.getBuffer1();
    for (unsigned int d = 0; d < dim; d++) {
      buffer[d] = b1[d];
    }
    buffer[dim] = cell.getBuffer2();
  }

  /** local rule to read from receive buffer and write data to macroscopic cell */
  virtual void writeToCell(const double* const buffer, TestCell<dim>& cell) {
    tarch::la::Vector<dim, double> b1(0.0);
    for (unsigned int d = 0; d < dim; d++) {
      b1[d] = buffer[d];
    }
    cell.setBuffer1(b1);
    cell.setBuffer2(buffer[dim]);
  }

  /** returns the number of doubles that are sent per macroscopic cell. */
  virtual unsigned int getDoublesPerCell() const { return dim + 1; }

  /** initialises the buffers required on the processes which by incrementing 1 are obtained from process coordinates (n,n,n).
   *
   */
  virtual void getBuffer4MacroscopicSolverCells(unsigned int& numberCells, std::vector<TestCell<dim>*>& cells, unsigned int*& globalCellIndices) {
    // for rank 0: expect data from last rank
    if (_indexConversion.getThisRank() == 0) {
      // determine number of cells on last rank (this only works for cubic sub-domains) and start coordinate of first
      // (non-ghost) global cell of last rank
      const tarch::la::Vector<dim, unsigned int> start =
          tarch::la::Vector<dim, unsigned int>(1) +
          (_indexConversion.getNumberProcesses() - tarch::la::Vector<dim, unsigned int>(1)) * _indexConversion.getAverageLocalNumberMacroscopicCells()[0];

      numberCells = (_indexConversion.getGlobalNumberMacroscopicCells() - (_indexConversion.getNumberProcesses() - tarch::la::Vector<dim, unsigned int>(1)) *
                                                                              _indexConversion.getAverageLocalNumberMacroscopicCells()[0])[0];
      for (unsigned int i = 0; i < numberCells; i++) {
        cells.push_back(new TestCell<dim>());
        if (cells[i] == NULL) {
          std::cout << "ERROR getBuffer4MacroscopicSolverCells: cells[i]==NULL!" << std::endl;
          exit(EXIT_FAILURE);
        }
      }
      globalCellIndices = new unsigned int[numberCells];
      if (globalCellIndices == NULL) {
        std::cout << "ERROR TestDataExchangeFromMD2Macro::getBuffer4ReceivedCells: NULL ptr returned!" << std::endl;
        exit(EXIT_FAILURE);
      }

      for (unsigned int i = 0; i < numberCells; i++) {
        globalCellIndices[i] = _indexConversion.getGlobalCellIndex(start + i * tarch::la::Vector<dim, unsigned int>(1));
      }
      // all other cases (rank != 0)
    } else {

      // determine if we obtain information from the neighbouring rank
      const unsigned int neighbourRank = _indexConversion.getThisRank() - 1;
      const tarch::la::Vector<dim, unsigned int> neighbourRankCoords = _indexConversion.getProcessCoordinates(neighbourRank);
      bool isSource = true;
      for (unsigned int d = 1; d < dim; d++) {
        isSource = isSource && (neighbourRankCoords[0] == neighbourRankCoords[d]);
      }

      if (isSource) {
        // we only plan with cubic sub-domains on the processes...
        const unsigned int numberLocalCells = _indexConversion.getAverageLocalNumberMacroscopicCells()[0];

        globalCellIndices = new unsigned int[numberLocalCells];
        for (unsigned int i = 0; i < numberLocalCells; i++) {
          cells.push_back(new TestCell<dim>());
          if (cells[i] == NULL) {
            std::cout << "ERROR getBuffer4MacroscopicSolverCells: cells[i]==NULL!" << std::endl;
            exit(EXIT_FAILURE);
          }
        }
        if (globalCellIndices == NULL) {
          std::cout << "ERROR TestDataExchangeFromMD2Macro::getBuffer4ReceivedCells: NULL ptrs returned!" << std::endl;
          exit(EXIT_FAILURE);
        }
        numberCells = numberLocalCells;

        // determine start coordinates of neighbouring local cells
        tarch::la::Vector<dim, unsigned int> start(1);
        start = start + neighbourRankCoords * numberLocalCells;

        for (unsigned int i = 0; i < numberLocalCells; i++) {
          globalCellIndices[i] = _indexConversion.getGlobalCellIndex(start + i * tarch::la::Vector<dim, unsigned int>(1));
        }
      }
    }
  }

protected:
  const coupling::IndexConversion<dim>& _indexConversion;
};

/** same as TestDataExchangeFromMD2Macro, but with inverted target and source ranks, i.e. the original target ranks are
 *  the current source ranks. Since we hold certain macroscopic cells on several Mamico-blocks (due to the ghost cells),
 *  we use a modified target ranks-vector (based on the indexConversion).
 */
template <unsigned int dim> class TestDataExchangeFromMacro2MD : public TestDataExchangeFromMD2Macro<dim> {
public:
  TestDataExchangeFromMacro2MD(unsigned int tag, const coupling::IndexConversion<dim>& indexConversion)
      : TestDataExchangeFromMD2Macro<dim>(tag, indexConversion) {}
  virtual ~TestDataExchangeFromMacro2MD() {}

  virtual std::vector<unsigned int> getTargetRanks(tarch::la::Vector<dim, unsigned int> globalCellIndex) {
    return TestDataExchangeFromMD2Macro<dim>::_indexConversion.getRanksForMacroscopicCell(globalCellIndex);
  }
  virtual std::vector<unsigned int> getSourceRanks(tarch::la::Vector<dim, unsigned int> globalCellIndex) {
    return TestDataExchangeFromMD2Macro<dim>::getTargetRanks(globalCellIndex);
  }
};

#endif // _MOLECULARDYNAMICS_COUPLING_TESTS_TESTDATAEXCHANGEFROMMD2MACRO_H_
