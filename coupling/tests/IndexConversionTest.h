// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_TESTS_INDEXCONVERSIONTEST_H_
#define _MOLECULARDYNAMICS_COUPLING_TESTS_INDEXCONVERSIONTEST_H_

#include "coupling/CouplingMDDefinitions.h"
#include "coupling/IndexConversion.h"
#include "coupling/tests/Test.h"

/** tests several properties of the IndexConversion. Currently, only XYZ is
 * tested in this scope.
 *  @author Philipp Neumann
 */
class IndexConversionTest : public Test {
public:
  IndexConversionTest() : Test("IndexConversionTest") {}
  virtual ~IndexConversionTest() {}

  virtual void run() {
    tarch::la::Vector<3, unsigned int> testGlobalNumberCells3D(10, 10, 10);
    tarch::la::Vector<3, unsigned int> testNumberProcesses3D(5, 5, 5);
    tarch::la::Vector<3, double> size3D(1.0);
    tarch::la::Vector<3, double> offset3D(0.0);
    std::cout << "Run 3D test..." << std::endl;
    test3D(testGlobalNumberCells3D, testNumberProcesses3D, size3D, offset3D);
    tarch::la::Vector<2, unsigned int> testGlobalNumberCells2D(10, 10);
    tarch::la::Vector<2, unsigned int> testNumberProcesses2D(5, 5);
    tarch::la::Vector<2, double> size2D(1.0);
    tarch::la::Vector<2, double> offset2D(0.0);
    std::cout << "Run 2D test..." << std::endl;
    test2D(testGlobalNumberCells2D, testNumberProcesses2D, size2D, offset2D);
    std::cout << "Run 3D test unique ranks..." << std::endl;
    testGlobalNumberCells3D = tarch::la::Vector<3, unsigned int>(5);
    testNumberProcesses3D = tarch::la::Vector<3, unsigned int>(2);
    testUniqueRanks<3>(testGlobalNumberCells3D, testNumberProcesses3D, size3D, offset3D);
    std::cout << "Run 2D test unique ranks..." << std::endl;
    testGlobalNumberCells2D = tarch::la::Vector<2, unsigned int>(10);
    testNumberProcesses2D = tarch::la::Vector<2, unsigned int>(3);
    testUniqueRanks<2>(testGlobalNumberCells2D, testNumberProcesses2D, size2D, offset2D);
  }

private:
  void test3D(const tarch::la::Vector<3, unsigned int>& testGlobalNumberCells, const tarch::la::Vector<3, unsigned int>& testNumberProcesses,
              const tarch::la::Vector<3, double>& size, const tarch::la::Vector<3, double>& offset) {
    // loop over all global number cell-configurations
    for (unsigned int zg = 1; zg < testGlobalNumberCells[2]; zg++) {
      for (unsigned int yg = 1; yg < testGlobalNumberCells[1]; yg++) {
        for (unsigned int xg = 1; xg < testGlobalNumberCells[0]; xg++) {
          const tarch::la::Vector<3, unsigned int> globalNumberMacroscopicCells(xg, yg, zg);

          // loop over all number of process configurations
          for (unsigned int znp = 1; znp < testNumberProcesses[2]; znp++) {
            for (unsigned int ynp = 1; ynp < testNumberProcesses[1]; ynp++) {
              for (unsigned int xnp = 1; xnp < testNumberProcesses[0]; xnp++) {
                const tarch::la::Vector<3, unsigned int> numberProcesses(xnp, ynp, znp);
                const unsigned int numberRanks = xnp * ynp * znp;

                // loop over all local process configurations
                for (unsigned int rank = 0; rank < numberRanks; rank++) {
                  coupling::IndexConversion<3> conversion(globalNumberMacroscopicCells, numberProcesses, rank, size, offset, coupling::paralleltopology::XYZ);
                  const tarch::la::Vector<3, unsigned int> wholeLocalDomain =
                      conversion.getLocalNumberMacroscopicCells() + tarch::la::Vector<3, unsigned int>(2);

                  // loop over all cells in the local domain
                  for (unsigned int z = 0; z < wholeLocalDomain[2]; z++) {
                    for (unsigned int y = 0; y < wholeLocalDomain[1]; y++) {
                      for (unsigned int x = 0; x < wholeLocalDomain[0]; x++) {

                        // perform tests
                        const tarch::la::Vector<3, unsigned int> localVectorCellIndex(x, y, z);
                        const unsigned int localCellIndex = conversion.getLocalCellIndex(localVectorCellIndex);
                        if (!conversion.isValidLocalCellIndex(localCellIndex)) {
                          std::cout << "ERROR IndexConversionTest: Local cell index out of range!" << std::endl;
                          std::cout << "Local vector index: " << localVectorCellIndex << ", local cell index: " << localCellIndex << std::endl;
                          return;
                        }

                        // check that conversion from and to vector/scalar index works
                        tarch::la::Vector<3, unsigned int> vec = conversion.getLocalVectorCellIndex(localCellIndex);
                        if ((vec[0] != localVectorCellIndex[0]) || (vec[1] != localVectorCellIndex[1]) || (vec[2] != localVectorCellIndex[2])) {
                          std::cout << "ERROR IndexConversionTest: Conversion local indices failed!" << std::endl;
                          std::cout << "Local vector index: " << localVectorCellIndex << ", local cell index: " << localCellIndex << std::endl;
                          std::cout << "Converted local vector index: " << vec << std::endl;
                          return;
                        }

                        vec = conversion.convertLocalToGlobalVectorCellIndex(vec);
                        if (!conversion.isValidGlobalVectorCellIndex(vec)) {
                          std::cout << "ERROR IndexConversionTest: Global vector cell index out of range!" << std::endl;
                          std::cout << "Local vector cell index: " << localVectorCellIndex << ", global vector cell index: " << vec << std::endl;
                          return;
                        }

                        vec = conversion.convertGlobalToLocalVectorCellIndex(vec);
                        if ((vec[0] != localVectorCellIndex[0]) || (vec[1] != localVectorCellIndex[1]) || (vec[2] != localVectorCellIndex[2])) {
                          std::cout << "ERROR IndexConversionTest: Conversion from global to local indices failed!" << std::endl;
                          std::cout << "Local vector index: " << localVectorCellIndex << ", local cell index: " << localCellIndex << std::endl;
                          std::cout << "Converted local vector index: " << vec << std::endl;
                          return;
                        }
                      }
                    }
                  } // loop over local domain
                }   // loop over local processes
              }
            }
          } // number of processes
        }
      }
    } // global number of macroscopic cells
  }

  void test2D(const tarch::la::Vector<2, unsigned int>& testGlobalNumberCells, const tarch::la::Vector<2, unsigned int>& testNumberProcesses,
              const tarch::la::Vector<2, double>& size, const tarch::la::Vector<2, double>& offset) {
    // loop over all global number cell-configurations
    for (unsigned int yg = 1; yg < testGlobalNumberCells[1]; yg++) {
      for (unsigned int xg = 1; xg < testGlobalNumberCells[0]; xg++) {
        const tarch::la::Vector<2, unsigned int> globalNumberMacroscopicCells(xg, yg);

        // loop over all number of process configurations
        for (unsigned int ynp = 1; ynp < testNumberProcesses[1]; ynp++) {
          for (unsigned int xnp = 1; xnp < testNumberProcesses[0]; xnp++) {
            const tarch::la::Vector<2, unsigned int> numberProcesses(xnp, ynp);
            const unsigned int numberRanks = xnp * ynp;

            // loop over all local process configurations
            for (unsigned int rank = 0; rank < numberRanks; rank++) {
              coupling::IndexConversion<2> conversion(globalNumberMacroscopicCells, numberProcesses, rank, size, offset, coupling::paralleltopology::XYZ);
              const tarch::la::Vector<2, unsigned int> wholeLocalDomain = conversion.getLocalNumberMacroscopicCells() + tarch::la::Vector<2, unsigned int>(2);

              // loop over all cells in the local domain
              for (unsigned int y = 0; y < wholeLocalDomain[1]; y++) {
                for (unsigned int x = 0; x < wholeLocalDomain[0]; x++) {

                  // perform tests
                  const tarch::la::Vector<2, unsigned int> localVectorCellIndex(x, y);
                  const unsigned int localCellIndex = conversion.getLocalCellIndex(localVectorCellIndex);
                  if (!conversion.isValidLocalCellIndex(localCellIndex)) {
                    std::cout << "ERROR IndexConversionTest: Local cell index out of range!" << std::endl;
                    std::cout << "Local vector index: " << localVectorCellIndex << ", local cell index: " << localCellIndex << std::endl;
                    return;
                  }

                  // check that conversion from and to vector/scalar index works
                  tarch::la::Vector<2, unsigned int> vec = conversion.getLocalVectorCellIndex(localCellIndex);
                  if ((vec[0] != localVectorCellIndex[0]) || (vec[1] != localVectorCellIndex[1])) {
                    std::cout << "ERROR IndexConversionTest: Conversion local indices failed!" << std::endl;
                    std::cout << "Local vector index: " << localVectorCellIndex << ", local cell index: " << localCellIndex << std::endl;
                    std::cout << "Converted local vector index: " << vec << std::endl;
                    return;
                  }

                  vec = conversion.convertLocalToGlobalVectorCellIndex(vec);
                  if (!conversion.isValidGlobalVectorCellIndex(vec)) {
                    std::cout << "ERROR IndexConversionTest: Global vector cell index out of range!" << std::endl;
                    std::cout << "Local vector cell index: " << localVectorCellIndex << ", global vector cell index: " << vec << std::endl;
                    return;
                  }

                  vec = conversion.convertGlobalToLocalVectorCellIndex(vec);
                  if ((vec[0] != localVectorCellIndex[0]) || (vec[1] != localVectorCellIndex[1])) {
                    std::cout << "ERROR IndexConversionTest: Conversion from global to local indices failed!" << std::endl;
                    std::cout << "Local vector index: " << localVectorCellIndex << ", local cell index: " << localCellIndex << std::endl;
                    std::cout << "Converted local vector index: " << vec << std::endl;
                    return;
                  }
                }
              } // loop over local domain
            }   // loop over local processes
          }
        } // number of processes
      }
    } // global number of macroscopic cells
  }

  template <unsigned int dim>
  void testUniqueRanks(const tarch::la::Vector<dim, unsigned int>& testGlobalNumberCells, const tarch::la::Vector<dim, unsigned int>& testNumberProcesses,
                       const tarch::la::Vector<dim, double>& size, const tarch::la::Vector<dim, double>& offset) {
    const coupling::IndexConversion<dim> indexConversion(testGlobalNumberCells, testNumberProcesses, 0, size, offset, coupling::paralleltopology::XYZ);
    const tarch::la::Vector<3, unsigned int> end = coupling::initRange<dim>(testGlobalNumberCells + tarch::la::Vector<dim, unsigned int>(2));
    tarch::la::Vector<3, unsigned int> loop(0);
    for (loop[2] = 0; loop[2] < end[2]; loop[2]++) {
      for (loop[1] = 0; loop[1] < end[1]; loop[1]++) {
        for (loop[0] = 0; loop[0] < end[0]; loop[0]++) {
          std::cout << "Unique rank for cell " << coupling::initDimVector<dim>(loop) << ": "
                    << indexConversion.getUniqueRankForMacroscopicCell(coupling::initDimVector<dim>(loop)) << std::endl;
        }
      }
    }
  }
};

#endif // _MOLECULARDYNAMICS_COUPLING_TESTS_INDEXCONVERSIONTEST_H_
