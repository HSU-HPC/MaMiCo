#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include <sstream>

#include "coupling/CouplingMDDefinitions.h"
#include "coupling/IndexConversion.h"

/** tests several properties of the IndexConversion. Currently, only XYZ is
 * tested in this scope.
 *  @author Philipp Neumann
 */
class IndexConversionTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(IndexConversionTest);
  CPPUNIT_TEST(test3D);
  CPPUNIT_TEST(test2D);
  CPPUNIT_TEST_SUITE_END();

public:
  void test3D() {
    tarch::la::Vector<3, unsigned int> testGlobalNumberCells3D(6, 4, 4);
    tarch::la::Vector<3, unsigned int> testNumberProcesses3D(3, 2, 2);
    tarch::la::Vector<3, double> size3D(1.0);
    tarch::la::Vector<3, double> offset3D(0.0);
    std::cout << "Run 3D test..." << std::endl;
    _test3D(testGlobalNumberCells3D, testNumberProcesses3D, size3D, offset3D);
  }

  void test2D() {
    tarch::la::Vector<2, unsigned int> testGlobalNumberCells2D(10, 6);
    tarch::la::Vector<2, unsigned int> testNumberProcesses2D(5, 3);
    tarch::la::Vector<2, double> size2D(1.0);
    tarch::la::Vector<2, double> offset2D(0.0);
    std::cout << "Run 2D test..." << std::endl;
    _test2D(testGlobalNumberCells2D, testNumberProcesses2D, size2D, offset2D);
  }

private:

  template <unsigned int dim>
  void _perform_tests(const coupling::IndexConversion<dim>& conversion, tarch::la::Vector<dim, unsigned int> localVectorCellIndex, unsigned int localCellIndex){
    if(!conversion.isValidLocalCellIndex(localCellIndex)){
      std::stringstream msg;
      msg << "ERROR IndexConversionTest: Local cell index out of range!" << std::endl;
      msg << "Local vector index: " << localVectorCellIndex << ", local cell index: " << localCellIndex << std::endl;
      CPPUNIT_FAIL(msg.str());
    }
    
    tarch::la::Vector<dim, unsigned int> vec = conversion.getLocalVectorCellIndex(localCellIndex);
    // check that conversion from and to vector/scalar index works
    if(vec != localVectorCellIndex){
      std::stringstream msg;
      msg << "ERROR IndexConversionTest: Conversion local indices failed!" << std::endl;
      msg << "Local vector index: " << localVectorCellIndex << ", local cell index: " << localCellIndex << std::endl;
      msg << "Converted local vector index: " << vec << std::endl;
      CPPUNIT_FAIL(msg.str());
    }

    vec = conversion.convertLocalToGlobalVectorCellIndex(vec);
    if(!conversion.isValidGlobalVectorCellIndex(vec)){
      std::stringstream msg;
      msg << "ERROR IndexConversionTest: Global vector cell index out of range!" << std::endl;
      msg << "Local vector cell index: " << localVectorCellIndex << ", global vector cell index: " << vec << std::endl;
      CPPUNIT_FAIL(msg.str());
    }    

    vec = conversion.convertGlobalToLocalVectorCellIndex(vec);
    if(vec != localVectorCellIndex){
      std::stringstream msg;
      msg << "ERROR IndexConversionTest: Conversion from global to local indices failed!" << std::endl;
      msg << "Local vector index: " << localVectorCellIndex << ", local cell index: " << localCellIndex << std::endl;
      msg << "Converted local vector index: " << vec << std::endl;
      CPPUNIT_FAIL(msg.str());
    }
  }

  void _test3D(const tarch::la::Vector<3, unsigned int>& testGlobalNumberCells, const tarch::la::Vector<3, unsigned int>& testNumberProcesses,
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
                      conversion.getLocalNumberCouplingCells() + tarch::la::Vector<3, unsigned int>(2);

                  // loop over all cells in the local domain
                  for (unsigned int z = 0; z < wholeLocalDomain[2]; z++) {
                    for (unsigned int y = 0; y < wholeLocalDomain[1]; y++) {
                      for (unsigned int x = 0; x < wholeLocalDomain[0]; x++) {

                        // perform tests
                        const tarch::la::Vector<3, unsigned int> localVectorCellIndex(x, y, z);
                        const unsigned int localCellIndex = conversion.getLocalCellIndex(localVectorCellIndex);

                        _perform_tests<3>(conversion, localVectorCellIndex, localCellIndex);
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

  void _test2D(const tarch::la::Vector<2, unsigned int>& testGlobalNumberCells, const tarch::la::Vector<2, unsigned int>& testNumberProcesses,
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
              const tarch::la::Vector<2, unsigned int> wholeLocalDomain = conversion.getLocalNumberCouplingCells() + tarch::la::Vector<2, unsigned int>(2);

              // loop over all cells in the local domain
              for (unsigned int y = 0; y < wholeLocalDomain[1]; y++) {
                for (unsigned int x = 0; x < wholeLocalDomain[0]; x++) {

                  // perform tests
                  const tarch::la::Vector<2, unsigned int> localVectorCellIndex(x, y);
                  const unsigned int localCellIndex = conversion.getLocalCellIndex(localVectorCellIndex);

                  _perform_tests<2>(conversion, localVectorCellIndex, localCellIndex);
                }
              } // loop over local domain
            }   // loop over local processes
          }
        } // number of processes
      }
    } // global number of macroscopic cells
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(IndexConversionTest);