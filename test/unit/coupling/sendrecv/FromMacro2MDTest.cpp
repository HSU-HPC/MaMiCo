#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

#include "coupling/CouplingMDDefinitions.h"
#include "coupling/IndexConversion.h"
#include "coupling/sendrecv/FromMacro2MD.h"
#include "test/unit/coupling/sendrecv/TestCell.h"
#include "test/unit/coupling/sendrecv/TestDataExchangeFromMD2Macro.h"
#include <unistd.h>
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
#include <cmath>
#include <mpi.h>
#endif

/** tests the communication from a selection of cells (macroscopic solver) to a block of cells (coupling tool).
 *  @author Philipp Neumann
 */
class FromMacro2MDTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(FromMacro2MDTest);
  // TODO While this code is named "test", it just prints debug output, but does not execute any tests
  // Thus there is no reason to execute it at all
  // CPPUNIT_TEST(test<2>);
  // CPPUNIT_TEST(test<3>);
  CPPUNIT_TEST_SUITE_END();

public:
  template <unsigned int dim> void test() {
    std::cout << "FromMacro2MDTest: Test " << dim << "D..." << std::endl;
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
    const tarch::la::Vector<dim, unsigned int> globalNumberCouplingCells(5 * oneDir + 2);
    unsigned int numberCellsInclBoundary = globalNumberCouplingCells[0] + 2;
    for (unsigned int d = 1; d < dim; d++) {
      numberCellsInclBoundary = numberCellsInclBoundary * (globalNumberCouplingCells[d] + 2);
    }
    const tarch::la::Vector<dim, unsigned int> numberProcesses(oneDir);
    const tarch::la::Vector<dim, double> mdDomainSize(1.0);
    const tarch::la::Vector<dim, double> mdDomainOffset(0.0);
    // output information
    if (rank == 0) {
      std::cout << "Global number coupling cells: " << globalNumberCouplingCells << std::endl;
      std::cout << "Number processes: " << numberProcesses << std::endl;
    }

    // init new indexing system
    coupling::indexing::IndexingService<dim>::getInstance().init(globalNumberCouplingCells, numberProcesses, coupling::paralleltopology::XYZ, 3,
                                                                 (unsigned int)rank);

    // define functional objects
    coupling::IndexConversion<dim> indexConversion(globalNumberCouplingCells, numberProcesses, rank, mdDomainSize, mdDomainOffset,
                                                   coupling::paralleltopology::XYZ);
    coupling::sendrecv::FromMacro2MD<TestCell<dim>, dim> fromMacro2MD;
    TestDataExchangeFromMacro2MD<dim> testDataExchangeFromMacro2MD(10, indexConversion);

    // initialise the test cells of MaMiCo and set all cell values to 0.0
    std::vector<TestCell<dim>*> cellsMamico;
    for (unsigned int i = 0; i < numberCellsInclBoundary; i++) {
      cellsMamico.push_back(new TestCell<dim>());
      if (cellsMamico[i] == NULL) {
        std::cout << "ERROR FromMD2MacroTest::testSequential3D: cellsMamico[i] ==NULL!" << std::endl;
        exit(EXIT_FAILURE);
      }
      cellsMamico[i]->setBuffer1(tarch::la::Vector<dim, double>(0.0));
      cellsMamico[i]->setBuffer2(0.0);
    }

    // get macroscopic solver buffers and initialise cells before sending
    std::vector<TestCell<dim>*> sentCells;
    unsigned int* sentGlobalIndices = NULL;
    unsigned int numberSentCells = 0;
    testDataExchangeFromMacro2MD.getBuffer4MacroscopicSolverCells(numberSentCells, sentCells, sentGlobalIndices);
    for (unsigned int i = 0; i < numberSentCells; i++) {
      sentCells[i]->setBuffer1(tarch::la::Vector<dim, double>(sentGlobalIndices[i]));
      sentCells[i]->setBuffer2(sentGlobalIndices[i]);
    }
    // INCLUDE THE FOLLOWING LINES FOR BETTER OUTPUT
    sleep(indexConversion.getThisRank());
    std::cout << "Send cells from rank " << indexConversion.getThisRank() << ": " << std::endl;
    for (unsigned int i = 0; i < numberSentCells; i++) {
      std::cout << "Cell " << sentGlobalIndices[i] << "(" << indexConversion.getGlobalVectorCellIndex(sentGlobalIndices[i])
                << "): " << sentCells[i]->getBuffer1() << " , " << sentCells[i]->getBuffer2() << std::endl;
    }
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    // exchange quantities from MD to macro
    fromMacro2MD.sendFromMacro2MD(indexConversion, testDataExchangeFromMacro2MD, cellsMamico, sentCells, sentGlobalIndices);

    // INCLUDE THE FOLLOWING LINES FOR BETTER OUTPUT
    sleep(indexConversion.getThisRank());
    std::cout << "Received information on rank " << indexConversion.getThisRank() << ": " << std::endl;
    for (unsigned int i = 0; i < numberCellsInclBoundary; i++) {
      if (cellsMamico[i]->getBuffer2() > 0.0) {
        std::cout << "Cell " << indexConversion.convertLocalToGlobalCellIndex(i) << "("
                  << indexConversion.getGlobalVectorCellIndex(indexConversion.convertLocalToGlobalCellIndex(i)) << "): " << cellsMamico[i]->getBuffer1()
                  << " , " << cellsMamico[i]->getBuffer2() << std::endl;
      }
    }

    // free memory
    for (unsigned int i = 0; i < (unsigned int)sentCells.size(); i++) {
      if (sentCells[i] != NULL) {
        delete sentCells[i];
        sentCells[i] = NULL;
      }
    }
    if (sentGlobalIndices != NULL) {
      delete[] sentGlobalIndices;
      sentGlobalIndices = NULL;
    }
    for (unsigned int i = 0; i < (unsigned int)cellsMamico.size(); i++) {
      if (cellsMamico[i] != NULL) {
        delete cellsMamico[i];
        cellsMamico[i] = NULL;
      }
    }
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(FromMacro2MDTest);