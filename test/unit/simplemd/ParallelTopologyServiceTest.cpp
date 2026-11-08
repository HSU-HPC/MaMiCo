#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

#include "simplemd/services/MoleculeService.h"

class ParallelTopologyServiceTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(ParallelTopologyServiceTest);
  CPPUNIT_TEST(testParallelTopologyService);
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp() {
    // TODO remove tarch::utils::RandomNumberService::getInstance().init(true);
    const simplemd::services::MolecularPropertiesService molecularPropertiesService(1 /* mass */, 1 /* epsilon */, 1 /* sigma */, 1 /* cutOffRadius */,
                                                                                    1 /*kB */
    );
    for (size_t i = 0; i < MD_DIM; i++) {
      _numCells[i] = _numCellsIf3D[i];
    }
  }

  void testParallelTopologyService() {
    const tarch::la::Vector<MD_DIM, double> meshWidth(1);
    const tarch::la::Vector<MD_DIM, unsigned int> numberProcesses =
#if (MD_PARALLEL == MD_YES)
#if (MD_DIM > 2)
    { 2,
      2,
      1 }
#else
    { 2,
      2 }
#endif
#else
        {1}
#endif
    ;
    const tarch::la::Vector<MD_LINKED_CELL_NEIGHBOURS, simplemd::BoundaryType> boundary(simplemd::BoundaryType::PERIODIC_BOUNDARY);
    simplemd::services::ParallelTopologyService parallelTopologyService(_numCells, _domainOffset, meshWidth, numberProcesses, boundary
#if (MD_PARALLEL == MD_YES)
                                                                        ,
                                                                        MPI_COMM_WORLD
#endif
    );

    // Trival test for getters
    CPPUNIT_ASSERT(parallelTopologyService.getGlobalDomainSize() == _numCells);
    CPPUNIT_ASSERT(parallelTopologyService.getGlobalDomainOffset() == _domainOffset);
  }

private:
  // use for persistent tests
  const tarch::la::Vector<3, size_t> _numCellsIf3D = {30, 30, 30};
  tarch::la::Vector<MD_DIM, double> _numCells;
  const tarch::la::Vector<MD_DIM, double> _domainOffset = {0, 0, 0};
};

CPPUNIT_TEST_SUITE_REGISTRATION(ParallelTopologyServiceTest);
