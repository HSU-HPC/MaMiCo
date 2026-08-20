#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

#include "simplemd/services/MoleculeService.h"

class MoleculeServiceTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(MoleculeServiceTest);
  CPPUNIT_TEST(testMoleculeIdIncrements);
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp() {
    tarch::utils::RandomNumberService::getInstance().init(true);
    const simplemd::services::MolecularPropertiesService molecularPropertiesService(1 /* mass */, 1 /* epsilon */, 1 /* sigma */, 1 /* cutOffRadius */,
                                                                                    1 /*kB */
    );
    tarch::la::Vector<MD_DIM, double> numCells(0);
    tarch::la::Vector<MD_DIM, unsigned int> numMolecules(0);
    for (size_t i = 0; i < MD_DIM; i++) {
      numCells[i] = _numCellsIf3D[i];
      numMolecules[i] = _numMoleculesIf3D[i];
    }
    const tarch::la::Vector<MD_DIM, double> domainOffset(0);
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
    simplemd::services::ParallelTopologyService parallelTopologyService(numCells, domainOffset, meshWidth, numberProcesses, boundary
#if (MD_PARALLEL == MD_YES)
                                                                        ,
                                                                        MPI_COMM_WORLD
#endif
    );

    _moleculeService = new simplemd::services::MoleculeService(parallelTopologyService.getGlobalDomainSize(), parallelTopologyService.getGlobalDomainOffset(),
                                                               numMolecules, 0 /* meanVelocity */, 1 /* kB */, 1 /* temperature */, 10 /* capacityFactor" */,
                                                               molecularPropertiesService, parallelTopologyService);
  }

  void tearDown() {
    if (_moleculeService != nullptr) {
      delete _moleculeService;
      _moleculeService = nullptr;
    }
  }

  void testMoleculeIdIncrements() {
    size_t offset = 1;
    for (size_t i = 0; i < MD_DIM; i++) {
      offset *= _numMoleculesIf3D[i];
    }
    for (size_t i = 0; i < 1000; i++) {
      const auto expected = i + offset;
      const auto actual = _moleculeService->getNextMoleculeID();
      CPPUNIT_ASSERT(actual == expected);
    }
  }

private:
  // use for persistent tests
  const tarch::la::Vector<3, size_t> _numCellsIf3D = {15, 15, 15};
  const tarch::la::Vector<3, unsigned int> _numMoleculesIf3D = {14, 14, 14};
  simplemd::services::MoleculeService* _moleculeService;
};

CPPUNIT_TEST_SUITE_REGISTRATION(MoleculeServiceTest);