#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include "coupling/paralleltopology/ParallelTopologyFactory.h"
#include <typeinfo>

class ParallelTopologyFactoryTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(ParallelTopologyFactoryTest);
  CPPUNIT_TEST(testGetParallelTopology);
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp() {
  }

  void tearDown() {
    delete _topology2;
    delete _topology3;
  }

   void testGetParallelTopology() {
    using namespace coupling::paralleltopology;
    _topology2 = ParallelTopologyFactory::getParallelTopology<2>(ParallelTopologyType::XYZ, tarch::la::Vector<2, unsigned int>{1});
    CPPUNIT_ASSERT(typeid(*_topology2)==typeid(XYZTopology<2>));
    delete _topology2;
    _topology2 = ParallelTopologyFactory::getParallelTopology<2>(ParallelTopologyType::ZYX, tarch::la::Vector<2, unsigned int>{1});
    CPPUNIT_ASSERT(typeid(*_topology2)==typeid(ZYXTopology<2>));
    delete _topology2;
    _topology2 = ParallelTopologyFactory::getParallelTopology<2>(ParallelTopologyType::UNDEFINED, tarch::la::Vector<2, unsigned int>{1});
    CPPUNIT_ASSERT(_topology2==NULL);
    _topology3 = ParallelTopologyFactory::getParallelTopology<3>(ParallelTopologyType::XYZ, tarch::la::Vector<3, unsigned int>{1});
    CPPUNIT_ASSERT(typeid(*_topology3)==typeid(XYZTopology<3>));
    delete _topology3;
    _topology3 = ParallelTopologyFactory::getParallelTopology<3>(ParallelTopologyType::ZYX, tarch::la::Vector<3, unsigned int>{1});
    CPPUNIT_ASSERT(typeid(*_topology3)==typeid(ZYXTopology<3>));
    delete _topology3;
    _topology3 = ParallelTopologyFactory::getParallelTopology<3>(ParallelTopologyType::UNDEFINED, tarch::la::Vector<3, unsigned int>{1});
    CPPUNIT_ASSERT(_topology3==NULL);
  }

private:
  coupling::paralleltopology::ParallelTopology<2>* _topology2;
  coupling::paralleltopology::ParallelTopology<3>* _topology3;

};

CPPUNIT_TEST_SUITE_REGISTRATION(ParallelTopologyFactoryTest);