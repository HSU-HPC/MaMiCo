#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

#include "coupling/interface/impl/SimpleMD/SimpleMDMolecule.h"

#include "simplemd/Molecule.h"
#include "tarch/la/Vector.h"

#include <sstream>

class SimpleMDMoleculeTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(SimpleMDMoleculeTest);
  CPPUNIT_TEST(testVelocity);
  CPPUNIT_TEST(testForce);
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp() {}
  void tearDown() {}
  void testVelocity() {
    simplemd::Molecule myMolecule;
    coupling::interface::SimpleMDMolecule simpleMolecule(&myMolecule);
    for (int i = -4000; i < 4001; i += 100) {
      for (int j = -4000; j < 4001; j += 100) {
        for (int k = -4000; k < 4001; k += 100) {
          tarch::la::Vector<3, double> velocity{i + 0.5, j + 0.5, k + 0.5};

          simpleMolecule.setVelocity(velocity);
          tarch::la::Vector<3, double> storedVelocity = simpleMolecule.getVelocity();
          CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("velocity assertion x", storedVelocity[0], velocity[0], 1e-6);
          CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("velocity assertion y", storedVelocity[1], velocity[1], 1e-6);
          CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("velocity assertion z", storedVelocity[2], velocity[2], 1e-6);
        }
      }
    }
  }
  void testForce() {
    simplemd::Molecule myMolecule;
    coupling::interface::SimpleMDMolecule simpleMolecule(&myMolecule);
    for (int i = -4000; i < 4001; i += 100) {
      for (int j = -4000; j < 4001; j += 100) {
        for (int k = -4000; k < 4001; k += 100) {
          tarch::la::Vector<3, double> force{i + 0.5, j + 0.5, k + 0.5};

          simpleMolecule.setForce(force);
          tarch::la::Vector<3, double> storedForce = simpleMolecule.getForce();
          CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("force assertion x", storedForce[0], force[0], 1e-6);
          CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("force assertion y", storedForce[1], force[1], 1e-6);
          CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("force assertion z", storedForce[2], force[2], 1e-6);
        }
      }
    }
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(SimpleMDMoleculeTest);