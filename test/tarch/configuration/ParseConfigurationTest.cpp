#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include "tarch/configuration/ParseConfiguration.h"

using namespace tarch;
using namespace configuration;

class ParseConfigurationTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(ParseConfigurationTest);
  CPPUNIT_TEST(test);
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp() {
  }

  void tearDown() {
  }

  void test() {
    std::string xmlFile = "<?xml version=\"1.0\"?>\n<mytag\n\tmydouble=\"0.034\"\n\tmyint=\"4\"\n\tmyvecint=\"1;3;1\"\n\tmybool=\"yes\"\n\tmyvecdbl=\"12.2 ; 4.3\">\n</mytag>";
    tinyxml2::XMLDocument file;
    file.Parse(xmlFile.c_str());
    tinyxml2::XMLElement* node = file.FirstChildElement("mytag");
    double mydouble;
    int myint;
    tarch::la::Vector<3, int> myIntVec;
    tarch::la::Vector<2, double> myDblVec;
    bool mybool = false;

    tarch::configuration::ParseConfiguration::readDoubleMandatory(mydouble, node, "mydouble");
    CPPUNIT_ASSERT(mydouble == 0.034);
    tarch::configuration::ParseConfiguration::readIntMandatory(myint, node, "myint");
    CPPUNIT_ASSERT(myint == 4);
    tarch::configuration::ParseConfiguration::readVector<2, double>(myDblVec, node, "myvecdbl");
    tarch::la::Vector<2, double> expectedMyDblVec{12.2, 4.3};
    CPPUNIT_ASSERT(myDblVec == expectedMyDblVec);
    tarch::configuration::ParseConfiguration::readVector<3, int>(myIntVec, node, "myvecint");
    tarch::la::Vector<3, int> expectedMyIntVec{1,3,1};
    CPPUNIT_ASSERT(myIntVec == expectedMyIntVec);
    tarch::configuration::ParseConfiguration::readBoolMandatory(mybool, node, "mybool");
    CPPUNIT_ASSERT(mybool);
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(ParseConfigurationTest);
