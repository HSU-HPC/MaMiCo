#include "tarch/configuration/ParseConfiguration.h"
#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

using namespace tarch;
using namespace configuration;

class ParseConfigurationTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(ParseConfigurationTest);
  CPPUNIT_TEST(testParseXML);
  CPPUNIT_TEST(testParseConfigurationWithRoot);
  CPPUNIT_TEST(testParseConfigurationWithoutRoot);
  CPPUNIT_TEST_SUITE_END();

  const char* const filenameConfigWithRoot = "/tmp/mamicotest_ParseConfigurationTest_testParseConfigurationWithRoot";
  const char* const filenameConfigWithoutRoot = "/tmp/mamicotest_ParseConfigurationTest_testParseConfigurationWithoutRoot";

public:
  void setUp() {
    std::ofstream(filenameConfigWithRoot)
        << "<?xml version=\"1.0\"?>\n<mamico-configuration>\n\t<foo></foo>\n\t<bar></bar>\n</mamico-configuration>";
    std::ofstream(filenameConfigWithoutRoot) << "<?xml version=\"1.0\"?>\n<foo></foo>\n<bar></bar>";
  }

  void tearDown() {
    std::remove(filenameConfigWithRoot);
    std::remove(filenameConfigWithoutRoot);
  }

  void testParseXML() {
    std::string xmlFile =
        "<?xml version=\"1.0\"?>\n<mytag\n\tmydouble=\"0.034\"\n\tmyint=\"4\"\n\tmyvecint=\"1;3;1\"\n\tmybool=\"yes\"\n\tmyvecdbl=\"12.2 ; 4.3\">\n</mytag>";
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
    tarch::configuration::ParseConfiguration::readVectorMandatory<2, double>(myDblVec, node, "myvecdbl");
    tarch::la::Vector<2, double> expectedMyDblVec{12.2, 4.3};
    CPPUNIT_ASSERT(myDblVec == expectedMyDblVec);
    tarch::configuration::ParseConfiguration::readVectorMandatory<3, int>(myIntVec, node, "myvecint");
    tarch::la::Vector<3, int> expectedMyIntVec{1, 3, 1};
    CPPUNIT_ASSERT(myIntVec == expectedMyIntVec);
    tarch::configuration::ParseConfiguration::readBoolMandatory(mybool, node, "mybool");
    CPPUNIT_ASSERT(mybool);
  }

  void testParseConfigurationWithRoot() {
    auto cfg = ParseConfiguration::XMLConfiguration::load(filenameConfigWithRoot);
    CPPUNIT_ASSERT(cfg.error == tinyxml2::XML_NO_ERROR);
    CPPUNIT_ASSERT(cfg.root->FirstChildElement("foo") != NULL);
  }

  void testParseConfigurationWithoutRoot() {
    auto cfg = ParseConfiguration::XMLConfiguration::load(filenameConfigWithoutRoot);
    CPPUNIT_ASSERT(cfg.error == tinyxml2::XML_NO_ERROR);
    CPPUNIT_ASSERT(cfg.root->FirstChildElement("foo") != NULL);
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(ParseConfigurationTest);
