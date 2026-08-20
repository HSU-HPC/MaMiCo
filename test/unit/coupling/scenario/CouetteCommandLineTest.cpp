#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

#include "coupling/scenario/CouetteCommandLine.h"

class CouetteCommandLineTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(CouetteCommandLineTest);
  CPPUNIT_TEST(testDefaultConfiguration);
  CPPUNIT_TEST(testCustomConfigurationAndRemainingArguments);
  CPPUNIT_TEST(testEqualsSyntax);
  CPPUNIT_TEST(testMissingConfigurationPath);
  CPPUNIT_TEST(testDuplicateConfiguration);
  CPPUNIT_TEST_SUITE_END();

public:
  void testDefaultConfiguration() {
    const CouetteCommandLineOptions options = parseCouetteCommandLine({"couette"});
    CPPUNIT_ASSERT_EQUAL(std::string("couette.xml"), options.configurationFile);
    CPPUNIT_ASSERT_EQUAL(std::size_t{1}, options.remainingArguments.size());
  }

  void testCustomConfigurationAndRemainingArguments() {
    const CouetteCommandLineOptions options = parseCouetteCommandLine({"couette", "--config", "configs/short-run.xml", "--kokkos-threads=2"});
    CPPUNIT_ASSERT_EQUAL(std::string("configs/short-run.xml"), options.configurationFile);
    CPPUNIT_ASSERT_EQUAL(std::size_t{2}, options.remainingArguments.size());
    CPPUNIT_ASSERT_EQUAL(std::string("--kokkos-threads=2"), options.remainingArguments[1]);
  }

  void testEqualsSyntax() {
    const CouetteCommandLineOptions options = parseCouetteCommandLine({"couette", "--config=configs/short-run.xml"});
    CPPUNIT_ASSERT_EQUAL(std::string("configs/short-run.xml"), options.configurationFile);
  }

  void testMissingConfigurationPath() {
    CPPUNIT_ASSERT_THROW(parseCouetteCommandLine({"couette", "--config"}), std::invalid_argument);
  }

  void testDuplicateConfiguration() {
    CPPUNIT_ASSERT_THROW(parseCouetteCommandLine({"couette", "--config", "couette.xml", "--config", "other.xml"}), std::invalid_argument);
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(CouetteCommandLineTest);
