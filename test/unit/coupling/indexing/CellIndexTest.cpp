#include "coupling/CouplingMDDefinitions.h"
#include "coupling/indexing/IndexingService.h"
#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

using namespace coupling::indexing;

enum class Testmode { conversion, operators, loop };

class CellIndexTest : public CppUnit::TestFixture {

  CPPUNIT_TEST_SUITE(CellIndexTest);
  CPPUNIT_TEST(testAllConversions);
  CPPUNIT_TEST(testAllOperators);
  CPPUNIT_TEST(testAllLoops);
  CPPUNIT_TEST(testOperatorStream);
  CPPUNIT_TEST_SUITE_END();

public:
  using Trait = IndexTrait;

  using T00 = CellIndex<3>;
  using T01 = CellIndex<3, Trait::vector>;
  using T02 = CellIndex<3, Trait::local>;
  using T03 = CellIndex<3, Trait::vector, Trait::local>;
  using T04 = CellIndex<3, Trait::md2macro>;
  using T05 = CellIndex<3, Trait::vector, Trait::md2macro>;
  using T06 = CellIndex<3, Trait::local, Trait::md2macro>;
  using T07 = CellIndex<3, Trait::vector, Trait::local, Trait::md2macro>;
  using T08 = CellIndex<3, Trait::noGhost>;
  using T09 = CellIndex<3, Trait::vector, Trait::noGhost>;
  using T10 = CellIndex<3, Trait::local, Trait::noGhost>;
  using T11 = CellIndex<3, Trait::vector, Trait::local, Trait::noGhost>;
  using T12 = CellIndex<3, Trait::md2macro, Trait::noGhost>;
  using T13 = CellIndex<3, Trait::vector, Trait::md2macro, Trait::noGhost>;
  using T14 = CellIndex<3, Trait::local, Trait::md2macro, Trait::noGhost>;
  using T15 = CellIndex<3, Trait::vector, Trait::local, Trait::md2macro, Trait::noGhost>;

  void setUp() {
    _rank = 0;
    _size = 1;
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
    MPI_Comm_size(MPI_COMM_WORLD, &_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
#endif

    tarch::la::Vector<3, unsigned int> numberProcesses{(unsigned int)_size, 1, 1};
    if (_size == 4)
      numberProcesses = {2, 2, 1};

    IndexingService<3>::getInstance().init({12}, numberProcesses, coupling::paralleltopology::XYZ, 3, (unsigned int)_rank);
  }

  void tearDown() { IndexingService<3>::getInstance().finalize(); }

  void testAllConversions() {
    _conversionsTested = 0;
    testAll<Testmode::conversion>();
    std::cout << "CellIndexTest::testAllConversions(): #conversions tested = " << _conversionsTested << std::endl;
  }

  void testAllOperators() {
    _conversionsTested = 0;
    testAll<Testmode::operators>();
    std::cout << "CellIndexTest::testAllOperators(): #operations tested = " << _conversionsTested << std::endl;
  }

  void testAllLoops() {
    _conversionsTested = 0;
    testAll<Testmode::loop>();
    std::cout << "CellIndexTest::testAllLoops(): #domains tested = " << _conversionsTested << std::endl;
  }

  template <Testmode mode> void testAll() {
    runtests<mode, T00>();
    runtests<mode, T01>();
    runtests<mode, T02>();
    runtests<mode, T03>();
    runtests<mode, T04>();
    runtests<mode, T05>();
    runtests<mode, T06>();
    runtests<mode, T07>();
    runtests<mode, T08>();
    runtests<mode, T09>();
    runtests<mode, T10>();
    runtests<mode, T11>();
    runtests<mode, T12>();
    runtests<mode, T13>();
    runtests<mode, T14>();
    runtests<mode, T15>();
  }

  template <Testmode mode, class T_in> void runtests() {
    if constexpr (mode == Testmode::conversion) {
      for (const auto& i : T_in()) {
        testConversion<T_in, T00>(i);
        testConversion<T_in, T01>(i);
        testConversion<T_in, T02>(i);
        testConversion<T_in, T03>(i);
        testConversion<T_in, T04>(i);
        testConversion<T_in, T05>(i);
        testConversion<T_in, T06>(i);
        testConversion<T_in, T07>(i);
        testConversion<T_in, T08>(i);
        testConversion<T_in, T09>(i);
        testConversion<T_in, T10>(i);
        testConversion<T_in, T11>(i);
        testConversion<T_in, T12>(i);
        testConversion<T_in, T13>(i);
        testConversion<T_in, T14>(i);
        testConversion<T_in, T15>(i);
      }
    }
    if constexpr (mode == Testmode::operators)
      testOperators<T_in>();
    if constexpr (mode == Testmode::loop)
      testLoop<T_in>();
  }

  template <class T_in, class T_out> void testConversion(T_in idx) {
    if (!T_out::contains(idx))
      return;

    // try "back and forth" conversion
    T_in converted = T_out{idx};
    CPPUNIT_ASSERT_EQUAL(converted, idx);
    _conversionsTested++;
  }

  template <class T> void testOperators() {
    for (const auto& i : T()) {
      T copy = i;
      ++copy;
      if (T::contains(copy)) {
        --copy;
        CPPUNIT_ASSERT_EQUAL(copy, i);
        _conversionsTested++;
      }
      copy = i;
      --copy;
      if (T::contains(copy)) {
        ++copy;
        CPPUNIT_ASSERT_EQUAL(copy, i);
        _conversionsTested++;
      }
    }
  }

  template <class T> void testLoop() {
    tearDown();
    IndexingService<3>::getInstance().init({12}, {(unsigned int)_size, 1, 1}, coupling::paralleltopology::XYZ, 3, (unsigned int)_rank);

    unsigned int count = 0;
    T foo;
    try {
      for (const auto& i : T()) {
        foo = i;
        count++;
      }
    } catch (const std::runtime_error& e) {
      std::cout << _rank << " ERRRROR testLoop() T = " << T::TNAME << std::endl;
      std::cout << _rank << " T::lowerBoundary.get() = " << T::lowerBoundary.get() << std::endl;
      std::cout << _rank << " T::upperBoundary.get() = " << T::upperBoundary.get() << std::endl;
      std::cout << _rank << " T::numberCellsInDomain = " << T::numberCellsInDomain << std::endl;
      std::cout << _rank << " T::linearNumberCellsInDomain = " << T::linearNumberCellsInDomain << std::endl;
      CPPUNIT_FAIL("CellIndexTest::testLoop(): got std::runtime_error!");
    }

    CPPUNIT_ASSERT_EQUAL(count, T::linearNumberCellsInDomain);
    _conversionsTested++;

    tearDown();
    setUp();
  }

  void testOperatorStream() {
    testOperatorsStream<T00, T01>();
    testOperatorsStream<T02, T03>();
    testOperatorsStream<T04, T05>();
    testOperatorsStream<T06, T07>();
    testOperatorsStream<T08, T09>();
    testOperatorsStream<T10, T11>();
    testOperatorsStream<T12, T13>();
    testOperatorsStream<T14, T15>();
  }

  template <class T_scalar, class T_vector> void testOperatorsStream() {
    std::stringstream ss, ss2, ss3;

    // scalar is expected to print the value
    T_scalar a{3};
    ss << a;
    CPPUNIT_ASSERT_EQUAL(ss.str(), std::string("3"));

    // vector is expected to print like tarch::la::Vector
    tarch::la::Vector<3, int> vec{1, 2, 3};
    T_vector b{vec};
    ss2 << b;
    ss3 << vec;
    CPPUNIT_ASSERT_EQUAL(ss2.str(), ss3.str());
  }

private:
  int _size, _rank;
  int _conversionsTested;
};

CPPUNIT_TEST_SUITE_REGISTRATION(CellIndexTest);