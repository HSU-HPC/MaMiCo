#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include "coupling/solvers/LBCouetteSolver.h"
#include <cstdlib>

#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
#include <mpi.h>
#endif

using namespace coupling::solvers;
using Solver = coupling::interface::PintableMacroSolver;
using State = coupling::interface::PintableMacroSolverState;

/** 
 *  @author Piet Jarmatz
 */
class PintableLBCouetteSolverTest : public CppUnit::TestFixture {

    CPPUNIT_TEST_SUITE(PintableLBCouetteSolverTest);
    CPPUNIT_TEST( testMode );
    CPPUNIT_TEST( testReturnToZero );
    CPPUNIT_TEST( testRunIntoSameState );
    CPPUNIT_TEST_SUITE_END();

public:
  void setUp() {
    _rank = 0;
    _size = 1;
    #if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
        MPI_Comm_size(MPI_COMM_WORLD, &_size);
        MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
    #endif

    coupling::indexing::IndexingService<3>::getInstance().init({14}, {1}, 
      coupling::paralleltopology::XYZ, 3, (unsigned int)_rank);

    F = std::make_unique<LBCouetteSolver>(50, tarch::la::Vector<3, double>{1.5,0,0}, 
        2.14, 2.5, 0.25, 0, "LBCouette", tarch::la::Vector<3, unsigned int>{1,1,1});
    auto supervisor = F->getSupervisor(5,2);
    G = std::unique_ptr<LBCouetteSolver>{dynamic_cast<LBCouetteSolver*>(supervisor.release())};
  }

  void testSetUp() {
    std::unique_ptr<LBCouetteSolver> F2 = std::make_unique<LBCouetteSolver>(50, tarch::la::Vector<3, double>{1.5,0,0}, 
        2.14, 2.5, 0.25, 0, "LBCouette", tarch::la::Vector<3, unsigned int>{1,1,1});

    CPPUNIT_ASSERT( F->getNumberProcesses() == G->getNumberProcesses());
    CPPUNIT_ASSERT( F->getNumberProcesses() == F2->getNumberProcesses());

    CPPUNIT_ASSERT( F->getAvgNumberLBCells() == G->getAvgNumberLBCells());
    CPPUNIT_ASSERT( F->getAvgNumberLBCells() == F2->getAvgNumberLBCells());
  }

  void tearDown() {}

  void testMode(){
    CPPUNIT_ASSERT( F->getMode() == Solver::Mode::coupling);
    CPPUNIT_ASSERT( G->getMode() == Solver::Mode::supervising);
  }

  void testReturnToZero() {
    std::unique_ptr<State> u0 = F->getState();
    F->advance(1.0);
    std::unique_ptr<State> u1 = F->getState();

    F->setState(u0, 0);
    CPPUNIT_ASSERT( *u0 == *(F->getState()) );
    F->advance(1.0);
    CPPUNIT_ASSERT( *u1 == *(F->getState()) );
  }

  void testRunIntoSameState() {
    F->advance(0.5);
    std::unique_ptr<State> u0 = F->getState();

    std::unique_ptr<LBCouetteSolver> F2 = std::make_unique<LBCouetteSolver>(50, tarch::la::Vector<3, double>{1.5,0,0}, 
        2.14, 2.5, 0.25, 0, "LBCouette", tarch::la::Vector<3, unsigned int>{1,1,1});
    F2->setState(u0, 0);

    sameDensityAndVelocity(F, F2);

    CPPUNIT_ASSERT( *u0 == *(F2->getState()) );
    F->advance(1.5);
    F2->advance(1.5);
    CPPUNIT_ASSERT( *(F->getState()) == *(F2->getState()) );

    sameDensityAndVelocity(F, F2);
  }

  void sameDensityAndVelocity(const std::unique_ptr<LBCouetteSolver>& F1, const std::unique_ptr<LBCouetteSolver>& F2){
    if(_rank > 0) return;
    for(int i = 0; i < 50; i++){
        tarch::la::Vector<3, double> pos = getRandomPos();
        auto v1 = F1->getVelocity(pos);
        auto v2 = F2->getVelocity(pos);
        if(v1 != v2){
            std::cout << pos << std::endl;
            std::cout << std::setprecision (17) << v1 << std::endl;
            std::cout << std::setprecision (17) << v2 << std::endl;
        }
        CPPUNIT_ASSERT_EQUAL( v1, v2 );
        double d1 = F1->getDensity(pos);
        double d2 = F2->getDensity(pos);
        if(d1 != d2){
            std::cout << std::setprecision (17) << d1 << std::endl;
            std::cout << std::setprecision (17) << d2 << std::endl;
        }
        CPPUNIT_ASSERT_EQUAL( d1, d2 );
    }
  }

  tarch::la::Vector<3, double> getRandomPos() {
    double a = (double)(rand() % 500) / 10.0;
    double b = (double)(rand() % 500) / 10.0;
    double c = (double)(rand() % 500) / 10.0;
    return {a,b,c};
  }

private:
    int _size, _rank;
    std::unique_ptr<LBCouetteSolver> F, G;

};

CPPUNIT_TEST_SUITE_REGISTRATION(PintableLBCouetteSolverTest);