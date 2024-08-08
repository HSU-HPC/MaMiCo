#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/CompilerOutputter.h>
#include <cppunit/TestResult.h>
#include <cppunit/TestResultCollector.h>
#include <cppunit/TextTestProgressListener.h>
#include "coupling/CouplingMDDefinitions.h"
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
#include <mpi.h>
#endif

class Listener : public CppUnit::TextTestProgressListener {
 public:
     void startTest(CppUnit::Test *test) override {
        CppUnit::stdCOut() << test->getName() << " ...\n";
     }
};

int main( int argc, char **argv)
{
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
    MPI_Init(&argc, &argv);
#endif
    std::string testPath = (argc > 1) ? std::string(argv[1]) : "";

    // Create the event manager and test controller
    CppUnit::TestResult controller;

    // Add a listener that colllects test result
    CppUnit::TestResultCollector result;
    controller.addListener( &result );        

    // Add a listener that print dots as test run.
    Listener progress;
    controller.addListener( &progress );      

    // Add the top suite to the test runner
    CppUnit::TestRunner runner;
    runner.addTest( CppUnit::TestFactoryRegistry::getRegistry().makeTest() );   
    try
    {
        std::cout << "Running  ..." << std::endl;
        runner.run( controller, testPath );

        std::cerr << std::endl;

        // Print test in a compiler compatible format.
        CppUnit::CompilerOutputter outputter( &result, std::cerr );
        outputter.write();                      
    }
        catch ( std::invalid_argument &e )  // Test path not resolved
    {
        std::cerr   <<  std::endl  
                    <<  "ERROR: "  <<  e.what()
                    << std::endl;
        return 0;
    }

#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
    MPI_Finalize();
#endif
    return result.wasSuccessful() ? 0 : 1;
}