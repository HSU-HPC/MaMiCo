#include "simplemd/MolecularDynamicsDefinitions.h"
#if (MD_PARALLEL == MD_YES)
#include <mpi.h>
#endif
#include "simplemd/MolecularDynamicsSimulation.h"
#include "simplemd/configurations/MolecularDynamicsConfiguration.h"
#include "tarch/configuration/ParseConfiguration.h"
#include "test/integration/Test.h"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sys/time.h>

class SimpleMDBench : public Test {
public:
  SimpleMDBench() : Test("SimpleMDBench") {}
  virtual ~SimpleMDBench() {}

  virtual void run() {
    init();
    bench();
    shutdown();
  }

private:
  void init() {}
  void shutdown() {}
  void bench() {}
};
