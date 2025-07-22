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
  void init() {
    _rank = 0;
#if (MD_PARALLEL == MD_PARALLEL)
    MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
#endif

    std::string fname = "mdconf.xml.tmp." + std::to_string(_rank);
    std::ofstream file(fname.c_str());
    if (!file.is_open()) {
      std::cout << "ERROR SimpleMDBench: Could not open file " << fname << "!" << std::endl;
      exit(EXIT_FAILURE);
    }
    file << R"mdconf(
  <scenario-configuration>
  <molecular-dynamics>
  <molecule-configuration
    mass="1.0"
    temperature="1.27"
    sigma="1.0"
    epsilon="1.0"
    mean-velocity="0.0 ; 0.0 ; 0.0"
  />
  <mpi-configuration number-of-processes="1 ; 1 ; 1" />
  <simulation-configuration
    dt="0.005"
    number-of-timesteps="50"
    reorganise-memory-every-timestep="20"
    compute-macroscopic-quantities-every-timestep="0"
    fix-seed="yes"
   />
  <vtk-configuration filename="Molecules" write-every-timestep="0"/>
  <checkpoint-configuration filename="CheckpointSimpleMD" write-every-timestep="0"/>
  <domain-configuration
    molecules-per-direction="56 ; 56 ; 56"
    domain-size="60.0 ; 60.0 ; 60.0"
    domain-offset="10.0 ; 10.0 ; 2.5"
    cutoff-radius="2.2"
    linked-cell-size="2.5 ; 2.5 ; 2.5"
    k_B="1.0"
    block-size="100"

    bottom-south-west="reflecting" bottom-south="reflecting" bottom-south-east="reflecting"
    bottom-west="reflecting"       bottom="reflecting"       bottom-east="reflecting"
    bottom-north-west="reflecting" bottom-north="reflecting" bottom-north-east="reflecting"
    south-west="reflecting"        south="reflecting"        south-east="reflecting"
    west="reflecting"                                          east="reflecting"
    north-west="reflecting"        north="reflecting"        north-east="reflecting"
    top-south-west="reflecting"    top-south="reflecting"    top-south-east="reflecting"
    top-west="reflecting"          top="reflecting"          top-east="reflecting"
    top-north-west="reflecting"    top-north="reflecting"    top-north-east="reflecting"
  />
  </molecular-dynamics>
  </scenario-configuration>)mdconf"
         << std::endl;
    file.close();
  }

  void shutdown() {}

  void bench() {}

  int _rank;
  simplemd::configurations::MolecularDynamicsConfiguration _simpleMDConfig;
};
