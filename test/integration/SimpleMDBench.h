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

enum class SimpleMDBenchSize { MD60, MD120, MD240, MD480 };

class BenchSim : public simplemd::MolecularDynamicsSimulation {
public:
  BenchSim(const simplemd::configurations::MolecularDynamicsConfiguration& configuration) : simplemd::MolecularDynamicsSimulation(configuration) {}
  virtual ~BenchSim() {}

  class ChecksumMapping {
  public:
    void beginMoleculeIteration() {}
    void handleMolecule(simplemd::Molecule& molecule) {
      for (unsigned int d = 0; d < MD_DIM; d++) {
        process(molecule.getConstPosition()[d]);
        process(molecule.getConstVelocity()[d]);
        process(molecule.getConstForceOld()[d]);
      }
    }
    void endMoleculeIteration() {}
    unsigned long long checksum() const { return sum; }

    static const bool IsParallel = false;
    static const bool IsReadonly = true;

  private:
    void process(const double& data) { sum ^= *((unsigned long long*)&data); }
    unsigned long long sum = 0;
  };

  bool tarchDebugIsOn() const { return simplemd::services::MoleculeService::tarchDebugIsOn(); }

  unsigned long long getChecksum() {
    ChecksumMapping mapping;
    _moleculeService->getContainer().iterateMolecules(mapping);
    return mapping.checksum();
  }
};

class SimpleMDBench : public Test {
public:
  SimpleMDBench(SimpleMDBenchSize benchsize) :
    Test("SimpleMDBench"), _MDSize(benchsize) {}
  virtual ~SimpleMDBench() {}

  virtual void run() {
    init();
    bench();
    check_result();
    shutdown();
  }

private:
  void init() {
    _rank = 0;
#if (MD_PARALLEL == MD_YES)
    MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
    int size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size != 1) {
      std::cout << "ERROR SimpleMDBench: Must be run with exactly one MPI process!" << std::endl;
      exit(EXIT_FAILURE);
    }
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
    capacity-factor="2.25"

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
    tarch::configuration::ParseConfiguration::parseConfiguration<simplemd::configurations::MolecularDynamicsConfiguration>(fname, "molecular-dynamics",
                                                                                                                           _simpleMDConfig);
    std::remove(fname.c_str());
    if (!_simpleMDConfig.isValid()) {
      std::cout << "ERROR SimpleMDBench: Invalid SimpleMD config!" << std::endl;
      exit(EXIT_FAILURE);
    }
    if(_MDSize == SimpleMDBenchSize::MD120){
      _simpleMDConfig.getDomainConfigurationNonConst().getMoleculesPerDirectionNonConst()
        = tarch::la::Vector<3, unsigned int>{112,112,112};
      _simpleMDConfig.getDomainConfigurationNonConst().getGlobalDomainSizeNonConst()
        = tarch::la::Vector<3, double>{120,120,120};
    } else if(_MDSize == SimpleMDBenchSize::MD240){
      _simpleMDConfig.getDomainConfigurationNonConst().getMoleculesPerDirectionNonConst()
        = tarch::la::Vector<3, unsigned int>{224,224,224};
      _simpleMDConfig.getDomainConfigurationNonConst().getGlobalDomainSizeNonConst()
        = tarch::la::Vector<3, double>{240,240,240};
    } else if(_MDSize == SimpleMDBenchSize::MD480){
      _simpleMDConfig.getDomainConfigurationNonConst().getMoleculesPerDirectionNonConst()
        = tarch::la::Vector<3, unsigned int>{448,448,448};
      _simpleMDConfig.getDomainConfigurationNonConst().getGlobalDomainSizeNonConst()
        = tarch::la::Vector<3, double>{480,480,480};
    }

    _simulation = std::make_unique<BenchSim>(_simpleMDConfig);
    _simulation->initServices();
    unsigned long long sum = _simulation->getChecksum();
    std::cout << "INFO SimpleMDBench: Initial Checksum is " << sum << std::endl;
  }

  void shutdown() { _simulation->shutdownServices(); }

  void bench() {
    // warm-up timestep, for more reliable benchmarking result
    _simulation->simulateOneTimestep(0);
    unsigned long long sum = _simulation->getChecksum();
    std::cout << "INFO SimpleMDBench: Warmup Checksum is " << sum << std::endl;

    timeval start, end;
    gettimeofday(&start, NULL);
    for (unsigned int t = 1; t < _simpleMDConfig.getSimulationConfiguration().getNumberOfTimesteps() + 1; t++) {
      _simulation->simulateOneTimestep(t);
    }
    Kokkos::fence();
    gettimeofday(&end, NULL);
    double runtime = (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec - start.tv_usec);
    std::cout << "Runtime: " << (int)(runtime / 1000) << "ms" << std::endl;
  }

  void check_result() {
#if defined(__FAST_MATH__)
    std::cout << "WARN SimpleMDBench: Result validity check FAILED: -ffast-math is active!" << std::endl;
    return;
#endif

#if (TARCH_DEBUG == TARCH_NO)
    std::cout << "WARN SimpleMDBench: Result validity check FAILED: TARCH_DEBUG is off!" << std::endl;
    return;
#endif

    if (!_simulation->tarchDebugIsOn()) {
      std::cout << "WARN SimpleMDBench: Result validity check FAILED: MoleculeContainer: TARCH_DEBUG is off!" << std::endl;
      return;
    }

    if (!tarch::utils::RandomNumberService::getInstance().tarchDebugIsOn()) {
      std::cout << "WARN SimpleMDBench: Result validity check FAILED: RandomNumberService: TARCH_DEBUG is off!" << std::endl;
      return;
    }

    unsigned long long sum = _simulation->getChecksum();
    std::cout << "INFO SimpleMDBench: Final XOR Checksum is " << sum << std::endl;
    unsigned long long correct;
    if(_MDSize == SimpleMDBenchSize::MD60)  correct = 9224833479527670225u;
    if(_MDSize == SimpleMDBenchSize::MD120) correct = 18377339084083654394u;
    if(_MDSize == SimpleMDBenchSize::MD240) correct = 18342705775936782749u;
    if(_MDSize == SimpleMDBenchSize::MD480) correct = 18349058309621237939u;
    if (sum == correct)
      std::cout << "INFO SimpleMDBench: SUCCESS Checksum is correct :-)" << std::endl;
    else {
      std::cout << "ERROR SimpleMDBench: ERROR Checksum is wrong!! " << std::endl;
      std::cout << "(Note that this is expected and OK in Release build mode,)" << std::endl;
      std::cout << "(as checksum is supposed to match only in Debug / DebugOptimized build mode.)" << std::endl;
      Kokkos::finalize();
      exit(EXIT_FAILURE);
    }
  }

  int _rank;
  simplemd::configurations::MolecularDynamicsConfiguration _simpleMDConfig;
  std::unique_ptr<BenchSim> _simulation;
  SimpleMDBenchSize _MDSize;
};
