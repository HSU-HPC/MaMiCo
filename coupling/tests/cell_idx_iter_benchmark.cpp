// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#include "coupling/CouplingMDDefinitions.h"
#include "coupling/configurations/MaMiCoConfiguration.h"
#include "coupling/indexing/IndexingService.cpp"
#include "coupling/solvers/CouetteSolverInterface.h"
#include "coupling/tests/Test.h"
#include "simplemd/configurations/MolecularDynamicsConfiguration.h"
#include "tarch/configuration/ParseConfiguration.h"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sys/time.h>
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
#include <mpi.h>
#endif

class CellIdxIterBench : public Test {
public:
  /** @brief simple constructor */
  CellIdxIterBench() : Test("CellIdxIterBench") {}
  /** @brief a dummy destructor */
  virtual ~CellIdxIterBench() {}

  virtual void run() {
    init();
    bench();
    shutdown();
  }

private:
  void init() {
    _rank = 0;
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
    MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
#endif

    std::string fname = std::tmpnam(nullptr);
    std::ofstream file(fname.c_str());
    if (!file.is_open()) {
      std::cout << "ERROR CellIdxIterBench: Could not open file " << fname
                << "!" << std::endl;
      exit(EXIT_FAILURE);
    }
    file << R"mdconf(
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
    fix-seed="no"
   />
  <vtk-configuration filename="Molecules" write-every-timestep="0"/>
  <checkpoint-configuration filename="CheckpointSimpleMD" write-every-timestep="0"/>
  <domain-configuration
    molecules-per-direction="28 ; 28 ; 28"
    domain-size="120.0 ; 120.0 ; 120.0"
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
  /> </molecular-dynamics>)mdconf"
         << std::endl;
    file.close();
    tarch::configuration::ParseConfiguration::parseConfiguration<
        simplemd::configurations::MolecularDynamicsConfiguration>(
        fname, "molecular-dynamics", _simpleMDConfig);
    std::remove(fname.c_str());
    if (!_simpleMDConfig.isValid()) {
      std::cout << "ERROR CellIdxIterBench: Invalid SimpleMD config!"
                << std::endl;
      exit(EXIT_FAILURE);
    }

    fname = std::tmpnam(nullptr);
    file.open(fname.c_str());
    if (!file.is_open()) {
      std::cout << "ERROR CellIdxIterBench: Could not open file " << fname
                << "!" << std::endl;
      exit(EXIT_FAILURE);
    }
    file << R"mamicoconf(
    <mamico>
      <macroscopic-cell-configuration
        cell-size="2.5 ; 2.5 ; 2.5"
        linked-cells-per-macroscopic-cell="1 ; 1 ; 1"
        write-every-microscopic-timestep="0"
        microscopic-filename="MacroscopicCell_micro"
        write-every-macroscopic-timestep="0"
        macroscopic-filename="MacroscopicCell_macro"
      />
      <particle-insertion type="usher" maximum-number-of-iterations="100" maximum-number-of-restarts="500" insert-every-timestep="10" tolerance="0.5" />
      <momentum-insertion type="nie-velocity-imposition" outermost-overlap-layer="2" innermost-overlap-layer="3" />
      <transfer-strategy type="nie-transfer" mass-flux-west="yes" mass-flux-east="yes" mass-flux-north="no" mass-flux-south="no" mass-flux-bottom="no" mass-flux-top="no" shift-by-timesteps="0.5"/>
      <boundary-force type="zhou-boundary-force" west="yes" east="yes" north="yes" south="yes" bottom="yes" top="yes" density="0.81" temperature="1.1" />
      <parallel-topology type="xyz" />
      <thermostat type='outerLayers' number-layers='1' />
    </mamico>
    )mamicoconf"
         << std::endl;
    file.close();
    tarch::configuration::ParseConfiguration::parseConfiguration<
        coupling::configurations::MaMiCoConfiguration<3>>(fname, "mamico",
                                                          _mamicoConfig);
    std::remove(fname.c_str());
    if (!_mamicoConfig.isValid()) {
      std::cout << "ERROR CellIdxIterBench: Invalid MaMiCo config!"
                << std::endl;
      exit(EXIT_FAILURE);
    }

    tarch::la::Vector<3, unsigned int> globalNumberMacroscopicCells(0);
    for (unsigned int d = 0; d < 3; d++)
      globalNumberMacroscopicCells[d] = (unsigned int)floor(
          _simpleMDConfig.getDomainConfiguration().getGlobalDomainSize()[d] /
              _mamicoConfig.getMacroscopicCellConfiguration()
                  .getMacroscopicCellSize()[d] +
          0.5);
    coupling::interface::MacroscopicSolverInterface<3> *couetteSolverInterface =
        new coupling::solvers::CouetteSolverInterface<3>(
            globalNumberMacroscopicCells,
            _mamicoConfig.getMomentumInsertionConfiguration()
                .getInnerOverlap());

    coupling::indexing::IndexingService<3>::getInstance().init(
        _simpleMDConfig, _mamicoConfig, couetteSolverInterface,
        (unsigned int)_rank);
    delete couetteSolverInterface;
  }

  void shutdown() {}

  void bench() {
    using namespace coupling::indexing;
    using CellIndex_T = CellIndex<3, IndexTrait::local, IndexTrait::md2macro,
                                  IndexTrait::noGhost>;
    using CellIndex_T_vec =
        CellIndex<3, IndexTrait::vector, IndexTrait::local,
                  IndexTrait::md2macro, IndexTrait::noGhost>;
    unsigned int numcells = CellIndex_T::linearNumberCellsInDomain;
    std::cout << "Number cells in test domain: " << numcells << std::endl;
    std::cout << "lowerBoundary = " << CellIndex_T::lowerBoundary.get()
              << std::endl;
    std::cout << "upperBoundary = " << CellIndex_T::upperBoundary.get()
              << std::endl;
    const int numcounts = 10000;
    timeval start, end;

    std::cout << std::endl << "Scalar benchmark ------------- " << std::endl;
    tarch::la::Vector<3, long> result(0);
    for (int count = 0; count < numcounts; count++)
      for (unsigned int i = 0; i < numcells; i++) {
        result[0] += i; // do something with i so that the compiler can not
                        // optimize this away ...
      }
    // (second try after warm-up phase)
    gettimeofday(&start, NULL);
    result = 0;
    for (int count = 0; count < numcounts; count++)
      for (unsigned int i = 0; i < numcells; i++) {
        result[0] += i; // do something with i so that the compiler can not
                        // optimize this away ...
      }
    std::cout << "Useless result: " << result[0] << std::endl;
    gettimeofday(&end, NULL);
    double runtime =
        (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec - start.tv_usec);
    std::cout << "Raw loop: " << (int)(runtime / 1000) << "ms" << std::endl;

    gettimeofday(&start, NULL);
    result = 0;
    for (int count = 0; count < numcounts; count++)
      for (auto idx : CellIndex_T()) {
        // do something with idx so that the compiler can not optimize this away
        // (and to check if we do the same here)
        result[0] += idx.get();
      }
    std::cout << "Useless result: " << result[0] << std::endl;
    gettimeofday(&end, NULL);
    runtime =
        (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec - start.tv_usec);
    std::cout << "Index range iterator: " << (int)(runtime / 1000) << "ms"
              << std::endl;

    std::cout << std::endl << "Vector benchmark ------------- " << std::endl;
    gettimeofday(&start, NULL);
    result = 0;
    for (int count = 0; count < numcounts; count++)
      for (unsigned int z = 0; z < CellIndex_T::numberCellsInDomain[2]; z++)
        for (unsigned int y = 0; y < CellIndex_T::numberCellsInDomain[1]; y++)
          for (unsigned int x = 0; x < CellIndex_T::numberCellsInDomain[0];
               x++) {
            tarch::la::Vector<3, int> idx(x, y, z);
            result[0] += idx[0];
            result[1] += idx[1];
            result[2] += idx[2];
          }
    std::cout << "Useless result: " << result << std::endl;
    gettimeofday(&end, NULL);
    runtime =
        (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec - start.tv_usec);
    std::cout << "Raw loop: " << (int)(runtime / 1000) << "ms" << std::endl;

    gettimeofday(&start, NULL);
    result = 0;
    for (int count = 0; count < numcounts; count++)
      for (auto idx : CellIndex_T_vec()) {
        result[0] += idx.get()[0];
        result[1] += idx.get()[1];
        result[2] += idx.get()[2];
      }
    std::cout << "Useless result: " << result << std::endl;
    gettimeofday(&end, NULL);
    runtime =
        (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec - start.tv_usec);
    std::cout << "Index range iterator: " << (int)(runtime / 1000) << "ms"
              << std::endl;
  }

  int _rank;
  simplemd::configurations::MolecularDynamicsConfiguration _simpleMDConfig;
  coupling::configurations::MaMiCoConfiguration<3> _mamicoConfig;
};

void runTest(Test *test) {
  if (test == NULL) {
    std::cout << "ERROR executeTest: test==NULL!" << std::endl;
    exit(EXIT_FAILURE);
  }
  test->run();
  delete test;
}

int main(int argc, char *argv[]) {
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
  MPI_Init(&argc, &argv);
#endif

  // run tests
  runTest(new CellIdxIterBench());

#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
  MPI_Finalize();
#endif

  return 0;
};

/* Sampe Output:

Run CellIdxIterBench...
Number cells in test domain: 74088
lowerBoundary = 4 , 4 , 4
upperBoundary = 45 , 45 , 45

Scalar benchmark -------------
Useless result: 27444788280000
Raw loop: 191ms
Useless result: 27444788280000
Index range iterator: 119ms

Vector benchmark -------------
Useless result: 15188040000 , 15188040000 , 15188040000
Raw loop: 147ms
Useless result: 15188040000 , 15188040000 , 15188040000
Index range iterator: 152ms
Shut down CellIdxIterBench

*/