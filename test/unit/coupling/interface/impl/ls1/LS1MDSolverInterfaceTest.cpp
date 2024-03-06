#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

#include "ls1/src/molecules/Molecule.h"
#include "ls1/src/particleContainer/ParticleContainer.h"
#include "ls1/src/particleContainer/RegionParticleIterator.h"

#include "coupling/interface/impl/ls1/LS1MDSolverInterface.h"
#include "coupling/interface/impl/ls1/LS1MoleculeIterator.h"
#include "coupling/interface/impl/ls1/LS1RegionWrapper.h"

#include "tarch/tinyxml2/tinyxml2.h"

#include <array>
#include <iostream>
#include <string>

#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
#include <mpi.h>
#endif

#define MY_LINKEDCELL ls1::LS1RegionWrapper

class LS1MDSolverInterfaceTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(LS1MDSolverInterfaceTest);
  CPPUNIT_TEST(testConstants);
  CPPUNIT_TEST(testAddAndDeleteParticle);
  CPPUNIT_TEST(testGetCell);
  CPPUNIT_TEST(testGetCellIterator);
  CPPUNIT_TEST(testMassSync);
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp() {
#ifdef MARDYN_AUTOPAS
    _ls1ConfigFileName = "../test/unit/coupling/interface/impl/ls1/autopasgridconfig.xml";
#else
    _ls1ConfigFileName = "../test/unit/coupling/interface/impl/ls1/ls1gridconfig.xml";
#endif
    global_log = new Log::Logger(Log::None);
    global_log->set_mpi_output_root(0);
    coupling::interface::LS1StaticCommData::getInstance().setBoxOffsetAtDim(0, 0);
    coupling::interface::LS1StaticCommData::getInstance().setBoxOffsetAtDim(1, 0);
    coupling::interface::LS1StaticCommData::getInstance().setBoxOffsetAtDim(2, 0);
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
    coupling::interface::LS1StaticCommData::getInstance().setLocalCommunicator(MPI_COMM_WORLD);
    _domainGridDecomposition = {2, 2, 1}; // tests run at -np=4, hence a decomp of 2x2x1, cell calculation is hacked based on this
    coupling::interface::LS1StaticCommData::getInstance().setDomainGridDecompAtDim(0, _domainGridDecomposition[0]);
    coupling::interface::LS1StaticCommData::getInstance().setDomainGridDecompAtDim(1, _domainGridDecomposition[1]);
    coupling::interface::LS1StaticCommData::getInstance().setDomainGridDecompAtDim(2, _domainGridDecomposition[2]);
#endif
    _testSimulation = new Simulation();
    global_simulation = _testSimulation;
    _testSimulation->disableFinalCheckpoint();
    _testSimulation->readConfigFile(_ls1ConfigFileName);
    //_testSimulation->getDomain()->thermostatOff();
    //_testSimulation->getDomain()->setExplosionHeuristics(false);
    // after this point the mamico plugin exists and is accessible
    _testSimulation->prepare_start();
    _testSimulation->preSimLoopSteps();
  }
  void tearDown() {
    if (_testSimulation != nullptr) {
      delete _testSimulation;
      _testSimulation = nullptr;
    }
  }
  void testConstants() {
    // create interface
    coupling::interface::LS1MDSolverInterface interface({5, 5, 5}, {1, 1, 1});
    // load config file to read consts
    tinyxml2::XMLDocument conffile;
    tinyxml2::XMLElement* siteInfo = nullptr;
    conffile.LoadFile(_ls1ConfigFileName.c_str());
    tinyxml2::XMLHandle fileHandle(conffile);

    siteInfo = fileHandle.FirstChildElement("mardyn")
                   .FirstChildElement("simulation")
                   .FirstChildElement("ensemble")
                   .FirstChildElement("components")
                   .FirstChildElement("moleculetype")
                   .FirstChildElement("site")
                   .ToElement();
    if (siteInfo) {
      double mass, sigma, epsilon;
      siteInfo->FirstChildElement("mass")->QueryDoubleText(&mass);
      siteInfo->FirstChildElement("sigma")->QueryDoubleText(&sigma);
      siteInfo->FirstChildElement("epsilon")->QueryDoubleText(&epsilon);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(mass, interface.getMoleculeMass(), 1e-6);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(sigma, interface.getMoleculeSigma(), 1e-6);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(epsilon, interface.getMoleculeEpsilon(), 1e-6);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(interface.getKB(), 1.0, 1e-6);
    } else {
      CPPUNIT_FAIL("XML node could not be found");
    }

    siteInfo = fileHandle.FirstChildElement("mardyn").FirstChildElement("simulation").FirstChildElement("integrator").FirstChildElement("timestep").ToElement();
    if (siteInfo) {
      double timestep;
      siteInfo->QueryDoubleText(&timestep);
      CPPUNIT_ASSERT(timestep == interface.getDt());
    } else {
      CPPUNIT_FAIL("XML node could not be found");
    }
  }
  void testAddAndDeleteParticle() {
    // choose a position
    tarch::la::Vector<3, double> position = {2, 2, 2};

    // init
    double bBoxMin[3];
    double bBoxMax[3];
    global_simulation->domainDecomposition().getBoundingBoxMinMax(global_simulation->getDomain(), bBoxMin, bBoxMax);
    ls1::LS1RegionWrapper wrapper(bBoxMin, bBoxMax, _testSimulation);
    // create interface
    coupling::interface::LS1MDSolverInterface interface({5, 5, 5}, {1, 1, 1});

    // verify that position is empty
    bool found = false;
    while (wrapper.iteratorValid()) {
      ::Molecule* temp = wrapper.getParticleAtIterator();
      if (temp->r(0) == position[0] && temp->r(1) == position[1] && temp->r(2) == position[2]) {
        found = true;
        break;
      }
      wrapper.iteratorNext();
    }
    CPPUNIT_ASSERT_MESSAGE("Found before insertion", !found);

    // insert particle at position
    ::Molecule blank;
    coupling::interface::LS1Molecule tempParticle(&blank);
    tempParticle.setPosition(position);
    interface.addMoleculeToMDSimulation(tempParticle);

    // verify that position is filled
    found = false; // sanity
    wrapper.iteratorReset();
    while (wrapper.iteratorValid()) {
      ::Molecule* temp = wrapper.getParticleAtIterator();
      if (temp->r(0) == position[0] && temp->r(1) == position[1] && temp->r(2) == position[2]) {
        found = true;
        break;
      }
      wrapper.iteratorNext();
    }
    CPPUNIT_ASSERT_MESSAGE("Found after insertion", found == wrapper.isInRegion(position));

    // deletion
    interface.deleteMoleculeFromMDSimulation(tempParticle, wrapper);

    // verify
    found = false; // sanity
    wrapper.iteratorReset();
    while (wrapper.iteratorValid()) {
      ::Molecule* temp = wrapper.getParticleAtIterator();
      if (temp->r(0) == position[0] && temp->r(1) == position[1] && temp->r(2) == position[2]) {
        found = true;
        break;
      }
      wrapper.iteratorNext();
    }
    CPPUNIT_ASSERT_MESSAGE("Found after deletion", !found);
  }

  void testGetCell() {
    // setup
    // create interface with coupling cell size 10 and 2 linked cells per mac.cell per dimension, hence linked cells are size 5,5,5
    tarch::la::Vector<3, double> couplingCellSize(10.0);
    tarch::la::Vector<3, unsigned int> linkedCellsPerCouplingCell(2);
    // index converter
    tarch::la::Vector<3, unsigned int> globalNumberCells(4); // hence size becomes 40x40x40
    tarch::la::Vector<3, unsigned int> numberProcesses(1);
    int rank = 0;
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    for (int i = 0; i < 3; i++)
      numberProcesses[i] = _domainGridDecomposition[i];
#endif
    coupling::indexing::IndexingService<3>::getInstance().init(globalNumberCells, numberProcesses, coupling::paralleltopology::ZYX, 3, (unsigned int)rank);
    coupling::interface::LS1MDSolverInterface interface(couplingCellSize, linkedCellsPerCouplingCell);
    const coupling::IndexConversion<3> indexConversion(globalNumberCells, numberProcesses, rank, interface.getGlobalMDDomainSize(),
                                                       interface.getGlobalMDDomainOffset(), coupling::paralleltopology::ZYX);

    // first: try to get ghost cell, make sure there is error
    CPPUNIT_ASSERT_THROW(interface.getLinkedCell({0, 0, 0}, {1, 1, 1}, {2, 2, 2}, indexConversion), std::runtime_error);

    // then: get an internal cell and check whether the bounds are as expected
    for (unsigned int macX = 1; macX <= globalNumberCells[0] / numberProcesses[0]; macX++) {
      for (unsigned int macY = 1; macY <= globalNumberCells[1] / numberProcesses[1]; macY++) {
        for (unsigned int macZ = 1; macZ <= globalNumberCells[2] / numberProcesses[2]; macZ++) {
          for (unsigned int linkX = 0; linkX < linkedCellsPerCouplingCell[0]; linkX++) {
            for (unsigned int linkY = 0; linkY < linkedCellsPerCouplingCell[1]; linkY++) {
              for (unsigned int linkZ = 0; linkZ < linkedCellsPerCouplingCell[2]; linkZ++) {
                ls1::LS1RegionWrapper cell = interface.getLinkedCell({macX, macY, macZ}, {linkX, linkY, linkZ}, linkedCellsPerCouplingCell, indexConversion);

                tarch::la::Vector<3, double> receivedStart = {cell.getStartRegionAtDim(0), cell.getStartRegionAtDim(1), cell.getStartRegionAtDim(2)};
                tarch::la::Vector<3, double> receivedEnd = {cell.getEndRegionAtDim(0), cell.getEndRegionAtDim(1), cell.getEndRegionAtDim(2)};
                tarch::la::Vector<3, double> actualStart = {
                    // hacked for 2x2x1
                    ((macX - 1) * couplingCellSize[0]) + (linkX * couplingCellSize[0] / linkedCellsPerCouplingCell[0] +
                                                          ((rank / numberProcesses[0]) * couplingCellSize[0] * globalNumberCells[0] / numberProcesses[0])),
                    ((macY - 1) * couplingCellSize[1]) + (linkY * couplingCellSize[1] / linkedCellsPerCouplingCell[1] +
                                                          ((rank % numberProcesses[0]) * couplingCellSize[1] * globalNumberCells[1] / numberProcesses[1])),
                    ((macZ - 1) * couplingCellSize[2]) + (linkZ * couplingCellSize[2] / linkedCellsPerCouplingCell[2])};
                tarch::la::Vector<3, double> actualEnd = {actualStart[0] + couplingCellSize[0] / linkedCellsPerCouplingCell[0],
                                                          actualStart[1] + couplingCellSize[1] / linkedCellsPerCouplingCell[1],
                                                          actualStart[2] + couplingCellSize[2] / linkedCellsPerCouplingCell[2]};
                CPPUNIT_ASSERT_DOUBLES_EQUAL(receivedStart[0], actualStart[0], 1e-6);
                CPPUNIT_ASSERT_DOUBLES_EQUAL(receivedStart[1], actualStart[1], 1e-6);
                CPPUNIT_ASSERT_DOUBLES_EQUAL(receivedStart[2], actualStart[2], 1e-6);
                CPPUNIT_ASSERT_DOUBLES_EQUAL(receivedEnd[0], actualEnd[0], 1e-6);
                CPPUNIT_ASSERT_DOUBLES_EQUAL(receivedEnd[1], actualEnd[1], 1e-6);
                CPPUNIT_ASSERT_DOUBLES_EQUAL(receivedEnd[2], actualEnd[2], 1e-6);

              } // linkZ
            }   // linkY
          }     // linkX
        }       // macZ
      }         // macY
    }           // macX
  }

  void testGetCellIterator() {
    // identical setup as testGetCell()
    // create interface with coupling cell size 10 and 2 linked cells per mac.cell per dimension, hence linked cells are size 5,5,5
    tarch::la::Vector<3, double> couplingCellSize(10.0);
    tarch::la::Vector<3, unsigned int> linkedCellsPerCouplingCell(2);
    // index converter
    tarch::la::Vector<3, unsigned int> globalNumberCells(4); // hence size becomes 40x40x40
    tarch::la::Vector<3, unsigned int> numberProcesses(1);
    int rank = 0;
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    for (int i = 0; i < 3; i++)
      numberProcesses[i] = _domainGridDecomposition[i];
#endif
    coupling::indexing::IndexingService<3>::getInstance().init(globalNumberCells, numberProcesses, coupling::paralleltopology::ZYX, 3, (unsigned int)rank);
    coupling::interface::LS1MDSolverInterface interface(couplingCellSize, linkedCellsPerCouplingCell);
    const coupling::IndexConversion<3> indexConversion(globalNumberCells, numberProcesses, rank, interface.getGlobalMDDomainSize(),
                                                       interface.getGlobalMDDomainOffset(), coupling::paralleltopology::ZYX);

    // first: try to get ghost cell, make sure there is error
    CPPUNIT_ASSERT_THROW(interface.getLinkedCell({0, 0, 0}, {1, 1, 1}, {2, 2, 2}, indexConversion), std::runtime_error);

    // then: get internal cells and verify cell contents
    for (unsigned int macX = 1; macX <= globalNumberCells[0] / numberProcesses[0]; macX++) {
      for (unsigned int macY = 1; macY <= globalNumberCells[1] / numberProcesses[1]; macY++) {
        for (unsigned int macZ = 1; macZ <= globalNumberCells[2] / numberProcesses[2]; macZ++) {
          for (unsigned int linkX = 0; linkX < linkedCellsPerCouplingCell[0]; linkX++) {
            for (unsigned int linkY = 0; linkY < linkedCellsPerCouplingCell[1]; linkY++) {
              for (unsigned int linkZ = 0; linkZ < linkedCellsPerCouplingCell[2]; linkZ++) {
                ls1::LS1RegionWrapper cell = interface.getLinkedCell({macX, macY, macZ}, {linkX, linkY, linkZ}, linkedCellsPerCouplingCell, indexConversion);

                double actualStart[3] = {
                    // hacked for 2x2x1
                    ((macX - 1) * couplingCellSize[0]) + (linkX * couplingCellSize[0] / linkedCellsPerCouplingCell[0] +
                                                          ((rank / numberProcesses[0]) * couplingCellSize[0] * globalNumberCells[0] / numberProcesses[0])),
                    ((macY - 1) * couplingCellSize[1]) + (linkY * couplingCellSize[1] / linkedCellsPerCouplingCell[1] +
                                                          ((rank % numberProcesses[0]) * couplingCellSize[1] * globalNumberCells[1] / numberProcesses[1])),
                    ((macZ - 1) * couplingCellSize[2]) + (linkZ * couplingCellSize[2] / linkedCellsPerCouplingCell[2])};
                double actualEnd[3] = {actualStart[0] + couplingCellSize[0] / linkedCellsPerCouplingCell[0],
                                       actualStart[1] + couplingCellSize[1] / linkedCellsPerCouplingCell[1],
                                       actualStart[2] + couplingCellSize[2] / linkedCellsPerCouplingCell[2]};
                coupling::interface::MoleculeIterator<ls1::LS1RegionWrapper, 3>* receivedIterator = interface.getMoleculeIterator(cell);
                receivedIterator->begin();
                RegionParticleIterator ls1Iterator =
                    global_simulation->getMoleculeContainer()->regionIterator(actualStart, actualEnd, ParticleIterator::ALL_CELLS);

                while (ls1Iterator.isValid() && receivedIterator->continueIteration()) {
                  const coupling::interface::Molecule<3>& temp1(receivedIterator->get());
                  ::Molecule temp2(*ls1Iterator);

                  CPPUNIT_ASSERT_DOUBLES_EQUAL(temp2.r(0), temp1.getPosition()[0], 1e-6);
                  CPPUNIT_ASSERT_DOUBLES_EQUAL(temp2.r(1), temp1.getPosition()[1], 1e-6);
                  CPPUNIT_ASSERT_DOUBLES_EQUAL(temp2.r(2), temp1.getPosition()[2], 1e-6);
                  CPPUNIT_ASSERT_DOUBLES_EQUAL(temp2.v(0), temp1.getVelocity()[0], 1e-6);
                  CPPUNIT_ASSERT_DOUBLES_EQUAL(temp2.v(1), temp1.getVelocity()[1], 1e-6);
                  CPPUNIT_ASSERT_DOUBLES_EQUAL(temp2.v(2), temp1.getVelocity()[2], 1e-6);
                  ++ls1Iterator;
                  receivedIterator->next();
                  CPPUNIT_ASSERT(ls1Iterator.isValid() == receivedIterator->continueIteration());
                }
              } // linkZ
            }   // linkY
          }     // linkX
        }       // macZ
      }         // macY
    }           // macX
  }

  void testMassSync() {
    // setup
    // create interface with coupling cell size 10 and 2 linked cells per mac.cell per dimension, hence linked cells are size 5,5,5
    tarch::la::Vector<3, double> couplingCellSize(10.0);
    tarch::la::Vector<3, unsigned int> linkedCellsPerCouplingCell(2);
    coupling::interface::LS1MDSolverInterface interface(couplingCellSize, linkedCellsPerCouplingCell);
    int rank = 0;
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
    // create a particle, ensure that position is empty
    tarch::la::Vector<3, double> position = {0.2, 0.2, 0.2};

    // init
    double bBoxMin[3];
    double bBoxMax[3];
    global_simulation->domainDecomposition().getBoundingBoxMinMax(global_simulation->getDomain(), bBoxMin, bBoxMax);
    for (int i = 0; i < 3; i++) {
      bBoxMin[i] = bBoxMin[i] - _testSimulation->getcutoffRadius();
      bBoxMax[i] = bBoxMax[i] + _testSimulation->getcutoffRadius();
    }

    ls1::LS1RegionWrapper wrapper(bBoxMin, bBoxMax, _testSimulation);

    // initial mass sync to populate halos
    interface.synchronizeMoleculesAfterMassModification();

    // verify that position is empty
    bool found = false;
    while (wrapper.iteratorValid()) {
      ::Molecule* temp = wrapper.getParticleAtIterator();
      if (temp->r(0) == position[0] && temp->r(1) == position[1] && temp->r(2) == position[2]) {
        found = true;
        break;
      }
      wrapper.iteratorNext();
    }
    CPPUNIT_ASSERT(!found);

    // insert the particle
    ::Molecule blank;
    coupling::interface::LS1Molecule tempParticle(&blank);
    tempParticle.setPosition(position);
    interface.addMoleculeToMDSimulation(tempParticle);

    // create all shifted positions
    std::array<tarch::la::Vector<3, double>, 7> shiftedPositions;
    tarch::la::Vector<3, double> domainSpan = interface.getGlobalMDDomainSize();
    for (int i = 1; i <= 7; i++) {
      int copy = i;
      shiftedPositions[i - 1][0] = position[0] + (copy & 1 ? domainSpan[0] : 0);
      copy >>= 1;
      shiftedPositions[i - 1][1] = position[1] + (copy & 1 ? domainSpan[1] : 0);
      copy >>= 1;
      shiftedPositions[i - 1][2] = position[2] + (copy & 1 ? domainSpan[2] : 0);
    }

    // verify that the halo particle doesn't yet exist
    for (tarch::la::Vector<3, double> curCheck : shiftedPositions) {
      wrapper.iteratorReset();
      bool found = false;
      while (wrapper.iteratorValid()) {
        ::Molecule* temp = wrapper.getParticleAtIterator();
        if (temp->r(0) == curCheck[0] && temp->r(1) == curCheck[1] && temp->r(2) == curCheck[2]) {
          found = true;
          break;
        }
        wrapper.iteratorNext();
      }
      CPPUNIT_ASSERT(!found);
    }

    // sync
    interface.synchronizeMoleculesAfterMassModification();

    // verify that the particle now exists
    for (tarch::la::Vector<3, double> curCheck : shiftedPositions) {
      wrapper.iteratorReset();
      bool found = false;
      while (wrapper.iteratorValid()) {
        ::Molecule* temp = wrapper.getParticleAtIterator();
        if (temp->r(0) == curCheck[0] && temp->r(1) == curCheck[1] && temp->r(2) == curCheck[2]) {
          found = true;
          break;
        }
        wrapper.iteratorNext();
      }
      CPPUNIT_ASSERT(found == wrapper.isInRegion(curCheck));
    }

    // delete the particle
    interface.deleteMoleculeFromMDSimulation(tempParticle, wrapper);
    interface.synchronizeMoleculesAfterMassModification(); // even though the interface no longer exists, sim does
  }

private:
  Simulation* _testSimulation;
  std::string _ls1ConfigFileName;
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
  std::array<int, 3> _domainGridDecomposition;
#endif
};

CPPUNIT_TEST_SUITE_REGISTRATION(LS1MDSolverInterfaceTest);