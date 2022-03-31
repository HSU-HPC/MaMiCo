// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_TESTS_MACROSCOPICCELLSTEST_H_
#define _MOLECULARDYNAMICS_COUPLING_TESTS_MACROSCOPICCELLSTEST_H_

#include "coupling/CouplingMDDefinitions.h"
#include "coupling/IndexConversion.h"
#include "coupling/datastructures/MacroscopicCells.h"
#include "coupling/interface/MDSolverInterface.h"
#include "coupling/tests/Test.h"
#include "simplemd/LinkedCell.h"
#include "tarch/la/Vector.h"

/** tests the access to the macroscopic cells in class MacroscopicCells.
 *  Therefore, we define two classes which just set the local vector cell
 * indices in the macroscopic cells
 *  and print the values.
 *  We loop over both ghost cells and non-ghost cells and can check the output.
 *  This test is carried out for 1D,2D,3D.
 *  @author Philipp Neumann
 */
class MacroscopicCellsTest : public Test {
public:
  MacroscopicCellsTest() : Test("MacroscopicCellsTest") {}
  virtual ~MacroscopicCellsTest() {}

  virtual void run() {
// this test is only required in sequential mode
#if (COUPLING_MD_PARALLEL == COUPLING_MD_NO)
    std::cout << "MacroscopicCellsTest: run 1D test" << std::endl;
    testMacroscopicCells<1>();
    std::cout << "MacroscopicCellsTest: run 2D test" << std::endl;
    testMacroscopicCells<2>();
    std::cout << "MacroscopicCellsTest: run 3D test" << std::endl;
    testMacroscopicCells<3>();
#endif
  }

private:
  /** test interface; only implements the getLinkedCell method (required in
   * initialisation of macroscopic cells) and
   *  strictly returns its only linked cell ;-) The MD domain is considered to
   * be the unit square/box.
   */
  template <unsigned int dim> class TestMDSolverInterface : public coupling::interface::MDSolverInterface<simplemd::LinkedCell, dim> {
  public:
    TestMDSolverInterface() : coupling::interface::MDSolverInterface<simplemd::LinkedCell, dim>() {
    } virtual ~TestMDSolverInterface() {
    }

    virtual simplemd::LinkedCell &getLinkedCell(const tarch::la::Vector<dim, unsigned int> &macroscopicCellIndex,
                                                const tarch::la::Vector<dim, unsigned int> &linkedCellInMacroscopicCell,
                                                const tarch::la::Vector<dim, unsigned int> &linkedCellsPerMacroscopicCell,
                                                const coupling::IndexConversion<dim> &indexConversion) {
      return _linkedcell;
    }

    /** returns the global size of the box-shaped MD domain */
    virtual tarch::la::Vector<dim, double> getGlobalMDDomainSize() const { return tarch::la::Vector<dim, double>(1.0); }

    /** returns the offset (i.e. lower,left corner) of MD domain */
    virtual tarch::la::Vector<dim, double> getGlobalMDDomainOffset() const { return tarch::la::Vector<dim, double>(0.0); }

    virtual double getMoleculeMass() const { return 1.0; }
    virtual double getKB() const { return 1.0; }
    virtual double getMoleculeSigma() const { return 1.0; }
    virtual double getMoleculeEpsilon() const { return 1.0; }
    virtual void getInitialVelocity(const tarch::la::Vector<dim, double> &meanVelocity, const double &kB, const double &temperature,
                                    tarch::la::Vector<dim, double> &initialVelocity) const {}
    virtual void deleteMoleculeFromMDSimulation(const coupling::interface::Molecule<dim> &molecule, simplemd::LinkedCell &cell) {}
    virtual void addMoleculeToMDSimulation(const coupling::interface::Molecule<dim> &molecule) {}
    virtual void setupPotentialEnergyLandscape(const tarch::la::Vector<dim, unsigned int> &indexOfFirstMacroscopicCell,
                                               const tarch::la::Vector<dim, unsigned int> &rangeMacroscopicCells,
                                               const tarch::la::Vector<dim, unsigned int> &linkedCellsPerMacroscopicCell) {}
    virtual tarch::la::Vector<dim, unsigned int> getLinkedCellIndexForMoleculePosition(const tarch::la::Vector<dim, double> &position) {
      return tarch::la::Vector<dim, unsigned int>(0);
    }
    virtual void calculateForceAndEnergy(coupling::interface::Molecule<dim> &molecule) {}
    virtual void synchronizeMoleculesAfterMassModification() {}
    virtual void synchronizeMoleculesAfterMomentumModification() {}
    virtual double getDt() { return 1.0; }
    virtual coupling::interface::MoleculeIterator<simplemd::LinkedCell, dim> *getMoleculeIterator(simplemd::LinkedCell &cell) {
      return NULL;
    }

    private : simplemd::LinkedCell _linkedcell;
  };

  /** test class for loop over macroscopic cells; writes the local vector cell
   * index into each cell's microscopic-momentum
   *  buffer.
   */
  template <unsigned int dim> class WriteIndex {
  public:
    WriteIndex(const coupling::IndexConversion<dim> &indexConversion) : _indexConversion(indexConversion) {}
    ~WriteIndex() {}
    void beginCellIteration() {}
    void endCellIteration() {}
    void apply(coupling::datastructures::MacroscopicCellWithLinkedCells<simplemd::LinkedCell, dim> &cell, const unsigned int &index) {
      tarch::la::Vector<dim, unsigned int> localIndex = _indexConversion.getLocalVectorCellIndex(index);
      tarch::la::Vector<dim, double> convertLocal(0.0);
      for (unsigned int d = 0; d < dim; d++) {
        convertLocal[d] = localIndex[d];
      }
      cell.setMicroscopicMomentum(convertLocal);
    }

  private:
    const coupling::IndexConversion<dim> &_indexConversion;
  };
  /** test class to print local cell index (-> microscopic momentum buffer) of
   * macroscopic cells */
  template <unsigned int dim> class PrintIndex {
  public:
    PrintIndex() {}
    ~PrintIndex() {}
    void beginCellIteration() {}
    void endCellIteration() {}
    void apply(coupling::datastructures::MacroscopicCellWithLinkedCells<simplemd::LinkedCell, dim> &cell, const unsigned int &index) {
      std::cout << cell.getMicroscopicMomentum() << std::endl;
    }
  };

  template <unsigned int dim> void testMacroscopicCells() {
    TestMDSolverInterface<dim> *testInterface = new TestMDSolverInterface<dim>();
    if (testInterface == NULL) {
      std::cout << "ERROR MacroscopicCellsTest::testMacroscopicCells(): "
                   "testInterface==NULL!"
                << std::endl;
      exit(EXIT_FAILURE);
    }
    const tarch::la::Vector<dim, unsigned int> globalNumberCells(3);
    const tarch::la::Vector<dim, unsigned int> numberProcesses(1);
    const unsigned int rank = 0;
    const coupling::IndexConversion<dim> indexConversion(globalNumberCells, numberProcesses, rank, testInterface->getGlobalMDDomainSize(),
                                                         testInterface->getGlobalMDDomainOffset(), coupling::paralleltopology::XYZ);
    const tarch::la::Vector<dim, unsigned int> numberLinkedCellsPerMacroscopicCell(2);

    coupling::datastructures::MacroscopicCells<simplemd::LinkedCell, dim> macroscopicCells(numberLinkedCellsPerMacroscopicCell, indexConversion, testInterface);
    WriteIndex<dim> writeIndex(indexConversion);
    PrintIndex<dim> printIndex;

    // write indices to inner macroscopic cells
    macroscopicCells.applyToLocalNonGhostMacroscopicCellsWithLinkedCells(writeIndex);
    macroscopicCells.applyToLocalNonGhostMacroscopicCellsWithLinkedCells(printIndex);

    // write indices to outer macroscopic cells
    macroscopicCells.applyToLocalGhostMacroscopicCellsWithLinkedCells(writeIndex);
    macroscopicCells.applyToLocalGhostMacroscopicCellsWithLinkedCells(printIndex);

    delete testInterface;
    testInterface = NULL;
  }
};

#endif // _MOLECULARDYNAMICS_COUPLING_TESTS_MACROSCOPICCELLSTEST_H_
