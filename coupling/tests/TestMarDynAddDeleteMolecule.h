// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef TESTMARDYNADDDELETEMOLECULE_H_
#define TESTMARDYNADDDELETEMOLECULE_H_

#include "TestMarDyn.h"

#include "particleContainer/LinkedCells.h"
#include "particleContainer/ParticleCell.h"

#include "coupling/interface/Molecule.h"
#include "coupling/interface/impl/MarDyn/MarDynCell.h"
#include "coupling/interface/impl/MarDyn/MarDynCoupledSimulation.h"

/**
 * 	test class for the MarDyn MD coupling interface
 * 	Tests the deletion and then insertion of a molecule into the
 *ParticleContainer
 *	@author Hanno Flohr
 */
class TestMarDynAddDeleteMolecule : public TestMarDyn {
public:
  TestMarDynAddDeleteMolecule(int argc, char **argv, std::string name)
      : TestMarDyn(argc, argv, name) {}

  virtual ~TestMarDynAddDeleteMolecule() {}

  virtual void run() {
    std::cout << "MarDyn coupling - Add/Delete Tests: " << std::endl;
    std::cout << "Initializing MarDyn MD solver interface.." << std::endl;
    this->loadMarDynTestConfiguration("mardyn_dummy_coupling.cfg", 1);
    this->loadMacroscopicSolverConfiguration();
    this->loadMamicoTestConfiguration();

    std::cout << "Test 1: delete molecule:" << std::endl;
    deleteMoleculeTest();

    std::cout << "Test 2: add molecule:" << std::endl;
    addMoleculeTest();
  }

private:
  void deleteMoleculeTest() {
    MarDynMDSolverInterface *mdsi =
        (MarDynMDSolverInterface *)coupling::interface::MamicoInterfaceProvider<
            MarDynCell, 3>::getInstance()
            .getMDSolverInterface();

    MarDynCoupledSimulation *sim = mdsi->getSimulation();

    // molecule count before deletion
    unsigned long int moleculeCountBefore =
        sim->getMolecules()->getNumberOfParticles();
    std::cout << "Test moleculeCount " << moleculeCountBefore << std::endl;

    LinkedCells *lc = (LinkedCells *)mdsi->getSimulation()->getMolecules();

    Molecule *molecule = lc->begin();

    // molecule values
    tarch::la::Vector<3, double> position(0.0);
    tarch::la::Vector<3, double> velocity(0.0);
    double force[3] = {0.0, 0.0, 0.0};
    for (int d = 0; d < 3; d++) {
      position[d] = molecule->r(d);
      velocity[d] = molecule->v(d);
      force[d] = molecule->F(d);
    }

    // create interface molecule
    Molecule internalMolecule;
    for (int d = 0; d < 3; d++) {
      internalMolecule.setr(d, position[d]);
      internalMolecule.setv(d, velocity[d]);
    }
    internalMolecule.setF(force);
    MarDynMoleculeWrapper testMolecule;
    testMolecule.setMolecule(&internalMolecule, sim->getcutoffRadius());

    coupling::services::MacroscopicCellService<3> *macroCellService =
        coupling::interface::MamicoInterfaceProvider<MarDynCell,
                                                     3>::getInstance()
            .getMacroscopicCellService();

    const tarch::la::Vector<3, unsigned int> linkedCellsPerMacroscopicCell(1);
    const tarch::la::Vector<3, unsigned int> linkedCellInMacroscopicCell(0);
    tarch::la::Vector<3, unsigned int> deleteCellIndexVector(1);

    // macroscopic cell index
    tarch::la::Vector<3, unsigned int> macroCellIndex(1);
    for (int d = 0; d < 3; d++)
      macroCellIndex[d] = 1 + (int)floor(position[d] / 2.5);

    std::cout << "macro cell index: " << macroCellIndex << std::endl;

    MarDynCell &mc = mdsi->getLinkedCell(
        macroCellIndex, linkedCellInMacroscopicCell,
        linkedCellsPerMacroscopicCell, macroCellService->getIndexConversion());

    mdsi->deleteMoleculeFromMDSimulation(testMolecule, mc);
    std::cout << "molecule count in cell after: "
              << mc.getParticleCell()->getParticlePointers().size()
              << std::endl;

    unsigned long int moleculeCountAfter =
        sim->getMolecules()->getNumberOfParticles();
    std::cout << "moleculeCount: " << moleculeCountBefore << " >? "
              << moleculeCountAfter << std::endl;
    if (moleculeCountBefore != (moleculeCountAfter + 1))
      std::cout << "Molecule deletion test failed!" << std::endl;
    else
      std::cout << "Molecule deletion test successful!" << std::endl;

    std::cout << "synchronize testing..." << std::endl;
    mdsi->synchronizeMoleculesAfterMassModification();
    unsigned long int moleculeCountAfterSync =
        sim->getMolecules()->getNumberOfParticles();
    if (moleculeCountAfter == moleculeCountAfterSync)
      std::cout << "synchronize test successfull!" << std::endl;
    else
      std::cout << "synchronize test failed: " << moleculeCountAfterSync
                << "in particleContainer after synchronize." << std::endl;
  }

  void addMoleculeTest() {
    // the md solver interface
    MarDynMDSolverInterface *mdsi =
        (MarDynMDSolverInterface *)coupling::interface::MamicoInterfaceProvider<
            MarDynCell, 3>::getInstance()
            .getMDSolverInterface();

    // molecule values
    tarch::la::Vector<3, double> position(0.0);
    double dummy[3] = {0.0, 0.0, 0.0};
    for (int d = 0; d < 3; d++)
      position[d] = 15.0;
    Molecule internalMolecule;
    for (int d = 0; d < 3; d++) {
      internalMolecule.setr(d, position[d]);
      internalMolecule.setv(d, dummy[d]);
    }
    internalMolecule.setF(dummy);

    // the molecule that will be inserted
    MarDynMoleculeWrapper testMolecule;
    testMolecule.setMolecule(&internalMolecule,
                             mdsi->getSimulation()->getcutoffRadius());

    unsigned long int moleculesBefore =
        mdsi->getSimulation()->getMolecules()->getNumberOfParticles();

    // add molecule to md sim via simulation interface
    mdsi->addMoleculeToMDSimulation(testMolecule);

    unsigned long int moleculesAfter =
        mdsi->getSimulation()->getMolecules()->getNumberOfParticles();

    std::cout << "moleculeCount before: " << moleculesBefore
              << " | after: " << moleculesAfter << std::endl;
    if (moleculesBefore != (moleculesAfter - 1))
      std::cout << "Molecule insertion test failed!" << std::endl;
    else
      std::cout << "Molecule insertion test successful!" << std::endl;
  }
};

#endif /* TESTMARDYNADDDELETEMOLECULE_H_ */
