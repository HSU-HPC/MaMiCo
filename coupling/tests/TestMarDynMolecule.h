// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef TESTMARDYNMOLECULE_H_
#define TESTMARDYNMOLECULE_H_

#include "TestMarDyn.h"

#include "coupling/interface/impl/MarDyn/MarDynMoleculeWrapper.h"
#include "coupling/interface/impl/MarDyn/MarDynMoleculeIterator.h"

#include "particleContainer/LinkedCells.h"

/**
 * Tests for the interfaces:
 * First a test on a single molecule is executed to test the get/set-methods of
 * the MarDynMoleculeWrapper.
 * Then a test using a whole cell (and the MarDynMoleculeIterator) is executed,
 * in which the values of the
 * molecules (position, velocity, force) obtained from the MarDyn simulation
 * directly are compared with the
 * values obtained through the interfaces.
 * The third test checks the results of the getPotentialEnergy() and
 * calculateForceAndEnergy()-methods
 * @author Hanno Flohr
 */
class TestMarDynMolecule : public TestMarDyn {
public:
  TestMarDynMolecule(int argc, char **argv, std::string name)
      : TestMarDyn(argc, argv, name) {}
  virtual ~TestMarDynMolecule() {}

  virtual void run() {
    std::cout << "Testing MarDynMoleculeWrapper interface:" << std::endl;
    std::cout << "Initializing MarDyn MD solver interface.." << std::endl;
    //init simulation and set md solver interface
    this->loadMarDynTestConfiguration("mardyn_dummy_coupling.cfg", 1);

    int errors = 0;

    std::cout << "Test 1: Get/set-methods: " << std::endl;
    getSetTest(true, errors);
    if (errors == 0)
      std::cout << "TestMarDynMolecule - Test 1 was successful! " << std::endl;
    else {
      std::cout << "TestMarDynMolecule - Test 1 not successful, " << errors
                << " errors occurred. Test 2 will not be executed."
                << std::endl;
      exit(1);
    }

    errors = 0;
    std::cout
        << "Test 2: Comparison of values from simulation and through interface:"
        << std::endl;
    int currentCellIndex = 0;
    while (currentCellIndex < 2744) {
      compareMoleculeList(false, currentCellIndex, errors);
    }
    if (errors == 0)
      std::cout << "TestMarDynMolecule - Test 2 was successful! " << std::endl;
    else
      std::cout << "TestMarDynMolecule - Test 2 was not successful, " << errors
                << " errors occurred." << std::endl;

    std::cout << "Test 3: Potential Energy Test: " << std::endl;
    errors = 0;
    potEnergyAndForceTest(true, errors);
    if (errors == 0)
      std::cout << "TestMarDynMolecule - Test 3 was successful! " << std::endl;
    else
      std::cout << "TestMarDynMolecule - Test 3 was not successful, " << errors
                << " errors occurred." << std::endl;
  }

protected:
  void potEnergyAndForceTest(bool detail, int &errors) {
    MarDynMDSolverInterface *mdsi =
        (MarDynMDSolverInterface *)coupling::interface::MamicoInterfaceProvider<
            MarDynCell, 3>::getInstance().getMDSolverInterface();
    double potentialEnergy = 0.0;

    //get particle container from simulation
    LinkedCells *cells = (LinkedCells *)mdsi->getSimulation()->getMolecules();

    Molecule *curMolecule;
    int counter = 0;
    for (curMolecule = cells->begin(); curMolecule != cells->end();
         curMolecule = cells->next()) {
      counter++;

      //			std::cout << "cellIndex is: " <<
      //cells->getCellIndexOfMolecule(curMolecule) << std::endl;
      //			std::cout << "pos       : ";
      //			for(int d=0; d<3; d++) std::cout << curMolecule->r(d) << ", ";
      //			std::cout << std::endl;

      //molecule values
      tarch::la::Vector<3, double> position(0.0);
      tarch::la::Vector<3, double> velocity(0.0);
      double force[3] = { 0.0, 0.0, 0.0 };
      for (int d = 0; d < 3; d++) {
        position[d] = curMolecule->r(d);
        velocity[d] = curMolecule->v(d);
        force[d] = curMolecule->F(d);
      }

      //create new MarDynMoleculeWrapper
      Molecule internalMolecule;
      for (unsigned int d = 0; d < 3; d++) {
        internalMolecule.setr(d, position[d]);
        internalMolecule.setv(d, velocity[d]);
      }
      internalMolecule.setF(force);
      MarDynMoleculeWrapper molecule;
      molecule.setMolecule(&internalMolecule,
                           mdsi->getSimulation()->getcutoffRadius());

      //compute force and potential energy
      mdsi->calculateForceAndEnergy(molecule);
      potentialEnergy = molecule.getPotentialEnergy();

      if (detail)
        std::cout << "molecule id: " << curMolecule->id() << std::endl;
      //std::cout << "molecule pos       : " << position << std::endl;
      if (detail)
        std::cout << "- potentialEnergy: " << potentialEnergy << std::endl;
      if (potentialEnergy == 0.0)
        errors++;
      if (detail)
        std::cout << "- force: " << molecule.getForce() << std::endl;
      if (molecule.getForce()[0] == 0.0 && molecule.getForce()[1] == 0.0 &&
          molecule.getForce()[2] == 0.0)
        errors++;

      //only test 1000 molecules
      if (counter == 1000)
        break;
    }

    //-------------------------------------------------------------------------
    //test getPotentialEnergy method for molecule that is not in the simulation
    std::cout << "Testing getPotentialEnergy() for a new molecule.."
              << std::endl;

    tarch::la::Vector<3, double> newPosition(23.0);
    tarch::la::Vector<3, double> newVelocity(1.0);
    double newForce[3] = { 1.0, 1.0, 1.0 };
    Molecule internalMolecule;
    for (unsigned int d = 0; d < 3; d++) {
      internalMolecule.setr(d, newPosition[d]);
      internalMolecule.setv(d, newVelocity[d]);
    }
    internalMolecule.setF(newForce);
    MarDynMoleculeWrapper newMolecule;
    newMolecule.setMolecule(&internalMolecule,
                            mdsi->getSimulation()->getcutoffRadius());

    potentialEnergy = newMolecule.getPotentialEnergy();
    if (detail)
      std::cout << "potentialEnergy of the new molecule: " << potentialEnergy
                << std::endl;
    if (potentialEnergy != potentialEnergy) {
      std::cout << "Potential energy fault for new molecule!" << std::endl;
      errors++;
    }
  }

  void getSetTest(bool detail, int &errors) {
    //get the first molecule from the simulation
    Molecule *molecule = this->_marDyn->getMolecules()->begin();

    //input values
    tarch::la::Vector<3, double> inputPosition(0.0);
    tarch::la::Vector<3, double> inputVelocity(0.0);
    tarch::la::Vector<3, double> inputForce(0.0);
    for (int d = 0; d < 3; d++) {
      inputPosition[d] = molecule->r(d);
      inputVelocity[d] = molecule->v(d);
      inputForce[d] = molecule->F(d);
    }

    if (detail)
      std::cout << "Values of first molecule:     " << inputPosition << " | "
                << inputVelocity << " | " << inputForce << std::endl;

    //create MarDynMoleculeWrapper
    tarch::la::Vector<3, double> position(0.0);
    tarch::la::Vector<3, double> velocity(0.0);
    double force[3] = { 0.0, 0.0, 0.0 };
    for (int d = 0; d < 3; d++) {
      position[d] = molecule->r(d);
      velocity[d] = molecule->v(d);
      force[d] = molecule->F(d);
    }
    Molecule internalMolecule;
    for (int d = 0; d < 3; d++) {
      internalMolecule.setr(d, position[d]);
      internalMolecule.setv(d, velocity[d]);
    }
    internalMolecule.setF(force);
    MarDynMoleculeWrapper mdmw;
    mdmw.setMolecule(&internalMolecule, this->_marDyn->getcutoffRadius());

    //test values
    tarch::la::Vector<3, double> testPosition = mdmw.getPosition();
    tarch::la::Vector<3, double> testVelocity = mdmw.getVelocity();
    tarch::la::Vector<3, double> testForce = mdmw.getForce();
    if (detail)
      std::cout << "Values of interface molecule: " << inputPosition << " | "
                << inputVelocity << " | " << inputForce << std::endl;

    if (testPosition != inputPosition) {
      std::cout << "ERROR: getPosition: " << testPosition << " should be "
                << inputPosition << std::endl;
      errors++;
    }
    if (testVelocity != inputVelocity) {
      std::cout << "ERROR: getVelocity: " << testVelocity << " should be "
                << inputVelocity << std::endl;
      errors++;
    }
    if (testForce != inputForce) {
      std::cout << "ERROR: getForce: " << testForce << " should be "
                << inputForce << std::endl;
      errors++;
    }

    //value to test the set-methods
    tarch::la::Vector<3, double> setValue(1.0);
    mdmw.setPosition(setValue);
    mdmw.setVelocity(setValue);
    mdmw.setForce(setValue);

    testPosition = mdmw.getPosition();
    testVelocity = mdmw.getVelocity();
    testForce = mdmw.getForce();

    if (detail)
      std::cout << "Molecule values after set:    " << testPosition << " | "
                << testVelocity << " | " << testForce << std::endl;

    if (testPosition != setValue) {
      std::cout << "ERROR: setPosition: " << testPosition << " should be "
                << setValue << std::endl;
      errors++;
    }
    if (testVelocity != setValue) {
      std::cout << "ERROR: setVelocity: " << testVelocity << " should be "
                << setValue << std::endl;
      errors++;
    }
    if (testForce != setValue) {
      std::cout << "ERROR: setForce: " << testForce << " should be " << setValue
                << std::endl;
      errors++;
    }
  }

  void compareMoleculeList(bool detail, int &currentCellIndex, int &errors) {
    //get particle container from simulation
    LinkedCells *cells = (LinkedCells *)this->_marDyn->getMolecules();

    //for testing get first cell and create MarDynCell
    ParticleCell pc;
    bool cellFound = false;
    for (int i = currentCellIndex; i < 2744; i++) {
      pc = cells->getCell(i);
      if (pc.getMoleculeCount() > 0) {
        std::cout << "Test 2 " << i << " " << pc.getMoleculeCount()
                  << std::endl;
        currentCellIndex = i + 1;
        cellFound = true;
        break;
      }
    }
    if (!cellFound) {
      currentCellIndex = 9999;
      return;
    }
    MarDynCell cell = MarDynCell(&pc, cells->getCutoff());

    //vector containing all molecules of the particle cell
    std::vector<Molecule *> moleculesInParticleCell = pc.getParticlePointers();

    //set up molecule iterator for testing
    MarDynMoleculeIterator mdmi(cell);

    //helper variables
    int moleculeCount = pc.getMoleculeCount();
    std::vector<tarch::la::Vector<3, double> > listMarDyn;
    std::vector<tarch::la::Vector<3, double> > listMarDynCouplingInterface;
    tarch::la::Vector<3, double> position(0.0);
    tarch::la::Vector<3, double> velocity(0.0);
    tarch::la::Vector<3, double> force(0.0);
    int counter = 0;

    //loop through all molecules in the cell
    for (mdmi.begin(); mdmi.continueIteration(); mdmi.next()) {
      if (detail)
        std::cout << "Molecule #" << counter << ":" << std::endl;
      //get position and velocity through interface
      position = mdmi.get().getPosition();
      velocity = mdmi.get().getVelocity();
      force = mdmi.get().getForce();
      if (detail)
        std::cout << "interface: " << position << " | " << velocity << " | "
                  << force << std::endl;
      listMarDynCouplingInterface.push_back(position);
      listMarDynCouplingInterface.push_back(velocity);
      listMarDynCouplingInterface.push_back(force);

      //get position and velocity from the ParticleCell for comparison
      for (int d = 0; d < 3; d++) {
        position[d] = moleculesInParticleCell[counter]->r(d);
        velocity[d] = moleculesInParticleCell[counter]->v(d);
        force[d] = moleculesInParticleCell[counter]->F(d);
      }
      if (detail)
        std::cout << "cell_____: " << position << " | " << velocity << " | "
                  << force << std::endl;
      listMarDyn.push_back(position);
      listMarDyn.push_back(velocity);
      listMarDyn.push_back(force);

      counter++;
    }

    if (listMarDyn != listMarDynCouplingInterface) {
      errors = 1;
      std::cout
          << "ERROR: the molecule values obtained from the MarDyn simulation "
             "and through the Coupling interface do not match! (Cell #"
          << currentCellIndex << ")" << std::endl;
    }
  }

};

#endif /* TESTMARDYNMOLECULE_H_ */
