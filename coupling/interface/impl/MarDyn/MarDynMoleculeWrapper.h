// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef MARDYNMOLECULEWRAPPER_H_
#define MARDYNMOLECULEWRAPPER_H_

// MarDyn
#include "Simulation.h"
#include "molecules/Molecule.h"
typedef Molecule MardynMolecule;
#include "particleContainer/LinkedCells.h"
#include "particleContainer/ParticleCell.h"
#include "particleContainer/adapter/LegacyCellProcessor.h"
// MaMiCo
#include "coupling/CouplingMDDefinitions.h"
#include "coupling/interface/Molecule.h"
#include "coupling/interface/impl/MarDyn/MarDynCell.h"
typedef coupling::interface::Molecule<3> MamicoMolecule;
#include "tarch/la/Vector.h"
#include <cmath>
#include <stdio.h>

/* 	Implementation of the MoleculeWrapper interface for MarDyn.
 *      We currently buffer position, velocity and force in tarch-vectors to
 * speed up the simulation. This is only valid, as long as the wrapper is indeed
 * a wrapper, i.e. the MardynMolecule pointed to by _myMolecule is only modified
 * through this wrapper and not by any other mechanism!!!
 * 	@author Hanno Flohr, Philipp Neumann
 */
class MarDynMoleculeWrapper : public coupling::interface::Molecule<3> {

public:
  // Constructor needed for the MarDynMoleculeIterator
  MarDynMoleculeWrapper() : _myMolecule(NULL), _cutoffRadius(0.0) {}

  virtual ~MarDynMoleculeWrapper() {}

  // sets the MarDynMolecule and the cutoffRadius (needed in the
  // MoleculeIterator interface)
  void setMolecule(MardynMolecule *molecule, double cutoffRadius) {
    _myMolecule = molecule;
    _cutoffRadius = cutoffRadius;
  }

  /* get/set the velocity of the molecule */
  virtual tarch::la::Vector<3, double> getVelocity() const {
    tarch::la::Vector<3, double> myVelocity(_myMolecule->v(0), _myMolecule->v(1), _myMolecule->v(2));
    return myVelocity;
  }
  virtual void setVelocity(const tarch::la::Vector<3, double> &velocity) {
    for (unsigned int d = 0; d < 3; d++)
      _myMolecule->setv(d, velocity[d]);
  }

  /* get/set the position of the molecule */
  virtual tarch::la::Vector<3, double> getPosition() const {
    tarch::la::Vector<3, double> myPosition(_myMolecule->r(0), _myMolecule->r(1), _myMolecule->r(2));
    return myPosition;
  }
  virtual void setPosition(const tarch::la::Vector<3, double> &position) {
    for (unsigned int d = 0; d < 3; d++)
      _myMolecule->setr(d, position[d]);
  }

  /* get/set the force on the molecule */
  virtual tarch::la::Vector<3, double> getForce() const {
    tarch::la::Vector<3, double> myForce(_myMolecule->F(0), _myMolecule->F(1), _myMolecule->F(2));
    return myForce;
  }
  virtual void setForce(const tarch::la::Vector<3, double> &force) {
    double Force[3] = {force[0], force[1], force[2]};
    _myMolecule->setF(Force);
  }

  /* get/set potential energy of the molecule */
  virtual double getPotentialEnergy() const {
    double potentialEnergy = 0.0;

    // pointer to the MarDyn simulation
    MarDynCoupledSimulation *mdSim = (MarDynCoupledSimulation *)global_simulation;
    LegacyCellProcessor legacyCellProcessor(_cutoffRadius, mdSim->getLJCutoff(), mdSim->getTersoffCutoff(), mdSim->getParticlePairsHandler());

    // if myMolecule is not in the particle container yet (check id):
    //		the molecule was initialized by the coupling part
    //		thus the correct molecule in the simulation has to be found
    // before
    // computing the energy else 		the molecule was initialized using a pointer
    // to the molecule in the simulation, no search needed
    if (_myMolecule->id() < 1 || _myMolecule->id() > global_simulation->getMaxID()) {
      int cellIndex;
      MardynMolecule *temp = NULL;
      bool moleculeFound = false;
      const tarch::la::Vector<3, double> myPosition(_myMolecule->r(0), _myMolecule->r(1), _myMolecule->r(2));

      // get size of domain in number of cells
      LinkedCells *lc = (LinkedCells *)mdSim->getMolecules();
      int *boxWidthInNumCells = lc->boxWidthInNumCells();

      // compute linked cell index vector
      tarch::la::Vector<3, double> cellLength(0.0);
      tarch::la::Vector<3, double> haloBoundingBoxMin(0.0);
      tarch::la::Vector<3, unsigned int> cellIndexVector(0.0);
      LinkedCells *c = (LinkedCells *)mdSim->getMolecules();
      for (int d = 0; d < 3; d++) {
        cellLength[d] = (c->cellLength())[d];
        haloBoundingBoxMin[d] = c->getBoundingBoxMin(d) - c->get_halo_L(d);
      }
      for (int d = 0; d < 3; d++)
        cellIndexVector[d] = (int)floor((myPosition[d] - haloBoundingBoxMin[d]) / cellLength[d]);

      // compute the cell index (cellIndex = (z * cellsY + y) * cellsX + x)
      cellIndex = (cellIndexVector[2] * (boxWidthInNumCells[1] + 2 * (int)lc->getHaloWidthNumCells()) + cellIndexVector[1]) *
                      (boxWidthInNumCells[0] + 2 * (int)lc->getHaloWidthNumCells()) +
                  cellIndexVector[0];

      // get cell in which the molecule is located
      ParticleCell pc = lc->getCell(cellIndex);
      // create MarDynCell for iteration of molecules in the particle cell
      MarDynCell marDynCell = MarDynCell(&pc, _cutoffRadius);

      // iterate through molecules in cell until right molecule is found
      for (marDynCell.begin(); marDynCell.continueIteration(); marDynCell.next()) {
        temp = marDynCell.get();
        if (temp->r(0) == myPosition[0] && temp->r(1) == myPosition[1] && temp->r(2) == myPosition[2]) {
          moleculeFound = true;
          break;
        }
      }
      if (!moleculeFound) {
        global_log->debug() << "MarDynMoleculeWrapper::setPotentialEnergy(): "
                               "Molecule not found, must be new, creating..."
                            << std::endl;
        MardynMolecule temp2(global_simulation->getMaxID() + 1, global_simulation->getEnsemble()->component(0), _myMolecule->r(0), _myMolecule->r(1),
                             _myMolecule->r(2), _myMolecule->v(0), _myMolecule->v(1), _myMolecule->v(2), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        potentialEnergy = mdSim->getMolecules()->getEnergy(mdSim->getParticlePairsHandler(), &temp2, legacyCellProcessor);
      } else {
        potentialEnergy = mdSim->getMolecules()->getEnergy(mdSim->getParticlePairsHandler(), temp, legacyCellProcessor);
      }
    } else { // the molecule was initialized using a pointer to the molecule in
             // the simulation, no search needed
      potentialEnergy = mdSim->getMolecules()->getEnergy(mdSim->getParticlePairsHandler(), _myMolecule, legacyCellProcessor);
    }

    return potentialEnergy;
  }
  virtual void setPotentialEnergy(const double &potentialEnergy) {
    // currently not needed, therefore not implemented
  }

private:
  MardynMolecule *_myMolecule;
  double _cutoffRadius;
};

#endif /* MARDYNMOLECULEWRAPPER_H_ */
