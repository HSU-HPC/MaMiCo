// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef MARDYNMDSOLVERINTERFACE_H_
#define MARDYNMDSOLVERINTERFACE_H_

// Mardyn
#include "Domain.h"
#include "Simulation.h"
#include "molecules/Molecule.h"
#include "particleContainer/AdaptiveSubCells.h"
#include "particleContainer/LinkedCells.h"
#include "particleContainer/ParticleContainer.h"
#include "particleContainer/adapter/LegacyCellProcessor.h"
typedef Molecule MardynMolecule;
#include "parallel/DomainDecompBase.h"

// MaMiCo
#include "coupling/CouplingMDDefinitions.h"
#include "coupling/IndexConversion.h"
#include "coupling/interface/MDSolverInterface.h"
#include "coupling/interface/impl/MarDyn/MarDynCell.h"
#include "coupling/interface/impl/MarDyn/MarDynCoupledSimulation.h"
#include "coupling/interface/impl/MarDyn/MarDynMoleculeIterator.h"

#include "tarch/la/ScalarOperations.h"
#include "tarch/la/Vector.h"
#include <cmath>
#include <iomanip>
#include <stdio.h>
#include <vector>

/*
 * 	Implementation of the MD solver interface for Mardyn
 * 	@author Hanno Flohr, Philipp Neumann
 */
class MarDynMDSolverInterface : public coupling::interface::MDSolverInterface<MarDynCell, 3> {

public:
  MarDynMDSolverInterface(MarDynCoupledSimulation *simulation)
      : coupling::interface::MDSolverInterface<MarDynCell, 3>(), _mySimulation(simulation), _cutoffRadius(_mySimulation->getLJCutoff()), _tolerance(1.0e-8) {
    if (_mySimulation->getContainerType() == 0) {
      LinkedCells *lc = (LinkedCells *)_mySimulation->getMolecules();
      _haloWidthInNumCells = 2 * (int)(lc->getHaloWidthNumCells());
    } else {
      std::cout << "MarDynMDSolver Constructor - ERROR: unsupported container type!" << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  virtual ~MarDynMDSolverInterface() {
    // clear the MarDyn cell storage
    for (unsigned int i = 0; i < _marDynCellPtrStorage.size(); i++) {
      if (_marDynCellPtrStorage[i] != NULL) {
        delete _marDynCellPtrStorage[i];
        _marDynCellPtrStorage[i] = NULL;
      }
    }
    _marDynCellPtrStorage.clear();
  }

  /** returns a particular linked cell inside a macroscopic cell.
   *  The macroscopic cells are currently located on the same process as the
   * respective linked cells. However, several linked cells may be part of a
   * macroscopic cell. The macroscopic cells also contain a ghost layer which
   * surrounds each local domain; the very first macroscopic cell inside the
   * global MD domain (or local MD domain) is thus given by coordinates (1,1,1)
   * (or (1,1) in 2D, respectively). The index linkedCellInMacroscopicCell
   * corresponds to the coordinates of the linked cell inside the given
   * macroscopic cell. These coordinates thus lie in a range
   * (0,linkedCellsPerMacroscopicCell-1).
   */
  virtual MarDynCell &getLinkedCell(const tarch::la::Vector<3, unsigned int> &macroscopicCellIndex,
                                    const tarch::la::Vector<3, unsigned int> &linkedCellInMacroscopic,
                                    const tarch::la::Vector<3, unsigned int> &linkedCellsPerMacroscopicCell,
                                    const coupling::IndexConversion<3> &indexConversion) {
    // no linked cells found in outer region
    for (unsigned int d = 0; d < 3; d++) {
      if (macroscopicCellIndex[d] == 0) {
        std::cout << "ERROR in MarDynMDSolverInterface::getLinkedCell(): "
                     "macroscopic cell index out of range for linked cells!"
                  << std::endl;
        exit(EXIT_FAILURE);
      }
    }

    MarDynCell *cell;
    if (_mySimulation->getContainerType() == 0) {
      LinkedCellsForCoupling *cells = (LinkedCellsForCoupling *)_mySimulation->getMolecules();

      // size of the macroscopic cells
      tarch::la::Vector<3, double> macroCellSize = indexConversion.getMacroscopicCellSize();
      // size of the md cells
      double *mdCellsize = cells->cellLength();
      // the requested 3D cell index vector
      tarch::la::Vector<3, unsigned int> requestedCellIndex(0);
      // number of cells in each dimension
      int *boxWidthInNumCells = cells->boxWidthInNumCells();

      // compute requested cell index vector
      for (int d = 0; d < 3; d++) {
        requestedCellIndex[d] = (int)floor(macroscopicCellIndex[d] * macroCellSize[d] / mdCellsize[d]);
        // adjust requested cell index if multiple linked cells per macroscopic
        // cell are used
        if (linkedCellsPerMacroscopicCell[d] > 1)
          requestedCellIndex[d] -= (linkedCellsPerMacroscopicCell[d] - 1) - linkedCellInMacroscopic[d];
      }
      // compute the Mardyn cell index of the requested cell
      unsigned int cellIndex =
          (requestedCellIndex[2] * (boxWidthInNumCells[1] + _haloWidthInNumCells) + requestedCellIndex[1]) * (boxWidthInNumCells[0] + _haloWidthInNumCells) +
          requestedCellIndex[0];
      // get the requested cell for the computed cell index
      cell = new MarDynCell(cells->getCellPointer(cellIndex), cells->getCutoff());
      // store the MarDynCell pointer in the pointer storage
      _marDynCellPtrStorage.push_back(cell);
    } else {
      std::cout << "ERROR (MarDynSolverInterface): Container type is not set "
                   "or supported! ("
                << _mySimulation->getContainerType() << ")" << std::endl;
      exit(1);
    }

    return *cell;
  }

  /* returns the global size of the box-shaped MD domain */
  virtual tarch::la::Vector<3, double> getGlobalMDDomainSize() const {
    tarch::la::Vector<3, double> globalMDsize(0.0);
    for (int i = 0; i < 3; i++) {
      globalMDsize[i] = _mySimulation->getEnsemble()->domain()->rmax()[i] - _mySimulation->getEnsemble()->domain()->rmin()[i];
    }
    return globalMDsize;
  }

  /* returns the offset (i.e. lower, left corner) of the MD domain */
  virtual tarch::la::Vector<3, double> getGlobalMDDomainOffset() const {
    tarch::la::Vector<3, double> globalMDoffset(0.0);
    for (int i = 0; i < 3; i++) {
      globalMDoffset[i] = _mySimulation->getEnsemble()->domain()->rmin()[i];
    }
    return globalMDoffset;
  }

  /* returns the mass of a single fluid molecule
   * the coupling right now is restricted to molecules of the same type
   * thus: component 0 contains the mass of all molecules
   */
  virtual double getMoleculeMass() const { return _mySimulation->getEnsemble()->component(0)->m(); }

  /* returns Boltzmann's constant
   * The Boltzmann's constant is not stored in MarDyn. When needed, the
   * value 1.0 is used. Therefore this method just returns 1.0 for now.
   */
  virtual double getKB() const { return 1.0; }

  /* returns the sigma parameter of the LJ potential */
  virtual double getMoleculeSigma() const { return _mySimulation->getEnsemble()->component(0)->getSigma(0); }

  /* returns the epsilon parameter of the LJ potential */
  virtual double getMoleculeEpsilon() const { return _mySimulation->getEnsemble()->component(0)->ljcenter(0).eps(); }

  /* sets a random velocity in the vector 'initialVelocity' */
  virtual void getInitialVelocity(const tarch::la::Vector<3, double> &meanVelocity, const double &kB, const double &temperature,
                                  tarch::la::Vector<3, double> &initialVelocity) const {
    // temperature based standard deviation of gaussian distribution
    const double standardDeviation = std::sqrt(kB * 3 * temperature / getMoleculeMass());

    // set a random number
    tarch::la::Vector<3, double> random(0.0);
    random[0] = tarch::utils::RandomNumberService::getInstance().getGaussianRandomNumber();
    for (unsigned int d = 1; d < 3; d++)
      random[d] = 2.0 * M_PI * tarch::utils::RandomNumberService::getInstance().getUniformRandomNumber();

    // set initial velocity with randomized values
    initialVelocity[0] = meanVelocity[0] + standardDeviation * (random[0] * std::sin(random[1]) * std::cos(random[2]));
    initialVelocity[1] = meanVelocity[1] + standardDeviation * (random[0] * std::sin(random[1]) * std::sin(random[2]));
    initialVelocity[2] = meanVelocity[2] + standardDeviation * (random[0] * std::cos(random[1]));
  }

  /* deletes the molecule from the MD simulation */
  virtual void deleteMoleculeFromMDSimulation(const coupling::interface::Molecule<3> &molecule, MarDynCell &cell) {
    // molecule position
    const tarch::la::Vector<3, double> position = molecule.getPosition();

    MardynMolecule *curMolecule;
    bool success = false;
    bool moleculeFound = true;
    tarch::la::Vector<3, double> curMoleculePosition(0.0);

    // loop through all molecules in the cell to find the desired molecule
    for (cell.begin(); cell.continueIteration(); cell.next()) {
      // set position of the current molecule
      curMolecule = cell.get();
      for (int d = 0; d < 3; d++)
        curMoleculePosition[d] = curMolecule->r(d);

      // compare positions
      for (unsigned int d = 0; d < 3; d++) {
        moleculeFound = moleculeFound && (tarch::la::equals(position[d], curMoleculePosition[d], _tolerance));
      }

      // if molecule is found: break
      if (moleculeFound) {
        success = true;
        break;
      } else
        moleculeFound = true;
    }

    // if a molecule to delete was found delete it, else through error message
    // and continue
    if (success) {
      std::cout << "deleting molecule with id: " << curMolecule->id() << std::endl;
      // delete molecule from the particle cell and the particle container
      cell.getParticleCell()->deleteMolecule(curMolecule->id());
      _mySimulation->removeMoleculeFromContainer(curMolecule->id());
    } else
      std::cout << "ERROR (delete molecule): no molecule found at desired "
                   "delete position!"
                << std::endl;
  }

  /* adds the molecule to the MD simulation */
  virtual void addMoleculeToMDSimulation(const coupling::interface::Molecule<3> &molecule) {
    // create a MarDynMolecule and set molecule id to maxiID+1
    MardynMolecule mdmw;
    const tarch::la::Vector<3, double> pos(molecule.getPosition());
    const tarch::la::Vector<3, double> vel(molecule.getVelocity());
    const tarch::la::Vector<3, double> f(molecule.getForce());
    double force[3] = {f[0], f[1], f[2]};

    mdmw.setid(_mySimulation->getMaxID() + 1);
    for (unsigned int d = 0; d < 3; d++) {
      mdmw.setr(d, pos[d]);
      mdmw.setv(d, vel[d]);
    }
    mdmw.setF(force);

    // add molecule to the particle container
    _mySimulation->getMolecules()->addParticle(mdmw);
  }

  /* sets up the potential energy landscape */
  virtual void setupPotentialEnergyLandscape(const tarch::la::Vector<3, unsigned int> &indexOfFirstMacroscopicCell,
                                             const tarch::la::Vector<3, unsigned int> &rangeMacroscopicCells,
                                             const tarch::la::Vector<3, unsigned int> &linkedCellsPerMacroscopicCell) {
    // The potential energy is calculated for each molecule individually using
    // the getPotentialEnergy() method of the MoleculeWrapper interface.
    // Therefore this method does nothing for now.
  }

  /* returns the local index vector */
  virtual tarch::la::Vector<3, unsigned int> getLinkedCellIndexForMoleculePosition(const tarch::la::Vector<3, double> &position) {
    tarch::la::Vector<3, double> cellLength(0.0);
    tarch::la::Vector<3, double> haloBoundingBoxMin(0.0);
    tarch::la::Vector<3, unsigned int> cellIndex(0);

    // get cell length and halo bounding box values for all dimensions
    if (_mySimulation->getContainerType() == 0) {
      LinkedCells *cells = (LinkedCells *)_mySimulation->getMolecules();
      for (int d = 0; d < 3; d++) {
        cellLength[d] = (cells->cellLength())[d];
        haloBoundingBoxMin[d] = cells->getBoundingBoxMin(d) - cells->get_halo_L(d);
      }
    } else {
      std::cout << "ERROR (MarDynSolverInterface): Container type is not "
                   "supported or unknown! ("
                << _mySimulation->getContainerType() << ")" << std::endl;
      exit(1);
    }

    // compute cell index for each dimension
    for (int d = 0; d < 3; d++)
      cellIndex[d] = (int)floor((position[d] - haloBoundingBoxMin[d]) / cellLength[d]);

    return cellIndex;
  }

  /** calculates force and energy for a molecule */
  virtual void calculateForceAndEnergy(coupling::interface::Molecule<3> &molecule) {
    // the requested force and energy
    tarch::la::Vector<3, double> force(0.0);
    double potentialEnergy = 0.0;
    // molecule position
    const tarch::la::Vector<3, double> moleculePosition = molecule.getPosition();
    bool moleculeFound = false;
    Molecule *moleculeInSim = NULL;
    // use a legacy cell processor  (the vectorized version seems not to work
    // for molecules that are not part of the sim)
    LegacyCellProcessor legacyCellProcessor(_cutoffRadius, _mySimulation->getLJCutoff(), _mySimulation->getTersoffCutoff(),
                                            _mySimulation->getParticlePairsHandler());

    if (_mySimulation->getContainerType() == 0) {
      // get cell index for the molecule
      LinkedCells *lc = (LinkedCells *)_mySimulation->getMolecules();
      Molecule refMolecule;
      for (int d = 0; d < 3; d++)
        refMolecule.setr(d, molecule.getPosition()[d]);
      long int cellIndex = lc->getCellIndexOfMolecule(&refMolecule);

      // get cell in which the molecule is located
      ParticleCell pc = lc->getCell(cellIndex);
      // create MarDynCell for iteration of molecules in the particle cell
      MarDynCell marDynCell = MarDynCell(&pc, _cutoffRadius);

      // iterate through molecules in cell until right molecule is found
      for (marDynCell.begin(); marDynCell.continueIteration(); marDynCell.next()) {
        moleculeInSim = marDynCell.get();
        // compare positions
        if (moleculeInSim->r(0) == moleculePosition[0] && moleculeInSim->r(1) == moleculePosition[1] && moleculeInSim->r(2) == moleculePosition[2]) {
          moleculeFound = true;
          break;
        }
      }
      // if no molecule at this position was found create new mardyn molecule
      // and compute force for it
      if (!moleculeFound) {
        global_log->debug() << "MarDynMDSolverInterface::calculateForceAndEnergy(): Molecule "
                               "not found, creating new molecule..."
                            << std::endl;
        // create new Mardyn molecule
        moleculeInSim = new Molecule(_mySimulation->getMaxID() + 1, _mySimulation->getEnsemble()->component(0), molecule.getPosition()[0],
                                     molecule.getPosition()[1], molecule.getPosition()[2], molecule.getVelocity()[0], molecule.getVelocity()[1],
                                     molecule.getVelocity()[2], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        // force calculation
        force = calculateForce(moleculeInSim, getLinkedCellIndexForMoleculePosition(moleculePosition));
      } else { // set force to the values of the molecule
        for (int d = 0; d < 3; d++)
          force[d] = moleculeInSim->F(d);
      }
    } else {
      std::cout << "ERROR (MarDynSolverInterface): Container type is not "
                   "supported or unknown! ("
                << _mySimulation->getContainerType() << ")" << std::endl;
      exit(1);
    }
    std::cout << "here5" << std::endl;
    // get potential energy of molecule
    potentialEnergy = _mySimulation->getMolecules()->getEnergy(_mySimulation->getParticlePairsHandler(), moleculeInSim, legacyCellProcessor);
    if (!moleculeFound) {
      delete moleculeInSim;
    }
    std::cout << "here6" << std::endl;
    // set potential energy and force of the molecule
    molecule.setPotentialEnergy(potentialEnergy);
    molecule.setForce(force);
  }

  /* called after insertion/deletion of molecules from the MD simulation by
   * MaMiCo */
  virtual void synchronizeMoleculesAfterMassModification() {
    /*	The following three steps are performed:
     * 	- clear the linked cells and fill them with updated molecule positions
     * 	- molecule exchange for the domain decomposition
     * 	- molecule cache is updated
     */
    _mySimulation->updateParticleContainerAndDecomposition();

    // clear halo cells in the MD simulation
    _mySimulation->getMolecules()->deleteOuterParticles();
  }

  /* called after momentum insertion into the MD simulation by MaMiCo */
  virtual void synchronizeMoleculesAfterMomentumModification() {
    // at the moment this method is not needed, therefore not implemented
  }

  /* returns the time step size of the MD simulation */
  virtual double getDt() { return _mySimulation->getTimestepLength(); }

  /* returns a new molecule iterator for a certain linked cell */
  virtual coupling::interface::MoleculeIterator<MarDynCell, 3> *getMoleculeIterator(MarDynCell &cell) {
    return new MarDynMoleculeIterator(cell);
  }

  // returns a pointer to the MarDynCoupledSimulation belonging to this
  // interface
  MarDynCoupledSimulation *
  getSimulation() {
    return _mySimulation;
  }

protected:
  // calculates the force of a Molecule that is not in the Mardyn simulation
  tarch::la::Vector<3, double> calculateForce(Molecule *molecule, tarch::la::Vector<3, unsigned int> cellIndexVector) {
    LinkedCellsForCoupling *cells = (LinkedCellsForCoupling *)_mySimulation->getMolecules();
    int *boxWidthInNumCells = cells->boxWidthInNumCells();
    tarch::la::Vector<3, double> force(0.0);
    tarch::la::Vector<3, double> position1(molecule->r(0), molecule->r(1), molecule->r(2));
    const double cutoffSquared = _cutoffRadius * _cutoffRadius;
    const double epsilon = getMoleculeEpsilon();
    const double sigma6 = pow(getMoleculeSigma(), 6.0);

    // initialize start and end coordinates for neighbor search
    tarch::la::Vector<3, unsigned int> start(0);
    tarch::la::Vector<3, unsigned int> end(1);
    tarch::la::Vector<3, unsigned int> loop(0);
    for (int d = 0; d < 3; d++) {
      start[d] = cellIndexVector[d] - 1;
      end[d] = start[d] + 3;
    }

    // loop through all neighboring cells
    for (loop[2] = start[2]; loop[2] < end[2]; loop[2]++) {
      for (loop[1] = start[1]; loop[1] < end[1]; loop[1]++) {
        for (loop[0] = start[0]; loop[0] < end[0]; loop[0]++) {
          unsigned int cellIndex =
              (loop[2] * (boxWidthInNumCells[1] + _haloWidthInNumCells) + loop[1]) * (boxWidthInNumCells[0] + _haloWidthInNumCells) + loop[0];
          MarDynCell cell(cells->getCellPointer(cellIndex), _cutoffRadius);
          coupling::interface::MoleculeIterator<MarDynCell, 3> *mdmi(getMoleculeIterator(cell));

          // loop through all molecules of the current cell
          for (mdmi->begin(); mdmi->continueIteration(); mdmi->next()) {
            const tarch::la::Vector<3, double> r = mdmi->getConst().getPosition() - position1;
            const double r2 = tarch::la::dot(r, r);

            // compute and add force contribution to overall force
            if (r2 <= cutoffSquared) {
              const double r6 = r2 * r2 * r2;
              force += 24.0 * epsilon / r2 * (sigma6 / r6) * (1.0 - 2.0 * (sigma6 / r6)) * r;
            }
          }
          delete mdmi;
        }
      }
    }
    return force;
  }

  // the simulation object
  MarDynCoupledSimulation *_mySimulation;
  // the amount of halo cells in each dimension
  int _haloWidthInNumCells;
  // the cutoff radius of the simulation
  const double _cutoffRadius;
  // the tolerance value for position comparisons
  const double _tolerance;

  /*	Used in the getLinkedCell() method to store the created MarDynCell
   * pointers to avoid memory leak the pointers are stored and this vector
   * cleared in the Destructor
   */
  std::vector<MarDynCell *> _marDynCellPtrStorage;
};

#endif /* MARDYNMDSOLVERINTERFACE_H_ */
