// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_INTERFACE_SIMPLEMDSOLVERINTERFACE_CPP_
#define _MOLECULARDYNAMICS_COUPLING_INTERFACE_SIMPLEMDSOLVERINTERFACE_CPP_
#include "simplemd/MolecularDynamicsDefinitions.h"
#include "simplemd/services/LinkedCellService.h"
#include "simplemd/services/ParallelTopologyService.h"
#include "simplemd/services/MolecularPropertiesService.h"
#include "simplemd/services/ExternalForceService.h"
#include "simplemd/cell-mappings/LennardJonesPotentialEnergyMapping.h"
#include "simplemd/cell-mappings/LennardJonesForceMapping.h"
#include "simplemd/cell-mappings/ResetPotentialEnergyMapping.h"
#include "simplemd/BoundaryTreatment.h"

#include "coupling/CouplingMDDefinitions.h"
#include "coupling/interface/MDSolverInterface.h"
#include "coupling/interface/impl/SimpleMD/SimpleMDMoleculeIterator.h"

namespace coupling {
namespace interface {

/** general MD solver interface for SimpleMD.
 *  @author Philipp Neumann
 */
class SimpleMDSolverInterface
    : public MDSolverInterface<simplemd::LinkedCell, MD_DIM> {
private:
  simplemd::services::ParallelTopologyService &_parallelTopologyService;
  simplemd::services::MoleculeService &_moleculeService;
  simplemd::services::LinkedCellService &_linkedCellService;
  const simplemd::services::MolecularPropertiesService &
      _molecularPropertiesService;
  /** used for the synchronization of molecules in boundary regions and between
   * different processes when
   *  mass was inserted/ deleted by the coupling.
   */
  simplemd::BoundaryTreatment &_boundaryTreatment;
  const tarch::la::Vector<MD_LINKED_CELL_NEIGHBOURS, simplemd::BoundaryType> &
      _localBoundaryInformation;

  /** number of linked cells in each macroscopic cell */
  const tarch::la::Vector<MD_DIM, unsigned int>
      _numberOfLinkedCellsPerMacroscopicCell;

  const double _sigma6;
  const double _epsilon;
  const double _cutoffRadiusSquared;
  const double _cutoffEnergy;

  const double _dt;

  tarch::la::Vector<MD_DIM, unsigned int>
  computeNumberOfLinkedCellsPerMacroscopicCells(const tarch::la::Vector<
      MD_DIM, unsigned int> &localNumberOfMacroscopicCells) const {
    tarch::la::Vector<MD_DIM, unsigned int> numberOfLinkedCells(0);
    for (unsigned int d = 0; d < MD_DIM; d++) {
      numberOfLinkedCells[d] = _linkedCellService.getLocalNumberOfCells()[d] /
                               localNumberOfMacroscopicCells[d];
#if (COUPLING_MD_ERROR == COUPLING_MD_YES)
      if (numberOfLinkedCells[d] * localNumberOfMacroscopicCells[d] !=
          _linkedCellService.getLocalNumberOfCells()[d]) {
        std::cout
            << "ERROR coupling::interface::SimpleMDSolverInterface: Number of "
               "linked cells is not a multiple of the macroscopic cells!"
            << std::endl;
        std::cout << "Number linked cells: " << numberOfLinkedCells
                  << ", macroscopic cells: " << localNumberOfMacroscopicCells
                  << std::endl;
        std::cout << "Local number of linked cells: "
                  << _linkedCellService.getLocalNumberOfCells() << std::endl;
        exit(EXIT_FAILURE);
      }
#endif
    }
    return numberOfLinkedCells;
  }

  /** Though this method is available in LennardJonesPotentialEnergyMapping, I
   * re-implemented it at
   *  this stage, as the latter always fetches the epsilon,sigma-parameters
   * first. This is nice for
   *  simulations where parameters are changed over time; however for static
   * simulations it just slows
   *  down...
   */
  double
  getPotentialEnergy(const tarch::la::Vector<MD_DIM, double> &position1,
                     const tarch::la::Vector<MD_DIM, double> &position2) const {
    const tarch::la::Vector<MD_DIM, double> distance(position2 - position1);
    const double r2 = tarch::la::dot(distance, distance);
    // if we are in the near distance, do computation, otherwise return 0.0
    if (r2 <= _cutoffRadiusSquared) {
      const double r6 = r2 * r2 * r2;
      const double energyBuffer = 4.0 * _epsilon * (_sigma6 / r6) *
                                      ((_sigma6 / r6) - 1.0) - _cutoffEnergy;
      return 0.5 * energyBuffer;
    } else {
      return 0.0;
    }
  }

  tarch::la::Vector<MD_DIM, double> getLennardJonesForce(
      const tarch::la::Vector<MD_DIM, double> &position1,
      const tarch::la::Vector<MD_DIM, double> &position2) const {
    const tarch::la::Vector<MD_DIM, double> rij(position2 - position1);
    const double rij2 = tarch::la::dot(rij, rij);
#if (MD_ERROR == MD_YES)
    if (rij2 == 0.0) {
      std::cout << "ERROR "
                   "cellmappings::LennardJonesForceMapping::getLennardJonesForc"
                   "e(): Particle positions are identical!" << std::endl;
      exit(EXIT_FAILURE);
    }
#endif

    if (rij2 <= _cutoffRadiusSquared) {
      const double rij6 = rij2 * rij2 * rij2;
      return 24.0 * _epsilon / rij2 * (_sigma6 / rij6) *
             (1.0 - 2.0 * (_sigma6 / rij6)) * rij;
    } else {
      return tarch::la::Vector<MD_DIM, double>(0.0);
    }
  }

public:
  SimpleMDSolverInterface(
      tarch::la::Vector<MD_DIM, unsigned int>
          numberLinkedCellsPerMacroscopicCell,
      simplemd::BoundaryTreatment &boundaryTreatment,
      simplemd::services::ParallelTopologyService &parallelTopologyService,
      simplemd::services::MoleculeService &moleculeService,
      simplemd::services::LinkedCellService &linkedCellService,
      const simplemd::services::MolecularPropertiesService &
          molecularPropertiesService,
      const tarch::la::Vector<MD_LINKED_CELL_NEIGHBOURS,
                              simplemd::BoundaryType> &localBoundaryInformation,
      const double &dt)
      : _parallelTopologyService(parallelTopologyService),
        _moleculeService(moleculeService),
        _linkedCellService(linkedCellService),
        _molecularPropertiesService(molecularPropertiesService),
        _boundaryTreatment(boundaryTreatment),
        _localBoundaryInformation(localBoundaryInformation),
        _numberOfLinkedCellsPerMacroscopicCell(
            numberLinkedCellsPerMacroscopicCell),
        _sigma6(molecularPropertiesService.getMolecularProperties().getSigma() *
                molecularPropertiesService.getMolecularProperties().getSigma() *
                molecularPropertiesService.getMolecularProperties().getSigma() *
                molecularPropertiesService.getMolecularProperties().getSigma() *
                molecularPropertiesService.getMolecularProperties().getSigma() *
                molecularPropertiesService.getMolecularProperties().getSigma()),
        _epsilon(
            molecularPropertiesService.getMolecularProperties().getEpsilon()),
        _cutoffRadiusSquared(molecularPropertiesService.getMolecularProperties()
                                 .getCutOffRadius() * molecularPropertiesService
                                 .getMolecularProperties().getCutOffRadius()),
        _cutoffEnergy(4.0 * _epsilon * _sigma6 /
                      (_cutoffRadiusSquared * _cutoffRadiusSquared *
                       _cutoffRadiusSquared) *
                      (_sigma6 / (_cutoffRadiusSquared * _cutoffRadiusSquared *
                                  _cutoffRadiusSquared) - 1.0)),
        _dt(dt) {}
  ~SimpleMDSolverInterface() {}

  simplemd::LinkedCell &getLinkedCell(
      const tarch::la::Vector<MD_DIM, unsigned int> &macroscopicCellIndex,
      const tarch::la::Vector<MD_DIM, unsigned int> &
          linkedCellInMacroscopicCell,
      const tarch::la::Vector<MD_DIM, unsigned int> &
          linkedCellsPerMacroscopicCell,
      const coupling::IndexConversion<MD_DIM> &indexConversion) {
    // no linked cells found in outer region!
    for (unsigned int d = 0; d < MD_DIM; d++) {
      if (macroscopicCellIndex[d] == 0) {
        std::cout << "ERROR SimpleMDSolverInterface::getLinkedCell(): "
                     "macroscopic cell outside range for linked cells!"
                  << std::endl;
        exit(EXIT_FAILURE);
      }
    }
    tarch::la::Vector<MD_DIM, unsigned int> index(
        _linkedCellService.getLocalIndexOfFirstCell());
    for (unsigned int d = 0; d < MD_DIM; d++) {
      index[d] = index[d] + (macroscopicCellIndex[d] - 1) *
                                linkedCellsPerMacroscopicCell[d] +
                 linkedCellInMacroscopicCell[d];
    }
    return _linkedCellService.getLinkedCell(index);
  }

  /** returns the global size of the box-shaped MD domain */
  tarch::la::Vector<MD_DIM, double> getGlobalMDDomainSize() const {
    return _parallelTopologyService.getGlobalDomainSize();
  }

  /** returns the offset (i.e. lower,left corner) of MD domain */
  tarch::la::Vector<MD_DIM, double> getGlobalMDDomainOffset() const {
    return _parallelTopologyService.getGlobalDomainOffset();
  }

  double getMoleculeMass() const {
    return _molecularPropertiesService.getMolecularProperties().getMass();
  }

  double getKB() const {
    return _molecularPropertiesService.getMolecularProperties().getKB();
  }

  double getMoleculeSigma() const {
    return _molecularPropertiesService.getMolecularProperties().getSigma();
  }

  double getMoleculeEpsilon() const {
    return _molecularPropertiesService.getMolecularProperties().getEpsilon();
  }

  void
  getInitialVelocity(const tarch::la::Vector<MD_DIM, double> &meanVelocity,
                     const double &kB, const double &temperature,
                     tarch::la::Vector<MD_DIM, double> &initialVelocity) const {
    _moleculeService.getInitialVelocity(meanVelocity, kB, temperature,
                                        _molecularPropertiesService,
                                        initialVelocity);
  }

  void deleteMoleculeFromMDSimulation(
      const coupling::interface::Molecule<MD_DIM> &molecule,
      simplemd::LinkedCell &cell) {
    std::list<simplemd::Molecule *>::iterator it = cell.begin();
    const tarch::la::Vector<MD_DIM, double> moleculePosition =
        molecule.getPosition();
    while (it != cell.end()) {
      const tarch::la::Vector<MD_DIM, double> &itPosition(
          (*it)->getConstPosition());
      if (moleculePosition == itPosition) {
        simplemd::Molecule *myMolecule = (*it);
        cell.deleteMolecule(myMolecule);
        _moleculeService.deleteMolecule(*myMolecule);
        return;
      }
      it++;
    }

    std::cout << "Could delete molecule at position " << moleculePosition << "!"
              << std::endl;
    exit(EXIT_FAILURE);
  }

  void addMoleculeToMDSimulation(
      const coupling::interface::Molecule<MD_DIM> &molecule) {
    const tarch::la::Vector<MD_DIM, double> position = molecule.getPosition();
    const tarch::la::Vector<MD_DIM, double> velocity = molecule.getVelocity();
    const tarch::la::Vector<MD_DIM, double> force = molecule.getForce();
    const double potentialEnergy = molecule.getPotentialEnergy();
    simplemd::Molecule newMolecule(position, velocity);
    newMolecule.setForce(force);
    newMolecule.setPotentialEnergy(potentialEnergy);
    // add molecule to MoleculeService and LinkedCellService
    simplemd::Molecule *myMolecule = _moleculeService.addMolecule(newMolecule);
    tarch::la::Vector<MD_DIM, unsigned int> linkedCellIndex(
        getLinkedCellIndexForMoleculePosition(myMolecule->getConstPosition()));
    _linkedCellService.addMoleculeToLinkedCell(*myMolecule, linkedCellIndex);
  }

  tarch::la::Vector<MD_DIM, unsigned int> getLinkedCellIndexForMoleculePosition(
      const tarch::la::Vector<MD_DIM, double> &position) {
    tarch::la::Vector<MD_DIM, double> domainOffset =
        _parallelTopologyService.getGlobalDomainOffset();
    tarch::la::Vector<MD_DIM, double> meshWidth =
        _parallelTopologyService.getMeshWidth();
    tarch::la::Vector<MD_DIM, unsigned int> globalIndexOfFirstCell =
        _parallelTopologyService.getGlobalIndexOfFirstCell();
    tarch::la::Vector<MD_DIM, unsigned int> localIndexOfFirstCell =
        _linkedCellService.getLocalIndexOfFirstCell();

#if (COUPLING_MD_ERROR == COUPLING_MD_YES)
    const tarch::la::Vector<MD_DIM, double> domainSize =
        _parallelTopologyService.getGlobalDomainSize();
    for (unsigned int d = 0; d < MD_DIM; d++) {
      if ((position[d] < domainOffset[d] - meshWidth[d]) ||
          (position[d] > domainOffset[d] + domainSize[d] + meshWidth[d])) {
        std::cout << "ERROR "
                     "coupling::interface::impl::MDSolverInterface::addMolecule"
                     "ToLinkedCell: Position ";
        std::cout << d << " is out of range!" << std::endl;
        std::cout << "Position: " << position << std::endl;
        exit(EXIT_FAILURE);
      }
    }
#endif
    tarch::la::Vector<MD_DIM, unsigned int> cellVectorIndex(0);

    // determine current cell index (in serial, i.e. 1-D, form)
    for (unsigned int d = 0; d < MD_DIM; d++) {
      // find global cell index
      int index = static_cast<int>(
          floor((position[d] - domainOffset[d]) / meshWidth[d]));
      // shift into local cell index
      index = index - globalIndexOfFirstCell[d] + localIndexOfFirstCell[d];
#if (COUPLING_MD_ERROR == COUPLING_MD_YES)
      if (index < 0) {
        std::cout << "ERROR "
                     "coupling::interfaces::impl::SimpleMD::SimpleMDSolverInter"
                     "face: index < 0: index=";
        std::cout << index << std::endl;
        exit(EXIT_FAILURE);
      }
#endif
      cellVectorIndex[d] = static_cast<unsigned int>(index);
    }

    return cellVectorIndex;
  }

  void setupPotentialEnergyLandscape(
      const tarch::la::Vector<MD_DIM, unsigned int> &
          indexOfFirstMacroscopicCell,
      const tarch::la::Vector<MD_DIM, unsigned int> &rangeMacroscopicCells,
      const tarch::la::Vector<MD_DIM, unsigned int> &
          linkedCellsPerMacroscopicCell) {
    simplemd::cellmappings::ResetPotentialEnergyMapping
        resetPotentialEnergyMapping;
    simplemd::cellmappings::LennardJonesPotentialEnergyMapping
        potentialEnergyMapping(_molecularPropertiesService);
    tarch::la::Vector<MD_DIM, unsigned int> rangeLinkedCellsExtended(0);
    tarch::la::Vector<MD_DIM, unsigned int> firstLinkedCell(0);
    // compute coordinates of the first linked cell and the range of linked
    // cells to be considered
    for (unsigned int d = 0; d < MD_DIM; d++) {
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
      if (indexOfFirstMacroscopicCell[d] == 0) {
        std::cout << "ERROR setupPotentialEnergyLandscape: "
                  << indexOfFirstMacroscopicCell[d] << std::endl;
        exit(EXIT_FAILURE);
      }
#endif

      firstLinkedCell[d] = (indexOfFirstMacroscopicCell[d] - 1) *
                           linkedCellsPerMacroscopicCell[d];
      // we need to loop over one more layer of linked cells since
      // firstLinkedCell(d) already starts one linked cell earlier
      // (e.g. already in linked cell ghost layer)
      rangeLinkedCellsExtended[d] =
          rangeMacroscopicCells[d] * linkedCellsPerMacroscopicCell[d] +
          _linkedCellService.getLocalIndexOfFirstCell()[d];
    }

    // reset potential energy first for the molecules in all relevant linked
    // cells
    _linkedCellService.iterateCells(resetPotentialEnergyMapping,
                                    firstLinkedCell, rangeLinkedCellsExtended,
                                    false);

    // compute potential energy in these cells
    _linkedCellService.iterateCellPairs(potentialEnergyMapping, firstLinkedCell,
                                        rangeLinkedCellsExtended, false);
  }

  void
  calculateForceAndEnergy(coupling::interface::Molecule<MD_DIM> &molecule) {
    const tarch::la::Vector<MD_DIM, double> position = molecule.getPosition();
    const tarch::la::Vector<MD_DIM, unsigned int> linkedCellIndex =
        getLinkedCellIndexForMoleculePosition(position);
#if (COUPLING_MD_ERROR == COUPLING_MD_YES)
    for (unsigned int d = 0; d < MD_DIM; d++) {
      if ((linkedCellIndex[d] < 1) ||
          (linkedCellIndex[d] >
           _linkedCellService.getLocalNumberOfCells()[d])) {
        std::cout << "ERROR "
                     "coupling::interface/impl/SimpleMD/"
                     "SimpleMDSolverInterface::calculateForceAndEnergy(): "
                     "linkedCellIndex out of range!" << std::endl;
        std::cout << "LinkedCellIndex=" << linkedCellIndex << std::endl;
        exit(EXIT_FAILURE);
      }
    }
#endif
    tarch::la::Vector<MD_DIM, double> forceOld(0.0);
    double energy = 0.0;
    molecule.setPotentialEnergy(0.0);

    tarch::la::Vector<MD_DIM, unsigned int> loopIndex(0);

// loop over all relevant linked cells and compute force and energy
// contributions
#if (MD_DIM > 2)
    for (loopIndex[2] = linkedCellIndex[2] - 1;
         loopIndex[2] < linkedCellIndex[2] + 2; loopIndex[2]++) {
#endif
#if (MD_DIM > 1)
      for (loopIndex[1] = linkedCellIndex[1] - 1;
           loopIndex[1] < linkedCellIndex[1] + 2; loopIndex[1]++) {
#endif
        for (loopIndex[0] = linkedCellIndex[0] - 1;
             loopIndex[0] < linkedCellIndex[0] + 2; loopIndex[0]++) {

          // loop over all molecules in each cell
          const simplemd::LinkedCell &thisCell =
              _linkedCellService.getLinkedCell(loopIndex);
          for (std::list<simplemd::Molecule *>::const_iterator it =
                   thisCell.constBegin();
               it != thisCell.constEnd(); it++) {
            forceOld +=
                getLennardJonesForce(position, (*it)->getConstPosition());
            energy += getPotentialEnergy(position, (*it)->getConstPosition());
          }
        }
#if (MD_DIM > 1)
      }
#endif
#if (MD_DIM > 2)
    }
#endif
    // set force in forceOld-entry
    molecule.setForce(forceOld);
    // set total potential energy
    molecule.setPotentialEnergy(energy);
  }

  void synchronizeMoleculesAfterMassModification() {
    _boundaryTreatment.emptyGhostBoundaryCells();
#if (TEST_TCHIPEV == MD_NO)
    _boundaryTreatment.fillBoundaryCells(_localBoundaryInformation,
                                         _parallelTopologyService);
#else
    _boundaryTreatment.putBoundaryParticlesToInnerCellsAndFillBoundaryCells(
        _localBoundaryInformation, _parallelTopologyService);
#endif
  }

  void synchronizeMoleculesAfterMomentumModification() {}

  double getDt() { return _dt; }

  coupling::interface::MoleculeIterator<simplemd::LinkedCell, MD_DIM> *
  getMoleculeIterator(simplemd::LinkedCell &cell) {
    return new coupling::interface::SimpleMDMoleculeIterator(cell);
  }
};

}
}
#endif // _MOLECULARDYNAMICS_COUPLING_INTERFACE_SIMPLEMDSOLVERINTERFACE_CPP_
