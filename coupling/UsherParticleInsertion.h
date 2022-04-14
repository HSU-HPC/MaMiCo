// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_USHERPARTICLEINSERTION_H_
#define _MOLECULARDYNAMICS_COUPLING_USHERPARTICLEINSERTION_H_

#include "coupling/CouplingMDDefinitions.h"
#include "coupling/ParticleInsertion.h"
#include "coupling/cell-mappings/ComputeMassMapping.h"
#include "coupling/cell-mappings/ComputeMomentumMapping.h"
#include "coupling/cell-mappings/ComputeTemperatureMapping.h"
#include "coupling/cell-mappings/DeleteParticleMapping.h"
#include "coupling/datastructures/MacroscopicCell.h"
#include "coupling/datastructures/Molecule.h"
#include "coupling/interface/MDSolverInterface.h"
#include "tarch/utils/RandomNumberService.h"

//#define USHER_DEBUG

namespace coupling {
template <class LinkedCell, unsigned int dim> class UsherParticleInsertion;
}

/** particles will be inserted or removed based on the microscopic mass of a
 * cell The algorithm is only applied to outer boundary cells of the md (non
 * ghost cells) We currently also allow particle insertion in energyy holes,
 * i.e. for energy(particle)==0.0.
 *  @brief handles particle insertion (via Usher algorithm) and random particle
 * deletion.
 *  @author Philipp Neumann
 *  @tparam LinkedCell the LinkedCell class is given by the implementation of linked cells in the molecular dynamics simulation
 *  @tparam dim  refers to the spacial dimension of the simulation, can be 1, 2, or 3 */
template <class LinkedCell, unsigned int dim> class coupling::UsherParticleInsertion : public coupling::ParticleInsertion<LinkedCell, dim> {
public:
  /** @brief a simple constructor
   *  @param insertDeleteMassAtTimestep the interval of time steps for the insertion/deletion of mass
   *  @param rSigmaCoeff coefficient to adapt r_sigma (see Usher algorithm for details,),
   *                     r_sigma = sigma * rSigmaCoeff
   *  @param meanPotentialEnergyFactor scales for the U_0 value of the Usher algorithm
   *  @param uOverlapCoeff scales the overlap energy U_ovlp (refers to the Usher algorithm)
   *  @param stepRefCoeff the maximal displacment within the usher method
   *  @param iterMax maximal number of usher iterations per try
   *  @param restartMax maximal number of usher restarts
   *  @param tolerance sets the tolerance, between the reached energy level and the goal energy level
   *  @param offsetFromOuterBoundary offset of the md domain (difference between the left, lower, front corner to 0,0,0)
   *  @param mdSolverInterface interface for the md solver */
  UsherParticleInsertion(unsigned int insertDeleteMassEveryTimestep, double rSigmaCoeff, double meanPotentialEnergyFactor, double uOverlapCoeff,
                         double stepRefCoeff, unsigned int iterMax, unsigned int restartMax, double tolerance, double offsetFromOuterBoundary,
                         coupling::interface::MDSolverInterface<LinkedCell, dim>* const mdSolverInterface);

  /** @brief a dummy destructor*/
  virtual ~UsherParticleInsertion() {}

  /** @todo why is this here? Already declared in the ParticleInsertion.*/
  // enum Action{NoAction=0,Insertion=1,Deletion=2};

  /** the function checks for the value of the microscopic mass and calls upon this
   *  either the inserteParticle() er deleteParticle() method, or just returns without an action
   *  @brief checks if a particle needs to be inserted or deleted in a cell
   *  @param cell the macroscopic cell to check if an action is necessary
   *  @param macroscopicCellPosition the postion of the macroscopic cell
   *  @param macroscopicCellSize the size of the macroscopic cell per direction
   *  @param meanVelocity the mean velocity of all particles in the cell
   *  @param temperature the temperature in the cell
   *  @param boundaryForceController an instance of the boundary force controller in application
   *  @returns the type of action which was applied (coupling::ParticleInsertion::Action) */
  virtual typename coupling::ParticleInsertion<LinkedCell, dim>::Action
  insertDeleteMass(coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell, dim>& cell,
                   const tarch::la::Vector<dim, double>& macroscopicCellPosition, const tarch::la::Vector<dim, double>& macroscopicCellSize,
                   const tarch::la::Vector<dim, double>& meanVelocity, const double& temperature,
                   const coupling::BoundaryForceController<LinkedCell, dim>& boundaryForceController);

  /** @brief since the Usher algorithm requires the potential energy landscape, the function returns true
   *  @returns true */
  virtual bool requiresPotentialEnergyLandscape() { return true; }

#ifdef USHER_DEBUG
  /** @brief the amount of energy inserted to the cell by the call to Usher */
  double _energyInserted;
  /** @brief the amount of energy removed of the cell by the call to Usher */
  double _energyRemoved;
  /** @brief the amount of Zhou energy inserted to the cell by the call to Usher */
  double _ZhouEnergyInserted;
  /** @brief the amount of Zhou energy removed of the cell by the call to Usher */
  double _ZhouEnergyRemoved;
  /** @brief the number of particles inserted into the cell*/
  int _particlesInserted;
  /** @brief the number of particles removed from the cell*/
  int _particlesRemoved;
#endif

  // private:
  /** Returns Insertion on success and NoAction otherwise.
   *  @brief tries to insert a particle in the respective macroscopic cell.
   *  @param cell the macroscopic cell to check if an action is necessary
   *  @param macroscopicCellPosition the postion of the macroscopic cell
   *  @param macroscopicCellSize the size of the macroscopic cell per direction
   *  @param meanVelocity the mean velocity of all particles in the cell
   *  @param temperature the temperature in the cell
   *  @param boundaryForceController an instance of the boundary force controller in application
   *  @returns NoAction or Insertion, depending on the applied action
   *  (coupling::ParticleInsertion::Action) */
  typename coupling::ParticleInsertion<LinkedCell, dim>::Action
  insertParticle(coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell, dim>& cell, const tarch::la::Vector<dim, double>& macroscopicCellPosition,
                 const tarch::la::Vector<dim, double>& macroscopicCellSize, const tarch::la::Vector<dim, double>& meanVelocity, const double& temperature,
                 const coupling::BoundaryForceController<LinkedCell, dim>& boundaryForceController) override;

private:
  /** Returns Deletion on success and NoAction otherwise.
   *  @brief tries to delete a particle from the macroscopic cell.
   *  @param cell the macroscopic cell to check if an action is necessary
   *  @param boundaryForceController an instance of the boundary force controller in application
   *  @returns NoAction or Deletion, depending on the action applied
   *  (coupling::ParticleInsertion::Action)*/
  typename coupling::ParticleInsertion<LinkedCell, dim>::Action
  deleteParticle(coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell, dim>& cell,
                 const coupling::BoundaryForceController<LinkedCell, dim>& boundaryForceController);

protected:
  /** the result is stored in the position entry of the molecule 'molecule'.
   *  The method returns Insertion if a position has been found and NoAction otherwise.
   *  @brief determines the position of a new particle within the macroscopic cell
   *  @param thisCell the macroscopic cell to check if an action is necessary
   *  @param macroscopicCellPosition the postion of the macroscopic cell
   *  @param macroscopicCellSize the size of the macroscopic cell per direction
   *  @param molecule the molecule which shall be inserted
   *  @param boundaryForceController an instance of the boundary force controller in application
   *  @returns NoAction or Insertion, depending on the result of the Usher method applied here
   *  (coupling::ParticleInsertion::Action)*/
  virtual typename coupling::ParticleInsertion<LinkedCell, dim>::Action
  findParticlePosition(coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell, dim>& thisCell,
                       const tarch::la::Vector<dim, double>& macroscopicCellPosition, const tarch::la::Vector<dim, double>& macroscopicCellSize,
                       coupling::datastructures::Molecule<dim>& molecule, const coupling::BoundaryForceController<LinkedCell, dim>& boundaryForceController);

  /** @brief collects all the parameters necessary for the Usher algorithm
   *  @author Philipp Neumann */
  class UsherParams {
  public:
    /** coefficient to adapt r_sigma (see Usher algorithm for details,),
     *                     r_sigma = sigma * rSigmaCoeff */
    double _rSigmaCoeff;
    /** scales for the U_0 value of the Usher algorithm */
    double _meanPotentialEnergyFactor;
    /** scales the overlap energy U_ovlp (refers to the Usher algorithm) */
    double _uOverlapCoeff;
    /** the maximal displacment */
    double _stepRefCoeff;
    /** maximal number of usher iterations per try */
    unsigned int _iterMax;
    /**  maximal number of usher restarts */
    unsigned int _restartMax;
    /** sets the tolerance, between the reached energy level and the goal energy level */
    double _tolerance;
    /** offset of the md domain (difference between the left, lower, front corner to 0,0,0) */
    double _offsetFromOuterBoundary;

    /** @brief a simple destructor */
    UsherParams(double rSigmaCoeff, double meanPotentialEnergyFactor, double uOverlapCoeff, double stepRefCoeff, unsigned int iterMax, unsigned int restartMax,
                double tolerance, double offsetFromOuterBoundary)
        : _rSigmaCoeff(rSigmaCoeff), _meanPotentialEnergyFactor(meanPotentialEnergyFactor), _uOverlapCoeff(uOverlapCoeff), _stepRefCoeff(stepRefCoeff),
          _iterMax(iterMax), _restartMax(restartMax), _tolerance(tolerance), _offsetFromOuterBoundary(offsetFromOuterBoundary) {}
    /** @brief a dummy destructor*/
    ~UsherParams() {}

    /** @brief returns the maximal displacement parameter for the usher method, ref. figure 7 in the paper
     *  @returns the maximal displacement*/
    double getStepRef(double numberDensity, double sigma) const {
      // if the coeff. was not parsed and thus set to -1.0, we choose the strategy
      // described in the usher paper...
      if (_stepRefCoeff == -1.0) {
        return 0.1 * pow((1.0 / numberDensity), 3.0 / 2.0) * sigma;
        // ... otherwise we fix the step size in units of sigma
      } else {
        return _stepRefCoeff * sigma;
      }
    }
  };

  /** @brief interface to the md solver */
  coupling::interface::MDSolverInterface<LinkedCell, dim>* const _mdSolverInterface;
  /** @brief instance of the UsherParams, stores all necessary parameters*/
  const UsherParams _usherParams;
};
#include "coupling/UsherParticleInsertion.cpph"
#endif // _MOLECULARDYNAMICS_COUPLING_USHERPARTICLEINSERTION_H_
