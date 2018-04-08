// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_USHERPARTICLEINSERTION_H_
#define _MOLECULARDYNAMICS_COUPLING_USHERPARTICLEINSERTION_H_

#include "coupling/cell-mappings/ComputeMassMapping.h"
#include "coupling/cell-mappings/ComputeMomentumMapping.h"
#include "coupling/cell-mappings/ComputeTemperatureMapping.h"
#include "coupling/cell-mappings/DeleteParticleMapping.h"
#include "coupling/datastructures/MacroscopicCell.h"
#include "coupling/datastructures/Molecule.h"
#include "coupling/interface/MDSolverInterface.h"
#include "coupling/CouplingMDDefinitions.h"
#include "coupling/ParticleInsertion.h"
#include "tarch/utils/RandomNumberService.h"


namespace coupling {
  template<class LinkedCell,unsigned int dim>
  class UsherParticleInsertion;
}


/** handles particle insertion (via Usher algorithm) and random particle deletion.
 *  We currently also allow particle insertion in energyy holes, i.e. for energy(particle)==0.0.
 *
 *  @author Philipp Neumann
 */
template<class LinkedCell,unsigned int dim>
class coupling::UsherParticleInsertion: public coupling::ParticleInsertion<LinkedCell,dim> {
  public:

    UsherParticleInsertion(
     unsigned int insertDeleteMassEveryTimestep,
     double rSigmaCoeff, double meanPotentialEnergyFactor,double uOverlapCoeff,
     double stepRefCoeff,unsigned int iterMax,unsigned int restartMax,double tolerance, double offsetFromOuterBoundary,
      coupling::interface::MDSolverInterface<LinkedCell,dim> * const mdSolverInterface
    );
    virtual ~UsherParticleInsertion(){}

    /** this state is returned by the insertDeleteMass() function and tells the user, if mass was inserted/deleted
     *  or if nothing happened at all.
     */
    enum Action{NoAction=0,Insertion=1,Deletion=2};

    virtual typename coupling::ParticleInsertion<LinkedCell,dim>::Action insertDeleteMass(
      coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell,dim>& cell,
      const tarch::la::Vector<dim,double>& macroscopicCellPosition,
      const tarch::la::Vector<dim,double>& macroscopicCellSize,
      const tarch::la::Vector<dim,double>& meanVelocity,
      const double &temperature
    ) const;

    virtual bool requiresPotentialEnergyLandscape(){ return true; }

  private:
    /** inserts a particle in the respective macroscopic cell. Returns Insertion on success and NoAction otherwise. */
    typename coupling::ParticleInsertion<LinkedCell,dim>::Action insertParticle(
      coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell,dim>& cell,
      const tarch::la::Vector<dim,double>& macroscopicCellPosition,
      const tarch::la::Vector<dim,double>& macroscopicCellSize,
      const tarch::la::Vector<dim,double>& meanVelocity,
      const double &temperature
    ) const;

    /** deletes a particle from the macroscopic cell. Returns Deletion on success and NoAction otherwise. */
    typename coupling::ParticleInsertion<LinkedCell,dim>::Action deleteParticle(
      coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell,dim>& cell
    ) const;


    /** determines the position of a new particle within the macroscopic cell thisCell and stores the result
     *  in the position entry of the molecule 'molecule'. The method returns Insertion if a position has been
     *  found and NoAction otherwise.
     */
  protected:
    virtual typename coupling::ParticleInsertion<LinkedCell,dim>::Action findParticlePosition(
      coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell,dim>& thisCell,
      const tarch::la::Vector<dim,double>& macroscopicCellPosition,
      const tarch::la::Vector<dim,double>& macroscopicCellSize,
      coupling::datastructures::Molecule<dim>& molecule
    ) const;


    class UsherParams {
      public:
        double _rSigmaCoeff;
        double _meanPotentialEnergyFactor;
        double _uOverlapCoeff;
        double _stepRefCoeff;
        unsigned int _iterMax;
        unsigned int _restartMax;
        double _tolerance;
        double _offsetFromOuterBoundary;

        UsherParams(double rSigmaCoeff, double meanPotentialEnergyFactor,double uOverlapCoeff,
            double stepRefCoeff,unsigned int iterMax,unsigned int restartMax,double tolerance,double offsetFromOuterBoundary):
        _rSigmaCoeff(rSigmaCoeff), _meanPotentialEnergyFactor(meanPotentialEnergyFactor), _uOverlapCoeff(uOverlapCoeff),
        _stepRefCoeff(stepRefCoeff),_iterMax(iterMax),_restartMax(restartMax),_tolerance(tolerance),_offsetFromOuterBoundary(offsetFromOuterBoundary){}
        ~UsherParams(){}


        double getStepRef(double numberDensity, double sigma) const {
          // if the coeff. was not parsed and thus set to -1.0, we choose the strategy
          // described in the usher paper...
          if (_stepRefCoeff==-1.0){
            return 0.1*pow((1.0/numberDensity),3.0/2.0)*sigma;
          // ... otherwise we fix the step size in units of sigma
          } else {
            return _stepRefCoeff*sigma;
          }
        }
    };

    coupling::interface::MDSolverInterface<LinkedCell,dim> * const _mdSolverInterface;
    const UsherParams _usherParams;
};
#include "coupling/UsherParticleInsertion.cpph"
#endif // _MOLECULARDYNAMICS_COUPLING_USHERPARTICLEINSERTION_H_
