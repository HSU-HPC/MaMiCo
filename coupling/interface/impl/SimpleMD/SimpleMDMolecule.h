// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_INTERFACE_SIMPLEMDMOLECULE_H_
#define _MOLECULARDYNAMICS_COUPLING_INTERFACE_SIMPLEMDMOLECULE_H_

#include "coupling/interface/Molecule.h"
#include "simplemd/Molecule.h"
#include "simplemd/MolecularDynamicsDefinitions.h"

namespace coupling {
  namespace interface{
    class SimpleMDMolecule;
  }
}


/** interface to access molecule of SimpleMD.
 *  @author Philipp Neumann
 */
class coupling::interface::SimpleMDMolecule: public coupling::interface::Molecule<MD_DIM> {
  public:
    SimpleMDMolecule(simplemd::Molecule *myMolecule): _myMolecule(myMolecule){}
    virtual ~SimpleMDMolecule(){}

    void setMolecule(simplemd::Molecule *newMolecule){ _myMolecule=newMolecule;}

    /** returns/ sets the velocity of the molecule */
    virtual tarch::la::Vector<MD_DIM,double> getVelocity() const {
      return _myMolecule->getConstVelocity();
    }
    virtual void setVelocity(const tarch::la::Vector<MD_DIM,double>& velocity){
      _myMolecule->setVelocity(velocity);
    }

    /** returns/ sets the position of the molecule */
    virtual tarch::la::Vector<MD_DIM,double> getPosition() const {
      return _myMolecule->getConstPosition();
    }

    virtual void setPosition(const tarch::la::Vector<MD_DIM,double>& position){
      _myMolecule->setPosition(position);
    }

    /** sets the force acting on this molecule. This function is only called in the USHER
     *  scheme so far If you want to set the force of a newly created molecule,
     *  you need to implement this function.
     */
    virtual void setForce(const tarch::la::Vector<MD_DIM,double>& force){
      _myMolecule->setForce(force);
    }
    virtual tarch::la::Vector<MD_DIM,double> getForce() const{
      return _myMolecule->getConstForce();
    }

    /** returns/ sets the potential energy of the molecule */
    virtual double getPotentialEnergy() const{
      return _myMolecule->getPotentialEnergy();
    }
    virtual void setPotentialEnergy(const double& potentialEnergy){
      _myMolecule->setPotentialEnergy(potentialEnergy);
    }

  private:
    simplemd::Molecule* _myMolecule;
};
#endif // _MOLECULARDYNAMICS_COUPLING_INTERFACE_SIMPLEMDMOLECULE_H_
