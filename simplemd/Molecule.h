// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_MOLECULE_H_
#define _MOLECULARDYNAMICS_MOLECULE_H_

#include "tarch/la/Vector.h"
#include "simplemd/MolecularDynamicsDefinitions.h"
#include "simplemd/MolecularProperties.h"

namespace simplemd {
  class Molecule;
}


/** describes a single molecule
 *  @author Philipp Neumann
 */
class simplemd::Molecule {
  public:
    /** initialise position and velocity of molecule */
    Molecule(
      const tarch::la::Vector<MD_DIM,double> &position,
      const tarch::la::Vector<MD_DIM,double> &velocity
    ): _position(position), _velocity(velocity), _force(0.0), _forceOld(0.0), _potentialEnergy(0.0),_id(0), _isFixed(false){}
    /** empty constructor */
    Molecule(): _position(0.0), _velocity(0.0), _force(0.0), _forceOld(0.0), _potentialEnergy(0.0), _id(0), _isFixed(false){}
    ~Molecule(){}

    /** fix position of particle. */
    void fix() {_isFixed = true; _velocity = tarch::la::Vector<MD_DIM, double>(0.0); }

    /** check if particle is fixed */
    const bool& isFixed() const {return _isFixed;}

    /** get/set ID */
    const unsigned int& getID() const { return _id;}
    void setID(const unsigned int& id){ _id = id; }

    /** get/ set position */
    tarch::la::Vector<MD_DIM,double>& getPosition(){ return _position;}
    const tarch::la::Vector<MD_DIM,double>& getConstPosition() const { return _position;}
    void setPosition(const tarch::la::Vector<MD_DIM,double>& position){ _position = position; }

    /** get/ set velocity */
    tarch::la::Vector<MD_DIM,double>& getVelocity(){ return _velocity;}
    const tarch::la::Vector<MD_DIM,double>& getConstVelocity() const { return _velocity;}
    void setVelocity(const tarch::la::Vector<MD_DIM,double>& velocity){ _velocity = velocity; }

    /** get/ set force */
    tarch::la::Vector<MD_DIM,double>& getForce(){ return _force;}
    const tarch::la::Vector<MD_DIM,double>& getConstForce() const { return _force;}
    void setForce(const tarch::la::Vector<MD_DIM,double>& force){ _force = force;}

    /** get/ set force of last timestep */
    tarch::la::Vector<MD_DIM,double>& getForceOld(){ return _forceOld;}
    const tarch::la::Vector<MD_DIM,double>& getConstForceOld() const { return _forceOld;}
    void setForceOld(const tarch::la::Vector<MD_DIM,double>& force){ _forceOld = force;}

    double& getPotentialEnergy(){ return _potentialEnergy; }
    const double& getConstPotentialEnergy() const { return _potentialEnergy; }
    void setPotentialEnergy(const double& potentialEnergy){ _potentialEnergy = potentialEnergy; }

  private:
    tarch::la::Vector<MD_DIM,double> _position;
    tarch::la::Vector<MD_DIM,double> _velocity;
    tarch::la::Vector<MD_DIM,double> _force;
    tarch::la::Vector<MD_DIM,double> _forceOld;
    double _potentialEnergy;
    unsigned int _id;
    bool _isFixed;
};

#endif

