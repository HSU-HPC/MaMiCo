// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_MOLECULE_H_
#define _MOLECULARDYNAMICS_MOLECULE_H_

#include "simplemd/MolecularDynamicsDefinitions.h"
#include "simplemd/MolecularProperties.h"
#include "tarch/la/Vector.h"

namespace simplemd {
class Molecule;
}

/** describes a single molecule
 *  @author Philipp Neumann
 */
class simplemd::Molecule {
public:
  /** initialise position and velocity of molecule */
  Molecule(const tarch::la::Vector<MD_DIM, double>& position, const tarch::la::Vector<MD_DIM, double>& velocity)
      : _position(position), _velocity(velocity),
#if (AD_RES == MD_YES)
        _twoBodyForce(0.0), _threeBodyForce(0.0),
#endif      
        _force(0.0), _forceOld(0.0),
#if (AD_RES == MD_YES)
        _twoBodyPotentialEnergy(0.0), _threeBodyPotentialEnergy(0.0),
#endif
        _potentialEnergy(0.0), _id(0), _isFixed(false) {}
  /** empty constructor */
  Molecule() : _position(0.0), _velocity(0.0),
#if (AD_RES == MD_YES)
               _twoBodyForce(0.0), _threeBodyForce(0.0),
#endif
               _force(0.0), _forceOld(0.0),
#if (AD_RES == MD_YES)
               _twoBodyPotentialEnergy(0.0), _threeBodyPotentialEnergy(0.0),
#endif
               _potentialEnergy(0.0), _id(0), _isFixed(false) {}
  ~Molecule() {}

  /** fix position of particle. */
  void fix() {
    _isFixed = true;
    _velocity = tarch::la::Vector<MD_DIM, double>(0.0);
  }

  /** check if particle is fixed */
  const bool& isFixed() const { return _isFixed; }

  /** get/set ID */
  const unsigned int& getID() const { return _id; }
  void setID(const unsigned int& id) { _id = id; }

  /** get/ set position */
  tarch::la::Vector<MD_DIM, double>& getPosition() { return _position; }
  const tarch::la::Vector<MD_DIM, double>& getConstPosition() const { return _position; }
  void setPosition(const tarch::la::Vector<MD_DIM, double>& position) { _position = position; }

  /** get/ set velocity */
  tarch::la::Vector<MD_DIM, double>& getVelocity() { return _velocity; }
  const tarch::la::Vector<MD_DIM, double>& getConstVelocity() const { return _velocity; }
  void setVelocity(const tarch::la::Vector<MD_DIM, double>& velocity) { _velocity = velocity; }

#if (AD_RES == MD_YES)
  /** get/ set two-body-force */
  tarch::la::Vector<MD_DIM, double>& getTwoBodyForce() { return _twoBodyForce; }
  const tarch::la::Vector<MD_DIM, double>& getConstTwoBodyForce() const { return _twoBodyForce; }
  void setTwoBodyForce(const tarch::la::Vector<MD_DIM, double>& force) { _twoBodyForce = force; }
  
  /** get/ set three-body-force */
  tarch::la::Vector<MD_DIM, double>& getThreeBodyForce() { return _threeBodyForce; }
  const tarch::la::Vector<MD_DIM, double>& getConstThreeBodyForce() const { return _threeBodyForce; }
  void setThreeBodyForce(const tarch::la::Vector<MD_DIM, double>& force) { _threeBodyForce = force; }
#endif

  /** get/ set force */
  tarch::la::Vector<MD_DIM, double>& getForce() { return _force; }
  const tarch::la::Vector<MD_DIM, double>& getConstForce() const { return _force; }
  void setForce(const tarch::la::Vector<MD_DIM, double>& force) { _force = force; }

  /** get/ set force of last timestep */
  tarch::la::Vector<MD_DIM, double>& getForceOld() { return _forceOld; }
  const tarch::la::Vector<MD_DIM, double>& getConstForceOld() const { return _forceOld; }
  void setForceOld(const tarch::la::Vector<MD_DIM, double>& force) { _forceOld = force; }

#if (AD_RES == MD_YES)
  /** get/ set two-body-potential */
  double& getTwoBodyPotentialEnergy() { return _twoBodyPotentialEnergy; }
  const double& getConstTwoBodyPotentialEnergy() const { return _twoBodyPotentialEnergy; }
  void setTwoBodyPotentialEnergy(const double& potentialEnergy) { _twoBodyPotentialEnergy = potentialEnergy; }
  
  /** get/ set three-body-potential */
  double& getThreeBodyPotentialEnergy() { return _threeBodyPotentialEnergy; }
  const double& getConstThreeBodyPotentialEnergy() const { return _threeBodyPotentialEnergy; }
  void setThreeBodyPotentialEnergy(const double& potentialEnergy) { _threeBodyPotentialEnergy = potentialEnergy; }
#endif

  double& getPotentialEnergy() { return _potentialEnergy; }
  const double& getConstPotentialEnergy() const { return _potentialEnergy; }
  void setPotentialEnergy(const double& potentialEnergy) { _potentialEnergy = potentialEnergy; }

private:
  tarch::la::Vector<MD_DIM, double> _position;
  tarch::la::Vector<MD_DIM, double> _velocity;
#if (AD_RES == MD_YES)
  tarch::la::Vector<MD_DIM, double> _twoBodyForce;
  tarch::la::Vector<MD_DIM, double> _threeBodyForce;
#endif
  tarch::la::Vector<MD_DIM, double> _force;
  tarch::la::Vector<MD_DIM, double> _forceOld;
#if (AD_RES == MD_YES)
  double _twoBodyPotentialEnergy;
  double _threeBodyPotentialEnergy;
#endif
  double _potentialEnergy;
  unsigned int _id;
  bool _isFixed;
};

#endif
