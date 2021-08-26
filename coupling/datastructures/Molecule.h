// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_DATASTRUCTURES_MOLECULE_H_
#define _MOLECULARDYNAMICS_COUPLING_DATASTRUCTURES_MOLECULE_H_

#include "tarch/la/Vector.h"
#include "coupling/interface/Molecule.h"

namespace coupling {
  namespace datastructures {
    template<unsigned int dim>
    class Molecule;
  }
}


/** molecule representation for coupling component.
 *
 *  @author Philipp Neumann
 */
template<unsigned int dim>
class coupling::datastructures::Molecule: public coupling::interface::Molecule<dim> {
  public:
    Molecule(
      const tarch::la::Vector<dim,double>& position,
      const tarch::la::Vector<dim,double>& velocity,
      const tarch::la::Vector<dim,double>& force,
      const double & potentialEnergy
    ): coupling::interface::Molecule<dim>(),
       _position(position),_velocity(velocity),_force(force),_potentialEnergy(potentialEnergy){}
    Molecule(): coupling::interface::Molecule<dim>(),
      _position(tarch::la::Vector<dim,double>(0.0)),_velocity(tarch::la::Vector<dim,double>(0.0)),
      _force(tarch::la::Vector<dim,double>(0.0)),_potentialEnergy(0.0){}
    virtual ~Molecule(){}

    /** returns/ sets the velocity of the molecule */
    tarch::la::Vector<dim,double> getVelocity() const{ return _velocity;}
    void setVelocity(const tarch::la::Vector<dim,double>& velocity){_velocity = velocity;}

    /** returns/ sets the position of the molecule */
    tarch::la::Vector<dim,double> getPosition() const {return _position;}
    void setPosition(const tarch::la::Vector<dim,double>& position){_position=position;}

    /** sets the force acting on this molecule. This function is only called in the USHER
     *  scheme so far If you want to set the force of a newly created molecule,
     *  you need to implement this function.
     */
    void setForce(const tarch::la::Vector<dim,double>& force){_force = force;}
    tarch::la::Vector<dim,double> getForce() const{ return _force;}

    /** returns/ sets the potential energy of the molecule */
    double getPotentialEnergy() const {return _potentialEnergy;}
    void setPotentialEnergy(const double& potentialEnergy){_potentialEnergy = potentialEnergy;}

  private:
    tarch::la::Vector<dim,double> _position;
    tarch::la::Vector<dim,double> _velocity;
    tarch::la::Vector<dim,double> _force;
    double _potentialEnergy;
};

#endif // _MOLECULARDYNAMICS_COUPLING_DATASTRUCTURES_MOLECULE_H_
