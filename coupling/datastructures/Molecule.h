// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_DATASTRUCTURES_MOLECULE_H_
#define _MOLECULARDYNAMICS_COUPLING_DATASTRUCTURES_MOLECULE_H_

#include "coupling/interface/Molecule.h"
#include "tarch/la/Vector.h"

namespace coupling {
namespace datastructures {
template <unsigned int dim> class Molecule;
}
} // namespace coupling

/**
 *	@brief molecule representation for coupling component. Dericed from the
 *class coupling::interface::Molecule
 *	@tparam dim Number of dimensions; it can be 1, 2 or 3
 *  @author Philipp Neumann
 */

template <unsigned int dim> class coupling::datastructures::Molecule : public coupling::interface::Molecule<dim> {
public:
  /** Constructor: initialises the molecule;
   *	@param position
   * 	@param velocity
   * 	@param force
   * 	@param potentialEnergy
   */
  Molecule(const tarch::la::Vector<dim, double>& position, const tarch::la::Vector<dim, double>& velocity, const tarch::la::Vector<dim, double>& force,
           const double& potentialEnergy)
      : coupling::interface::Molecule<dim>(), _position(position), _velocity(velocity), _force(force), _potentialEnergy(potentialEnergy) {}
  Molecule()
      : coupling::interface::Molecule<dim>(), _position(tarch::la::Vector<dim, double>(0.0)), _velocity(tarch::la::Vector<dim, double>(0.0)),
        _force(tarch::la::Vector<dim, double>(0.0)), _potentialEnergy(0.0) {}

  /** Destructor */
  virtual ~Molecule() {}

  /** returns the velocity of the molecule
   * @return _velocity Velocity*/
  tarch::la::Vector<dim, double> getVelocity() const { return _velocity; }
  /** sets the velocity of the molecule
   * @param velocity Velocity*/
  void setVelocity(const tarch::la::Vector<dim, double>& velocity) { _velocity = velocity; }

  /** returns the velocity of the molecule
   * @return _position Position*/
  tarch::la::Vector<dim, double> getPosition() const { return _position; }
  /** sets the velocity of the molecule
   * @param position Position*/
  void setPosition(const tarch::la::Vector<dim, double>& position) { _position = position; }

  /** sets the force acting on this molecule. This function is called in the
   * USHER scheme so far only if the force of a newly created molecule should be
   * set
   *  @param force Force
   *  @todo Philipp  When the force should be set? when not? you need to
   * implement this function.
   */
  void setForce(const tarch::la::Vector<dim, double>& force) { _force = force; }
  /** returns the force of the molecule
   *  @return _force Force*/
  tarch::la::Vector<dim, double> getForce() const { return _force; }

  /** returns potential energy of the molecule
   * @return _potentialEnergy Potential energy of the molecule */
  double getPotentialEnergy() const { return _potentialEnergy; }
  /** sets potential energy of the molecule
   * @param _potentialEnergy Potential energy of the molecule */
  void setPotentialEnergy(const double& potentialEnergy) { _potentialEnergy = potentialEnergy; }

private:
  /** Position of the molecule */
  tarch::la::Vector<dim, double> _position;
  /** Velocity vector of the molecule */
  tarch::la::Vector<dim, double> _velocity;
  /** Force vector of the molecule */
  tarch::la::Vector<dim, double> _force;
  /** Potential energy of the molecule */
  double _potentialEnergy;
};

#endif // _MOLECULARDYNAMICS_COUPLING_DATASTRUCTURES_MOLECULE_H_
