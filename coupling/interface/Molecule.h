// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_INTERFACE_MOLECULE_H_
#define _MOLECULARDYNAMICS_COUPLING_INTERFACE_MOLECULE_H_

namespace coupling {
  namespace interface {
    template<unsigned int dim>
    class Molecule;
  }
}

/** This class provides interface to access a single molecule in the MD simulation. The MaMiCo molecule in datastructures also inherits from this class.
 *	@brief interface to access a single molecule in the MD simulation.
 *	@tparam dim Number of dimensions; it can be 1, 2 or 3
 *  @author Philipp Neumann
 */
template<unsigned int dim>
class coupling::interface::Molecule {
  public:
  /** Destructor */
    virtual ~Molecule(){}
	
    /** returnsthe velocity of the molecule */
    virtual tarch::la::Vector<dim,double> getVelocity() const = 0;
	
	/** sets the velocity of the molecule
	 *	@param velocity
	 */
    virtual void setVelocity(const tarch::la::Vector<dim,double>& velocity) = 0;

    /** returns the position of the molecule */
    virtual tarch::la::Vector<dim,double> getPosition() const = 0;
	
	/** sets the position of the molecule
	 *	@param position
	 */
    virtual void setPosition(const tarch::la::Vector<dim,double>& position) = 0;

    /** sets the force acting on this molecule. This function is only called in the USHER
     *  scheme so far If you want to set the force of a newly created molecule,
     *  you need to implement this function.
	 *	@param force
     */
    virtual void setForce(const tarch::la::Vector<dim,double>& force) = 0;
	
	/** returns the force acting on this molecule.
	 *	@sa setForce
	 */
    virtual tarch::la::Vector<dim,double> getForce() const = 0;

    /** returns the potential energy of the molecule */
    virtual double getPotentialEnergy() const = 0;
	
	/** sets the potential energy of the molecule
	 *	@param potentialEnergy
	 */
    virtual void setPotentialEnergy(const double& potentialEnergy) = 0;
};
#endif // _MOLECULARDYNAMICS_COUPLING_INTERFACE_MOLECULE_H_
