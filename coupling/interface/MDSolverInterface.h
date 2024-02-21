// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_MDSOLVERINTERFACE_H_
#define _MOLECULARDYNAMICS_COUPLING_MDSOLVERINTERFACE_H_

#include "coupling/indexing/IndexTypes.h"
#include "coupling/interface/Molecule.h"
#include "coupling/interface/MoleculeIterator.h"
#include "tarch/la/Vector.h"
#include <list>

namespace coupling {
namespace interface {
template <class LinkedCell, unsigned int dim> class MDSolverInterface;
}
} // namespace coupling

/** This class provides
 *	@brief interface to the MD simulation
 *	@tparam dim Number of dimensions; it can be 1, 2 or 3
 *  @author Philipp Neumann
 */
template <class LinkedCell, unsigned int dim> class coupling::interface::MDSolverInterface {
public:
  using CellIndex_T = I11;

  /** Destructor */
  virtual ~MDSolverInterface() {}

  /** This function specifies a particular linked cell inside a coupling
   *cell. The coupling cells are currently located on the same process as the
   *respective linked cells. However, several linked cells may be part of a
   *coupling cell. The coupling cells also contain a ghost layer which
   *surrounds each local domain; the very first coupling cell inside the
   *global MD domain (or local MD domain) is thus given by coordinates (1,1,1)
   *(or (1,1) in 2D, respectively). The index linkedCellInCouplingCell
   *corresponds to the coordinates of the linked cell inside the given
   *coupling cell. These coordinates thus lie in a range
   *(0,linkedCellsPerCouplingCell-1).
   *	@param couplingCellIndex
   *	@param linkedCellInCouplingCell
   *	@param linkedCellsPerCouplingCell
   *	@param indexConversion
   *	@returns a particular linked cell inside a coupling cell.
   */
  virtual LinkedCell& getLinkedCell(const CellIndex_T& couplingCellIndex, const tarch::la::Vector<dim, unsigned int>& linkedCellInCouplingCell,
                                    const tarch::la::Vector<dim, unsigned int>& linkedCellsPerCouplingCell) = 0;

  /** This function specifies the global size of the box-shaped MD domain
   *  @returns the global size of the box-shaped MD domain */
  virtual tarch::la::Vector<dim, double> getGlobalMDDomainSize() const = 0;

  /** This function determines the offset (i.e. lower,left corner) of MD domain
   *  @returns the offset (i.e. lower,left corner) of MD domain */
  virtual tarch::la::Vector<dim, double> getGlobalMDDomainOffset() const = 0;

  /** This function specifies the mass of a single fluid molecule
   *  @returns the mass of a single fluid molecule */
  virtual double getMoleculeMass() const = 0;

  /**
   *  @returns Boltzmann's constant */
  virtual double getKB() const = 0;

  /**
   *  @returns the sigma parameter of the Lennard-Jones potential
   */
  virtual double getMoleculeSigma() const = 0;

  /**
   *  @returns the epsilon parameter of the Lennard-Jones potential */
  virtual double getMoleculeEpsilon() const = 0;

  /** This function sets a random velocity in the vector 'initialVelocity'. This
   *velocity is sampled from a Maxwellian assuming a mean flow velocity
   *'meanVelocity' and a temperature 'temperature' of the fluid.
   *	@param meanVelocity
   *	@param kB
   *	@param temperature
   *	@param initialVelocity
   */
  virtual void getInitialVelocity(const tarch::la::Vector<dim, double>& meanVelocity, const double& kB, const double& temperature,
                                  tarch::la::Vector<dim, double>& initialVelocity) const = 0;

  /** This function deletes the molecule from the MD simulation
   *	@param molecule
   *	@param cell
   */
  virtual void deleteMoleculeFromMDSimulation(const coupling::interface::Molecule<dim>& molecule, LinkedCell& cell) = 0;

  /** This function adds the molecule to the MD simulation.
   *	@param molecule
   */
  virtual void addMoleculeToMDSimulation(const coupling::interface::Molecule<dim>& molecule) = 0;

  /** This function sets up the potential energy landscape over the domain
   *spanned by indexOfFirstCouplingCell and rangeCoordinates. The first
   *vector denotes the position of the lower,left,front corner of the domain,
   *rangeCoordinates the number of coupling cells in each spatial direction
   *in which the potential energy needs to be computed.
   *	@param indexOfFirstCouplingCell
   *	@param rangeCouplingCells
   *	@param linkedCellsPerCouplingCell
   */
  virtual void setupPotentialEnergyLandscape(const tarch::la::Vector<dim, unsigned int>& indexOfFirstCouplingCell,
                                             const tarch::la::Vector<dim, unsigned int>& rangeCouplingCells,
                                             const tarch::la::Vector<dim, unsigned int>& linkedCellsPerCouplingCell) = 0;

  /** This function specifies the local index vector (w.r.t. lexicographic
   *ordering of the linked cells in the MD simulation, for the linked cell that
   *the position 'position' belongs to. This method is crucial for the USHER
   *particle insertion scheme as we need to loop over all neighboured linked
   *cells and the cell itself to determine potential energy and forces acting on
   *the molecule at position 'position'.
   *	@param position
   *	@return the local index vector for the linked cell that the position
   *'position' belongs to
   *	@sa getLinkedCell()
   */
  virtual tarch::la::Vector<dim, unsigned int> getLinkedCellIndexForMoleculePosition(const tarch::la::Vector<dim, double>& position) = 0;

  /** This function assumes that a molecule is placed somewhere inside the
   *linked cell at index 'linkedCellIndex' and computes the force and potential
   *energy contributions from all molecules in the same linked cell and the
   *neighboured linked cells onto this molecule. The molecule "molecule" is
   *considered to NOT be part of the simulation domain and thus the linked
   *cells. Therefore, another molecule inside the linked cells may even coincide
   *with "molecule". The results are stored within the molecule.
   *	@param molecule
   */
  virtual void calculateForceAndEnergy(coupling::interface::Molecule<dim>& molecule) = 0;

  /** This function is called each time when MaMiCo tried to insert/ delete
   * molecules from the MD simulation. As a consequence, a synchronization
   * between processes or with local boundary data might be necessary. Example:
   * If in the lower left cell of a 2D MD simulation a molecule was inserted and
   * we use periodic boundary conditions, we need to provide the inserted
   * molecule (amongst others) to the upper right cell by adding this molecule
   * to the ghost boundary layer (when using the builtin MD simulation). For the
   * builtin MD simulation, the implementation of this method clears all
   * molecules from the ghost layers and re-fills the ghost layers again.
   */
  virtual void synchronizeMoleculesAfterMassModification() = 0;

  /** This function is called each time when MaMiCo tried to insert momentum in
   * the MD simulation. For the builtin MD simulation, this method is empty as
   * the simulation does not need to synchronize changing momentum over the
   * processes within a timestep (only the positions and numbers of molecules in
   * possible ghost layers matter!). However, other solvers might need to
   * implement something in here.
   */
  virtual void synchronizeMoleculesAfterMomentumModification() = 0;

  /**
   *  @returns the timestep of the MD simulation
   */
  virtual double getDt() = 0;

  /**
   *	@param cell
   *  @returns a new molecule iterator for a certain linked cell
   */
  virtual coupling::interface::MoleculeIterator<LinkedCell, dim>* getMoleculeIterator(LinkedCell& cell) = 0;
};

#endif // _MOLECULARDYNAMICS_COUPLING_MDSOLVERINTERFACE_H_
