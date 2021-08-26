// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_MDSOLVERINTERFACE_H_
#define _MOLECULARDYNAMICS_COUPLING_MDSOLVERINTERFACE_H_

#include "tarch/la/Vector.h"
#include "coupling/IndexConversion.h"
#include "coupling/interface/Molecule.h"
#include "coupling/interface/MoleculeIterator.h"
#include <list>

namespace coupling {
  namespace interface {
    template<class LinkedCell,unsigned int dim>
    class MDSolverInterface;
  }
}


/** interface to the moleculardynamics simulation
 *
 *  @author Philipp Neumann
 */
template<class LinkedCell,unsigned int dim>
class coupling::interface::MDSolverInterface {
  public:
    virtual ~MDSolverInterface(){}

    /** returns a particular linked cell inside a macroscopic cell.
     *  The macroscopic cells are currently located on the same process as the respective linked cells.
     *  However, several linked cells may be part of a macroscopic cell.
     *  The macroscopic cells also contain a ghost layer which surrounds each local domain; the very
     *  first macroscopic cell inside the global MD domain (or local MD domain) is thus given by coordinates
     *  (1,1,1) (or (1,1) in 2D, respectively).
     *  The index linkedCellInMacroscopicCell corresponds to the coordinates of the linked cell inside the
     *  given macroscopic cell. These coordinates thus lie in a range (0,linkedCellsPerMacroscopicCell-1).
     */
    virtual LinkedCell& getLinkedCell(
      const tarch::la::Vector<dim,unsigned int>& macroscopicCellIndex,
      const tarch::la::Vector<dim,unsigned int>& linkedCellInMacroscopicCell,
      const tarch::la::Vector<dim,unsigned int>& linkedCellsPerMacroscopicCell,
      const coupling::IndexConversion<dim> &indexConversion
    ) = 0;

    /** returns the global size of the box-shaped MD domain */
    virtual tarch::la::Vector<dim,double> getGlobalMDDomainSize() const = 0;

    /** returns the offset (i.e. lower,left corner) of MD domain */
    virtual tarch::la::Vector<dim,double> getGlobalMDDomainOffset() const = 0;

    /** returns the mass of a single fluid molecule */
    virtual double getMoleculeMass() const  = 0;

    /** returns Boltzmann's constant */
    virtual double getKB() const = 0;

    /** returns the sigma parameter of the LJ potential */
    virtual double getMoleculeSigma() const = 0;

    /** returns the epsilon parameter of the LJ potential */
    virtual double getMoleculeEpsilon() const = 0;

    /** sets a random velocity in the vector 'initialVelocity'. This velocity is sampled from
     *  a Maxwellian assuming a mean flow velocity 'meanVelocity' and a temperature 'temperature'
     *  of the fluid.
     */
    virtual void getInitialVelocity(
      const tarch::la::Vector<dim,double>& meanVelocity, const double &kB, const double& temperature, tarch::la::Vector<dim,double>& initialVelocity
    ) const = 0;

    /** deletes the molecule from the MD simulation */
    virtual void deleteMoleculeFromMDSimulation(const coupling::interface::Molecule<dim>& molecule,LinkedCell& cell) = 0;

    /** adds the molecule to the MD simulation.
     */
    virtual void addMoleculeToMDSimulation(const coupling::interface::Molecule<dim>& molecule) = 0;

    /** sets up the potential energy landscape over the domain spanned by indexOfFirstMacroscopicCell and
     *  rangeCoordinates. The first vector denotes the position of the lower,left,front corner of the
     *  domain, rangeCoordinates the number of macroscopic cells in each spatial direction in which the
     *  potential energy needs to be computed.
     */
    virtual void setupPotentialEnergyLandscape(
      const tarch::la::Vector<dim,unsigned int>& indexOfFirstMacroscopicCell,
      const tarch::la::Vector<dim,unsigned int>& rangeMacroscopicCells,
      const tarch::la::Vector<dim,unsigned int>& linkedCellsPerMacroscopicCell
    ) = 0;

    /** returns the local index vector (w.r.t. lexicographic ordering of the linked cells in the MD simulation,
     *  see getLinkedCell()) for the linked cell that the position 'position' belongs to. This method is
     *  crucial for the USHER particle insertion scheme as we need to loop over all neighboured linked cells
     *  and the cell itself to determine potential energy and forces acting on the molecule at position 'position'.
     */
    virtual tarch::la::Vector<dim,unsigned int> getLinkedCellIndexForMoleculePosition(
      const tarch::la::Vector<dim,double>& position
    ) = 0;


    /** assumes that a molecule is placed somewhere inside the linked cell at index
     *  'linkedCellIndex' and computes the force and potential energy contributions from all molecules
     *  in the same linked cell and the neighboured linked cells onto this molecule.
     *  The molecule "molecule" is considered to NOT be part of the simulation domain and thus the linked cells.
     *  Therefore, another molecule inside the linked cells may even coincide with "molecule".
     *  The results are stored within the molecule.
     */
    virtual void calculateForceAndEnergy(
      coupling::interface::Molecule<dim> &molecule
    ) = 0;

    /** is called each time when MaMiCo tried to insert/ delete molecules from the MD simulation. As a consequence,
     *  a synchronization between processes or with local boundary data might be necessary.
     *  Example: If in the lower left cell of a 2D MD simulation a molecule was inserted and we use periodic
     *  boundary conditions, we need to provide the inserted molecule (amongst others) to the upper right cell
     *  by adding this molecule to the ghost boundary layer (when using the builtin MD simulation).
     *  For the builtin MD simulation, the implementation of this method clears all molecules from the ghost
     *  layers and re-fills the ghost layers again.
     */
    virtual void synchronizeMoleculesAfterMassModification() = 0;

    /** is called each time when MaMiCo tried to insert momentum in the MD simulation. For the builtin MD simulation,
     *  this method is empty as the simulation does not need to synchronize changing momentum over the processes
     *  within a timestep (only the positions and numbers of molecules in possible ghost layers matter!).
     *  However, other solvers might need to implement something in here.
     */
    virtual void synchronizeMoleculesAfterMomentumModification() = 0;

    /** returns the timestep of the MD simulation */
    virtual double getDt() = 0;

    /** returns a new molecule iterator for a certain linked cell */
    virtual coupling::interface::MoleculeIterator<LinkedCell,dim>* getMoleculeIterator(LinkedCell& cell) = 0;
};

#endif  // _MOLECULARDYNAMICS_COUPLING_MDSOLVERINTERFACE_H_

