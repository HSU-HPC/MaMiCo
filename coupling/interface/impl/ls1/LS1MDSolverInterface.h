#ifndef LS1_MD_SOLVER_INTERFACE_H_
#define LS1_MD_SOLVER_INTERFACE_H_

#include "coupling/interface/MDSolverInterface.h"
#include "coupling/interface/impl/ls1/LS1RegionWrapper.h"
#include "coupling/interface/impl/ls1/LS1MoleculeIterator.h"

#include "ls1/src/ensemble/EnsembleBase.h"
#include "ls1/src/integrators/Integrator.h"
#include "ls1/src/parallel/DomainDecompBase.h"

namespace coupling
{
    namespace interface
    {
        class LS1MDSolverInterface;
    }
}

class coupling::interface::LS1MDSolverInterface : public coupling::interface::MDSolverInterface<ls1::LS1RegionWrapper, 3>
{
    /** returns a particular linked cell inside a macroscopic cell.
     *  The macroscopic cells are currently located on the same process as the respective linked cells.
     *  However, several linked cells may be part of a macroscopic cell.
     *  The macroscopic cells also contain a ghost layer which surrounds each local domain; the very
     *  first macroscopic cell inside the global MD domain (or local MD domain) is thus given by coordinates
     *  (1,1,1) (or (1,1) in 2D, respectively).
     *  The index linkedCellInMacroscopicCell corresponds to the coordinates of the linked cell inside the
     *  given macroscopic cell. These coordinates thus lie in a range (0,linkedCellsPerMacroscopicCell-1).
     */
    virtual ls1::LS1RegionWrapper& getLinkedCell(
      const tarch::la::Vector<3,unsigned int>& macroscopicCellIndex,
      const tarch::la::Vector<3,unsigned int>& linkedCellInMacroscopicCell,
      const tarch::la::Vector<3,unsigned int>& linkedCellsPerMacroscopicCell,
      const coupling::IndexConversion<3> &indexConversion
    )
    {
      //only one cell per macroscopic cell
      double bBoxMin[3];
			double bBoxMax[3];
      ::global_simulation->domainDecomposition().getBoundingBoxMinMax(::global_simulation->getDomain(), bBoxMin, bBoxMax);
      ls1::LS1RegionWrapper cell(bBoxMin, bBoxMax);
      return cell;
    }

    /** returns the global size of the box-shaped MD domain */
    virtual tarch::la::Vector<3,double> getGlobalMDDomainSize() const
    {
      auto up = ::global_simulation->getEnsemble()->domain()->rmax();
      auto down = ::global_simulation->getEnsemble()->domain()->rmin();
      tarch::la::Vector<3, double> globalSize = {up[0]-down[0], up[1]-down[1], up[2]-down[2]};
      return globalSize;
    }

    /** returns the offset (i.e. lower,left corner) of MD domain */
    virtual tarch::la::Vector<3,double> getGlobalMDDomainOffset() const
    {
      auto down = ::global_simulation->getEnsemble()->domain()->rmin();
      tarch::la::Vector<3, double> globalOffset(down[0], down[1], down[2]);
      return globalOffset;
    }

    /** returns the mass of a single fluid molecule */
    virtual double getMoleculeMass() const { return ::global_simulation->getEnsemble()->getComponent(0)->m(); }

    /** returns Boltzmann's constant */
    virtual double getKB() const { return 1.0; }

    /** returns the sigma parameter of the LJ potential */
    virtual double getMoleculeSigma() const { return ::global_simulation->getEnsemble()->getComponent(0)->getSigma(0);  }

    /** returns the epsilon parameter of the LJ potential */
    virtual double getMoleculeEpsilon() const { return ::global_simulation->getEnsemble()->getComponent(0)->ljcenter(0).eps(); }

    /** sets a random velocity in the vector 'initialVelocity'. This velocity is sampled from
     *  a Maxwellian assuming a mean flow velocity 'meanVelocity' and a temperature 'temperature'
     *  of the fluid.
     */
    virtual void getInitialVelocity(
      const tarch::la::Vector<3,double>& meanVelocity, const double &kB, const double& temperature, tarch::la::Vector<3,double>& initialVelocity
    ) const = 0;

    /** deletes the molecule from the MD simulation */
    virtual void deleteMoleculeFromMDSimulation(const coupling::interface::Molecule<3>& molecule,ls1::LS1RegionWrapper& cell) { cell.deleteMolecule(molecule); }

    /** adds the molecule to the MD simulation.
     */
    virtual void addMoleculeToMDSimulation(const coupling::interface::Molecule<3>& molecule) 
    { 
      double bBoxMin[3];
			double bBoxMax[3];
      ::global_simulation->domainDecomposition().getBoundingBoxMinMax(::global_simulation->getDomain(), bBoxMin, bBoxMax);
      ls1::LS1RegionWrapper cell(bBoxMin, bBoxMax);
      cell.addMolecule(molecule);
    }

    /** sets up the potential energy landscape over the domain spanned by indexOfFirstMacroscopicCell and
     *  rangeCoordinates. The first vector denotes the position of the lower,left,front corner of the
     *  domain, rangeCoordinates the number of macroscopic cells in each spatial direction in which the
     *  potential energy needs to be computed.
     */
    virtual void setupPotentialEnergyLandscape(
      const tarch::la::Vector<3,unsigned int>& indexOfFirstMacroscopicCell,
      const tarch::la::Vector<3,unsigned int>& rangeMacroscopicCells,
      const tarch::la::Vector<3,unsigned int>& linkedCellsPerMacroscopicCell
    ) = 0;

    /** returns the local index vector (w.r.t. lexicographic ordering of the linked cells in the MD simulation,
     *  see getLinkedCell()) for the linked cell that the position 'position' belongs to. This method is
     *  crucial for the USHER particle insertion scheme as we need to loop over all neighboured linked cells
     *  and the cell itself to determine potential energy and forces acting on the molecule at position 'position'.
     */
    virtual tarch::la::Vector<3,unsigned int> getLinkedCellIndexForMoleculePosition(
      const tarch::la::Vector<3,double>& position
    ) = 0;


    /** assumes that a molecule is placed somewhere inside the linked cell at index
     *  'linkedCellIndex' and computes the force and potential energy contributions from all molecules
     *  in the same linked cell and the neighboured linked cells onto this molecule.
     *  The molecule "molecule" is considered to NOT be part of the simulation domain and thus the linked cells.
     *  Therefore, another molecule inside the linked cells may even coincide with "molecule".
     *  The results are stored within the molecule.
     */
    virtual void calculateForceAndEnergy(
      coupling::interface::Molecule<3> &molecule
    )
    {
      tarch::la::Vector<3,double> force (0.0);
      double potentialEnergy = 0.0;

      //molecule position
      const tarch::la::Vector<3,double> moleculePosition = molecule.getPosition();
      tarch::la::Vector<3,double> tempMoleculePosition;

      //calculate force
      //find all molecules within cutoff
      double cutoff = ::global_simulation->getcutoffRadius();
      Ensemble* ensemble = ::global_simulation->getEnsemble();
      const double sigma = ensemble->getComponent(0)->getSigma(0);
      const double epsilon = ensemble->getComponent(0)->ljcenter(0).eps();

      const double sigma2 = sigma * sigma;
      const double sigma6 = sigma2 * sigma2 * sigma2;

      double startRegion[] = {moleculePosition[0] - cutoff, moleculePosition[1] - cutoff, moleculePosition[2] - cutoff};
      double endRegion[] = {moleculePosition[0] + cutoff, moleculePosition[1] + cutoff, moleculePosition[2] + cutoff};

      ls1::LS1RegionWrapper region(startRegion, endRegion);
      double cutoff2 = cutoff * cutoff;

      //calculate lennard jones energy
      while(region.iteratorValid())
      {
          ::Molecule* temp = region.getParticleAtIterator();
          tempMoleculePosition = { temp->r(0), temp->r(1), temp->r(2) };
          const auto r = tempMoleculePosition - moleculePosition;
          const double r2 = tarch::la::dot(r, r);
          if(r2 < cutoff2)
          {
              const double r6 = r2 * r2 * r2;
              const auto contrib =  (24.0 * epsilon / r2 * (sigma6/r6) * (1.0 - 2.0 * (sigma6/r6))) * r;
              force += contrib;
          }

          region.iteratorNext();
      }
      //calculate energy (copied from coupling::interface, assuming that the molecule used here is a coupling::datastructures)
      region.iteratorReset();
      while(region.iteratorValid())
        {
            ::Molecule* temp = region.getParticleAtIterator();
            tempMoleculePosition = { temp->r(0), temp->r(1), temp->r(2) };
            const auto r = tempMoleculePosition - moleculePosition;
            const double r2 = tarch::la::dot(r, r);
            if(r2 < cutoff2)
            {
                const double r6 = r2 * r2 * r2;
                const double contrib =  2.0* epsilon * (sigma6/r6) * ((sigma6/r6) - 1.0);
                potentialEnergy += contrib;
            }

            region.iteratorNext();
        }

      molecule.setForce(force);
      molecule.setPotentialEnergy(potentialEnergy);
    }

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
    virtual double getDt() { return ::global_simulation->getIntegrator()->getTimestepLength();}

    /** returns a new molecule iterator for a certain linked cell */
    virtual coupling::interface::MoleculeIterator<ls1::LS1RegionWrapper,3>* getMoleculeIterator(ls1::LS1RegionWrapper& cell)
    {
        return new coupling::interface::LS1MoleculeIterator(cell);
    }
};
#endif