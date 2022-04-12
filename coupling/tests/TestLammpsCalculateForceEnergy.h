// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _TESTLAMMPSCALCULATEFORCEENERGY_H_
#define _TESTLAMMPSCALCULATEFORCEENERGY_H_

#include "TestLammps.h"
#include "coupling/datastructures/Molecule.h"
#include "mamico_lammps_md_solver_interface.h"

/** tests force and energy computation
 *  @author Philipp Neumann
 */
template <unsigned int dim> class TestLammpsCalculateForceEnergy : public TestLammps<dim> {
public:
  TestLammpsCalculateForceEnergy(int argc, char** argv, std::string name) : TestLammps<dim>(argc, argv, name) {}
  virtual ~TestLammpsCalculateForceEnergy() {}

  virtual void run() {
    // set a molecule into the very middle of the existing molecules
    coupling::datastructures::Molecule<dim> molecule;
    const tarch::la::Vector<dim, double> thisPosition(3.0);
    molecule.setPosition(thisPosition);

    // initialise all interfaces and simulation parts
    if (dim == 2) {
      TestLammps<dim>::loadLammpsTestConfiguration("inputpositionsonly2D.xyz", 4);
    } else {
      TestLammps<dim>::loadLammpsTestConfiguration("inputpositionsonly3D.xyz", 8);
    }
    TestLammps<dim>::loadMacroscopicSolverConfiguration();
    TestLammps<dim>::loadMamicoTestConfiguration();

    // run one time step to sort molecules into linked lists
    LAMMPS_NS::LAMMPS* lammps = TestLammps<dim>::_lammps;
    lammps->input->one("run 1");

    // get MD solver interface (required for sorting)
    LAMMPS_NS::MamicoLammpsMDSolverInterface<dim>* mdSolverInterface =
        (LAMMPS_NS::MamicoLammpsMDSolverInterface<dim>*)coupling::interface::MamicoInterfaceProvider<LAMMPS_NS::MamicoCell, dim>::getInstance()
            .getMDSolverInterface();
    if (mdSolverInterface == NULL) {
      std::cout << "ERROR TestLammpsGhost: could not cast MD Solver interface!" << std::endl;
      exit(EXIT_FAILURE);
    }
    const coupling::IndexConversion<dim>& indexConversion =
        coupling::interface::MamicoInterfaceProvider<LAMMPS_NS::MamicoCell, dim>::getInstance().getMacroscopicCellService()->getIndexConversion();

    // compute LJ parameters (copy-paste from calculateForceAndEnergy)
    const double sigma6 = pow(mdSolverInterface->getMoleculeSigma(), 6.0);
    const double epsilon = mdSolverInterface->getMoleculeEpsilon();
    const double cutOffRadius = 2.5;
    const double cutOffRadiusSquared = cutOffRadius * cutOffRadius;
    const double sigma6OverCutoff6 = sigma6 / (cutOffRadiusSquared * cutOffRadiusSquared * cutOffRadiusSquared);
    const double cutOffEnergy = 4.0 * epsilon * sigma6OverCutoff6 * (sigma6OverCutoff6 - 1.0);
    // set tolerance in measurements
    const double tolerance = 1.0e-8;

    // compute LJ potential energy and force reference values
    const double rij2 = tarch::la::dot(tarch::la::Vector<dim, double>(4.0) - thisPosition, tarch::la::Vector<dim, double>(4.0) - thisPosition);
    const double rij6 = rij2 * rij2 * rij2;
    const double potEnergyRef = pow(2.0, dim) * 0.5 * (4.0 * epsilon * (sigma6 / rij6) * ((sigma6 / rij6) - 1.0) - cutOffEnergy);
    const tarch::la::Vector<dim, double> forceRef(0.0);

    // compute force/energy onto this molecule -> only do this on rank 0, since
    // this is the rank containing the cell with the molecule at 3.0x3.0x3.0
    if (indexConversion.getThisRank() == 0) {
      mdSolverInterface->calculateForceAndEnergy(molecule);

      tarch::la::Vector<dim, double> force = molecule.getForce();
      double potEnergy = molecule.getPotentialEnergy();

      if (fabs(potEnergy - potEnergyRef) > tolerance) {
        std::cout << "ERROR TestLammpsCalculateForceEnergy: potential Energies "
                     "do not match1! Potential energy ref.="
                  << potEnergyRef << ", potential energy=" << potEnergy << std::endl;
        exit(EXIT_FAILURE);
      }
      if (fabs(tarch::la::norm2(force - forceRef) > tolerance)) {
        std::cout << "ERROR TestLammpsCalculcateForceEnergy: forces do not "
                     "match1! Force ref.="
                  << forceRef << ", force=" << force << std::endl;
        exit(EXIT_FAILURE);
      }

      // add molecule to simulation
      std::cout << "Add molecule to MD simulation..." << std::endl;
      mdSolverInterface->addMoleculeToMDSimulation(molecule);
    }

    // synchronise molecules over all processes
    mdSolverInterface->synchronizeMoleculesAfterMassModification();

    // search for this molecule and try to extract potential energy from it;
    // this energy should again return the same value
    if (indexConversion.getThisRank() == 0) {
      std::cout << "Molecule synchronised; do next test..." << std::endl;
      const int nlocal = lammps->atom->nlocal;
      int n = -1;
      for (int i = 0; i < nlocal; i++) {
        bool found = true;
        for (int d = 0; d < dim; d++) {
          found = found && (thisPosition[d] == lammps->atom->x[i][d]);
        }
        if (found) {
          n = i;
          break;
        }
      }
      if (n < 0) {
        std::cout << "ERROR TestLammpsCalculateForceEnergy: Could not find "
                     "molecule after adding it!"
                  << std::endl;
        exit(EXIT_FAILURE);
      } else {
        std::cout << "Molecule added at local position " << n << std::endl;
      }

      // construct a lammps molecule
      LAMMPS_NS::MamicoLammpsMolecule<dim> lammpsMolecule(lammps->atom->x, lammps->atom->v, lammps->atom->f, n, cutOffRadius);
      double potEnergy = lammpsMolecule.getPotentialEnergy();

      if (fabs(potEnergy - potEnergyRef) > tolerance) {
        std::cout << "ERROR TestLammpsCalculateForceEnergy: potential Energies "
                     "do not match2! Potential energy ref.="
                  << potEnergyRef << ", potential energy=" << potEnergy << std::endl;
        exit(EXIT_FAILURE);
      }
    }
  }
};
#endif // _TESTLAMMPSCALCULATEFORCEENERGY_H_
