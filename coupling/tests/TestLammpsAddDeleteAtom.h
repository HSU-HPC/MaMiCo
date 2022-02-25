// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_TESTS_TESTLAMMPSADDDELETEATOM_H_
#define _MOLECULARDYNAMICS_COUPLING_TESTS_TESTLAMMPSADDDELETEATOM_H_

#include "TestLammps.h"
#include "coupling/datastructures/Molecule.h"
#include "tarch/la/ScalarOperations.h"

template <unsigned int dim>
class TestLammpsAddDeleteAtom : public TestLammps<dim> {
public:
  TestLammpsAddDeleteAtom(int argc, char **argv, std::string name)
      : TestLammps<dim>(argc, argv, name) {}
  virtual ~TestLammpsAddDeleteAtom() {}

  virtual void run() {
    // LOAD COUPLED SCENARIO
    // ----------------------------------------------------------------------
    if (dim == 2) {
      TestLammps<dim>::loadLammpsTestConfiguration("inputpositionsonly2D.xyz",
                                                   4);
    } else {
      TestLammps<dim>::loadLammpsTestConfiguration("inputpositionsonly3D.xyz",
                                                   8);
    }
    TestLammps<dim>::loadMacroscopicSolverConfiguration();
    TestLammps<dim>::loadMamicoTestConfiguration();
    std::cout << "All configs, solvers and coupling loaded and initialised..."
              << std::endl;

    // DEFINE PROPERTIES OF ATOM TO BE ADDED/DELETED
    // ---------------------------------------------- helper variables
    const tarch::la::Vector<dim, unsigned int> linkedCellsPerMacroscopicCell(1);
    const tarch::la::Vector<dim, unsigned int> linkedCellInMacroscopicCell(0);
    const coupling::IndexConversion<dim> &indexConversion =
        TestLammps<dim>::_macroscopicCellService->getIndexConversion();
    const double tolerance =
        1.0e-8; // this value should be the same as tolerance defined in
                // MDSolverInterface of LAMMPS
    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // position and dummy values of atom to be deleted
    tarch::la::Vector<dim, double> posDeleteAtom(0.0);
    posDeleteAtom[0] = 4.0;
    posDeleteAtom[1] = 2.0;
    if (dim != 2)
      posDeleteAtom[2] = 2.0;
    const tarch::la::Vector<dim, double> dummyVec(0.0);
    const double dummyEnergy = 0.0;
    const coupling::datastructures::Molecule<dim> deleteMolecule(
        posDeleteAtom, dummyVec, dummyVec, dummyEnergy);
    tarch::la::Vector<dim, unsigned int> deleteCellIndex(1);
    deleteCellIndex[0] = 2; // global cell index for deletion
    // true, if this is the respective rank that should contain this molecule
    const bool performDeletion = (size == 1) || ((size == 4) && (rank == 0)) ||
                                 ((size == 16) && (rank == 1)) ||
                                 ((size == 8) && (rank == 0)) ||
                                 ((size == 64) && (rank == 1));

    // position, velocity, force of atom to be inserted
    tarch::la::Vector<dim, double> posAddAtom(3.6);
    const tarch::la::Vector<dim, double> velAddAtom(5.0);
    const tarch::la::Vector<dim, double> fAddAtom(10.0);
    const coupling::datastructures::Molecule<dim> addMolecule(
        posAddAtom, velAddAtom, fAddAtom, dummyEnergy);
    tarch::la::Vector<dim, unsigned int> addCellIndex(
        2); // global cell index for adding atom
    // true, if this rank holds the respective domain containing the new
    // molecule's position
    const bool performAdd = (size == 1) || ((size == 4) && (rank == 0)) ||
                            ((size == 16) && (rank == 5)) ||
                            ((size == 8) && (rank == 0)) ||
                            ((size == 64) && (rank == 21));

    // set a very small time step (so that molecules hardly move and positions
    // thus remain quasi identical) and run one time step to sort molecules into
    // mamico cells
    // TestLammps<dim>::_lammps->input->one("timestep 1e-14");
    TestLammps<dim>::_lammps->input->one("run 1");

    // for testing purposes, set cut-off radius to 2.5 (=mamico-cell-size); we
    // can then check whether the newly inserted/deleted atoms are transferred
    // to the neighbouring processes in the
    // synchronizeMoleculesAfterMassModification()-method.
    // TestLammps<dim>::_lammps->input->one("pair_style lj/cut 2.5");
    // TestLammps<dim>::_lammps->input->one("pair_coeff 1 1 1.0 1.0 2.5");

    // print out all molecules for each rank (synchronised)
    TestLammps<dim>::printMolecules();

    // PERFORM DELETION TEST
    // -----------------------------------------------------------------------------------
    // delete atom
    if (performDeletion) {
      const int numberAtoms = TestLammps<dim>::_lammps->atom->nlocal;
      // convert global to local coordinate
      deleteCellIndex =
          indexConversion.convertGlobalToLocalVectorCellIndex(deleteCellIndex);
      LAMMPS_NS::MamicoCell &deleteCell =
          coupling::interface::MamicoInterfaceProvider<LAMMPS_NS::MamicoCell,
                                                       dim>::getInstance()
              .getMDSolverInterface()
              ->getLinkedCell(deleteCellIndex, linkedCellInMacroscopicCell,
                              linkedCellsPerMacroscopicCell, indexConversion);
      coupling::interface::MoleculeIterator<LAMMPS_NS::MamicoCell,
                                            dim> *iterator =
          coupling::interface::MamicoInterfaceProvider<LAMMPS_NS::MamicoCell,
                                                       dim>::getInstance()
              .getMDSolverInterface()
              ->getMoleculeIterator(deleteCell);
      int cellCounterBeforeDeletion = 0;
      int cellCounterAfterDeletion = 0;

      // determine number of molecules in this cell before deletion and plot
      // molecule positions
      std::cout << "Molecules before deletion:" << std::endl;
      for (iterator->begin(); iterator->continueIteration(); iterator->next()) {
        std::cout << iterator->getConst().getPosition() << std::endl;
        cellCounterBeforeDeletion++;
      };

      // delete molecule
      coupling::interface::MamicoInterfaceProvider<LAMMPS_NS::MamicoCell,
                                                   dim>::getInstance()
          .getMDSolverInterface()
          ->deleteMoleculeFromMDSimulation(deleteMolecule, deleteCell);

      // check if number of atoms has been decremented
      if (numberAtoms - 1 != TestLammps<dim>::_lammps->atom->nlocal) {
        std::cout << "ERROR TestLammpsAddDeleteAtom::run(): Same number of "
                     "atoms after deletion as before!"
                  << std::endl;
        exit(EXIT_FAILURE);
      }
      // check if deleted atom can be found in mamico cell
      std::cout << "Molecules after deletion:" << std::endl;
      for (iterator->begin(); iterator->continueIteration(); iterator->next()) {
        std::cout << iterator->getConst().getPosition() << std::endl;
        cellCounterAfterDeletion++;
      }
      if (cellCounterAfterDeletion != cellCounterBeforeDeletion - 1) {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        std::cout << "Rank " << rank
                  << ": ERROR TestLammpsAddDeleteAtom::run(): Did not reduce "
                     "number of molecules in the respective mamico cell!"
                  << std::endl;
        exit(EXIT_FAILURE);
      }

      // check if deleted atom can be found anywhere in LAMMPS
      for (int i = 0; i < TestLammps<dim>::_lammps->atom->nlocal; i++) {
        bool found = true;
        for (unsigned int d = 0; d < dim; d++) {
          found = found &&
                  tarch::la::equals(TestLammps<dim>::_lammps->atom->x[i][d],
                                    posDeleteAtom[d], tolerance);
        }
        if (found) {
          std::cout << "ERROR TestLammpsAddDeleteAtom::run(): Found atom in "
                       "LAMMPS though it should have been deleted!"
                    << std::endl;
          exit(EXIT_FAILURE);
        }
      }

      delete iterator;
    } // finish deletion tests

    if (rank == 0) {
      std::cout << "Print molecules after deletion test..." << std::endl;
    }
    TestLammps<dim>::printMolecules();

    // TEST ADDING A MOLECULE
    // ---------------------------------------------------------
    if (performAdd) {
      const int numberAtoms = TestLammps<dim>::_lammps->atom->nlocal;
      // convert global to local cell index
      addCellIndex =
          indexConversion.convertGlobalToLocalVectorCellIndex(addCellIndex);
      LAMMPS_NS::MamicoCell &addCell =
          coupling::interface::MamicoInterfaceProvider<LAMMPS_NS::MamicoCell,
                                                       dim>::getInstance()
              .getMDSolverInterface()
              ->getLinkedCell(addCellIndex, linkedCellInMacroscopicCell,
                              linkedCellsPerMacroscopicCell, indexConversion);
      coupling::interface::MoleculeIterator<LAMMPS_NS::MamicoCell,
                                            dim> *iterator =
          coupling::interface::MamicoInterfaceProvider<LAMMPS_NS::MamicoCell,
                                                       dim>::getInstance()
              .getMDSolverInterface()
              ->getMoleculeIterator(addCell);
      int cellCounterBeforeAdding = 0;
      int cellCounterAfterAdding = 0;

      // determine number of molecules in this cell before adding new one and
      // plot molecule positions
      std::cout << "Molecules before adding:" << std::endl;
      for (iterator->begin(); iterator->continueIteration(); iterator->next()) {
        std::cout << iterator->getConst().getPosition() << ", "
                  << iterator->getConst().getVelocity() << ", "
                  << iterator->getConst().getForce() << std::endl;
        cellCounterBeforeAdding++;
      };

      // add molecule
      coupling::interface::MamicoInterfaceProvider<LAMMPS_NS::MamicoCell,
                                                   dim>::getInstance()
          .getMDSolverInterface()
          ->addMoleculeToMDSimulation(addMolecule);

      // check if number of atoms is now equal again
      if (numberAtoms + 1 != TestLammps<dim>::_lammps->atom->nlocal) {
        std::cout << "ERROR TestLammpsAddDeleteAtom::run(): Same number of "
                     "atoms after adding new molecule as before!"
                  << std::endl;
        exit(EXIT_FAILURE);
      }
      // check if deleted atom can be found in mamico cell
      std::cout << "Molecules after adding new one:" << std::endl;
      for (iterator->begin(); iterator->continueIteration(); iterator->next()) {
        std::cout << iterator->getConst().getPosition() << ", "
                  << iterator->getConst().getVelocity() << ", "
                  << iterator->getConst().getForce() << std::endl;
        cellCounterAfterAdding++;
      }
      if (cellCounterAfterAdding != cellCounterBeforeAdding + 1) {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        std::cout << "Rank " << rank
                  << ": ERROR TestLammpsAddDeleteAtom::run(): Did not increase "
                     "number of molecules in the respective mamico cell!"
                  << std::endl;
        exit(EXIT_FAILURE);
      }

      // check if new atom can be found in LAMMPS
      bool foundMolecule = false;
      for (int i = 0; i < TestLammps<dim>::_lammps->atom->nlocal; i++) {
        bool found = true;
        for (unsigned int d = 0; d < dim; d++) {
          found = found &&
                  tarch::la::equals(TestLammps<dim>::_lammps->atom->x[i][d],
                                    posAddAtom[d], tolerance);
          found = found &&
                  tarch::la::equals(TestLammps<dim>::_lammps->atom->v[i][d],
                                    velAddAtom[d], tolerance);
          found = found &&
                  tarch::la::equals(TestLammps<dim>::_lammps->atom->f[i][d],
                                    fAddAtom[d], tolerance);
        }
        if (found) {
          foundMolecule = true;
          break;
        }
      }
      if (!foundMolecule) {
        std::cout << "ERROR TestLammpsAddDeleteAtom::run(): Could not find new "
                     "atom in LAMMPS!"
                  << std::endl;
        exit(EXIT_FAILURE);
      }

      delete iterator;
    } // finish adding molecule test

    // test synchronisation of mass across processes
    {
      const int totalNumberAtoms = TestLammps<dim>::_lammps->atom->natoms;
      if (dim == 2) {
        if (totalNumberAtoms != 4) {
          std::cout << "ERROR TestLammpsAddDeleteAtom: total number of atoms="
                    << totalNumberAtoms << " instead of 4" << std::endl;
          exit(EXIT_FAILURE);
        }
      }
      if (dim == 3) {
        if (totalNumberAtoms != 8) {
          std::cout << "ERROR TestLammpsAddDeleteAtom: total number of atoms="
                    << totalNumberAtoms << " instead of 8" << std::endl;
          exit(EXIT_FAILURE);
        }
      }
      coupling::interface::MamicoInterfaceProvider<LAMMPS_NS::MamicoCell,
                                                   dim>::getInstance()
          .getMDSolverInterface()
          ->synchronizeMoleculesAfterMassModification();
      const int totalNumberAtomsAfterInsertion =
          TestLammps<dim>::_lammps->atom->natoms;
      if (totalNumberAtomsAfterInsertion != totalNumberAtoms) {
        std::cout << "ERROR TestLammpsAddDeleteAtom: total number of atoms has "
                     "changed; new total number of atoms: "
                  << totalNumberAtomsAfterInsertion << "!" << std::endl;
        exit(EXIT_FAILURE);
      }
    }

    if (rank == 0) {
      std::cout << "Print molecules after insertion test..." << std::endl;
    }
    TestLammps<dim>::printMolecules();

    // do one more time step
    TestLammps<dim>::_lammps->input->one("run 1");
  }
};

#endif // _MOLECULARDYNAMICS_COUPLING_TESTS_TESTLAMMPSADDDELETEATOM_H_
