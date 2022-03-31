// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _TESTLAMMPSGHOST_H_
#define _TESTLAMMPSGHOST_H_

#include "TestLammps.h"
#include "mamico_lammps_md_solver_interface.h"

/** prints the inner and ghost atoms.
 *  @author Philipp Neumann
 */
template <unsigned int dim> class TestLammpsGhost : public TestLammps<dim> {
public:
  TestLammpsGhost(int argc, char **argv, std::string name)
      : TestLammps<dim>(argc, argv, name) {}
  virtual ~TestLammpsGhost() {}

  virtual void run() {
    // initialise all interfaces and simulation parts
    if (dim == 2) {
      TestLammps<dim>::loadLammpsTestConfiguration(
          "inputpositionsonly2D.xyz",
          4); //"inputpositions2D_moleculeiterator.xyz",24);
    } else {
      TestLammps<dim>::loadLammpsTestConfiguration(
          "inputpositionsonly3D.xyz",
          8); //"inputpositions3D_moleculeiterator.xyz",96);
    }
    // extend cut-off radius
    // TestLammps<dim>::_lammps->input->one("pair_coeff 1 1 1.0 1.0 2.5");
    TestLammps<dim>::loadMacroscopicSolverConfiguration();
    TestLammps<dim>::loadMamicoTestConfiguration();

    // get MD solver interface (required for sorting)
    LAMMPS_NS::MamicoLammpsMDSolverInterface<dim> *mdSolverInterface =
        (LAMMPS_NS::MamicoLammpsMDSolverInterface<
            dim> *)coupling::interface::MamicoInterfaceProvider<
            LAMMPS_NS::MamicoCell, dim>::getInstance().getMDSolverInterface();
    if (mdSolverInterface == NULL) {
      std::cout << "ERROR TestLammpsGhost: could not cast MD Solver interface!"
                << std::endl;
      exit(EXIT_FAILURE);
    }
    const coupling::IndexConversion<dim> &indexConversion =
        coupling::interface::MamicoInterfaceProvider<LAMMPS_NS::MamicoCell,
                                                     dim>::getInstance()
            .getMacroscopicCellService()->getIndexConversion();

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    TestLammps<dim>::_lammps->input->one("run 1");
    // print molecules
    for (int i = 0; i < size; i++) {
      MPI_Barrier(MPI_COMM_WORLD);
      if (rank == i) {
        std::cout << "Print molecules in inner cells:" << std::endl;
        mdSolverInterface->printMolecules(
            indexConversion, LAMMPS_NS::Sorting<dim>::PRINT_INNER_CELLS);
        std::cout << "Print molecules in ghost cells: " << std::endl;
        mdSolverInterface->printMolecules(
            indexConversion, LAMMPS_NS::Sorting<dim>::PRINT_GHOST_CELLS);
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }

};
#endif // _TESTLAMMPSGHOST_H_
