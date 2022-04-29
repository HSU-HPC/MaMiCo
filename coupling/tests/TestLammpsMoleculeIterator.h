// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_TESTS_TESTLAMMPSMOLECULEITERATOR_H_
#define _MOLECULARDYNAMICS_COUPLING_TESTS_TESTLAMMPSMOLECULEITERATOR_H_

#include "coupling/cell-mappings/MoleculeExtractor.h"
#include "coupling/interface/impl/LAMMPS/USER-MAMICO/mamico_cell.h"
#include "coupling/interface/impl/LAMMPS/USER-MAMICO/sorting.h"
#include "coupling/tests/TestLammps.h"
#include <vector>

/** this class is used to test the iterator over molecules in LAMMPS.
 *  We use two input files for molecule coordinates
 * inputpositionsdD_moleculeiterator.xyz where d=2 or d=3, respectively. In 2D,
 * the molecules are sorted into the cells as follows (remember, we use a domain
 * that starts at (1,1,1) and goes to (11,11,11) with cells with dx=2.5):
 *  -------------------------
 *  |  x  | xxx |  x  | xxx |
 *  |     |     |     |     |
 *  -------------------------
 *  |     |     |     |     |
 *  |     | x x |     | x x |
 *  -------------------------
 *  |  x  | xxx |  x  | xxx |
 *  |     |     |     |     |
 *  -------------------------
 *  |     |     |     |     |
 *  |     | x x |     | x x |
 *  -------------------------
 *  Each x marks a molecule (with approx. positions ;-), but correct cell
 * indexing). For 3D, we use the same grid on each z-layer of cells (i.e. we
 * have a 4x4x4 grid where each 4x4x(z=Z) plane corresponds to this picture). We
 * now use the molecule iterator to traverse the atoms and check if we obtain
 * the correct particle numbers, positions, etc.
 *
 *  @author Philipp Neumann
 */
template <unsigned int dim> class TestLammpsMoleculeIterator : public TestLammps<dim> {
public:
  TestLammpsMoleculeIterator(int argc, char** argv, std::string name) : TestLammps<dim>(argc, argv, name) {}
  virtual ~TestLammpsMoleculeIterator() {}

  virtual void run() {
    // initialise all interfaces and simulation parts
    if (dim == 2) {
      TestLammps<dim>::loadLammpsTestConfiguration("inputpositions2D_moleculeiterator.xyz", 24);
    } else {
      TestLammps<dim>::loadLammpsTestConfiguration("inputpositions3D_moleculeiterator.xyz", 96);
    }
    TestLammps<dim>::loadMacroscopicSolverConfiguration();
    TestLammps<dim>::loadMamicoTestConfiguration();

    // extract information on current process -> required for exact evaluation of quantities in respective,
    // subsequently called test methods
    coupling::services::MacroscopicCellServiceImpl<LAMMPS_NS::MamicoCell, dim>* macroscopicCellService =
        (coupling::services::MacroscopicCellServiceImpl<LAMMPS_NS::MamicoCell, dim>*)
            coupling::interface::MamicoInterfaceProvider<LAMMPS_NS::MamicoCell, dim>::getInstance()
                .getMacroscopicCellService();
    if (macroscopicCellService == NULL) {
      std::cout << "ERROR: Could not cast pointer to MacroscopicCellServiceImpl!" << std::endl;
      exit(EXIT_FAILURE);
    }
    const coupling::IndexConversion<dim>& indexConversion = macroscopicCellService->getIndexConversion();
    coupling::datastructures::MacroscopicCells<LAMMPS_NS::MamicoCell, dim>& macroscopicCells = macroscopicCellService->getMacroscopicCells();
    const std::vector<tarch::la::Vector<2, unsigned int>> numberMoleculesPerMacroscopicCell = initNumberMoleculesPerMacroscopicCell(indexConversion);

    // instantiate a molecule extractor
    coupling::cellmappings::MoleculeExtractor<LAMMPS_NS::MamicoCell, dim> moleculeExtractor(
        coupling::interface::MamicoInterfaceProvider<LAMMPS_NS::MamicoCell, dim>::getInstance().getMDSolverInterface());

    // run one simulation time step; this is required so that molecules are sorted into the respective cells once as sort-of start-up phase
    std::cout << "Run one time step..." << std::endl;
    TestLammps<dim>::_lammps->input->one("run 1");
    std::cout << "Check cell sorting via molecule extraction..." << std::endl;

    // test all macroscopic cells and contained molecules on correctness
    testMolecules(macroscopicCells, moleculeExtractor, indexConversion, numberMoleculesPerMacroscopicCell);
  }

private:
  /** returns the pairs (global macroscopic cell index, number molecules in this cell) in a vector. The number is hard-coded, oriented at the setup described
   * above. */
  std::vector<tarch::la::Vector<2, unsigned int>> initNumberMoleculesPerMacroscopicCell(const coupling::IndexConversion<dim>& indexConversion) const {
    std::vector<tarch::la::Vector<2, unsigned int>> moleculesPerMacroscopicCell;
    tarch::la::Vector<2, unsigned int> buffer(0);
    if (dim == 2) {
      buffer[0] = 7;
      buffer[1] = 0;
      moleculesPerMacroscopicCell.push_back(buffer);
      buffer[0] = 9;
      moleculesPerMacroscopicCell.push_back(buffer);
      buffer[0] = 19;
      moleculesPerMacroscopicCell.push_back(buffer);
      buffer[0] = 21;
      moleculesPerMacroscopicCell.push_back(buffer);

      buffer[0] = 8;
      buffer[1] = 2;
      moleculesPerMacroscopicCell.push_back(buffer);
      buffer[0] = 10;
      moleculesPerMacroscopicCell.push_back(buffer);
      buffer[0] = 20;
      moleculesPerMacroscopicCell.push_back(buffer);
      buffer[0] = 22;
      moleculesPerMacroscopicCell.push_back(buffer);

      buffer[0] = 13;
      buffer[1] = 1;
      moleculesPerMacroscopicCell.push_back(buffer);
      buffer[0] = 15;
      moleculesPerMacroscopicCell.push_back(buffer);
      buffer[0] = 25;
      moleculesPerMacroscopicCell.push_back(buffer);
      buffer[0] = 27;
      moleculesPerMacroscopicCell.push_back(buffer);

      buffer[0] = 14;
      buffer[1] = 3;
      moleculesPerMacroscopicCell.push_back(buffer);
      buffer[0] = 16;
      moleculesPerMacroscopicCell.push_back(buffer);
      buffer[0] = 26;
      moleculesPerMacroscopicCell.push_back(buffer);
      buffer[0] = 28;
      moleculesPerMacroscopicCell.push_back(buffer);
    } else {
      for (int layer = 0; layer < 4; layer++) {
        buffer[0] = 43 + layer * 36;
        buffer[1] = 0;
        moleculesPerMacroscopicCell.push_back(buffer);
        buffer[0] = 45 + layer * 36;
        moleculesPerMacroscopicCell.push_back(buffer);
        buffer[0] = 55 + layer * 36;
        moleculesPerMacroscopicCell.push_back(buffer);
        buffer[0] = 57 + layer * 36;
        moleculesPerMacroscopicCell.push_back(buffer);

        buffer[0] = 44 + layer * 36;
        buffer[1] = 2;
        moleculesPerMacroscopicCell.push_back(buffer);
        buffer[0] = 46 + layer * 36;
        moleculesPerMacroscopicCell.push_back(buffer);
        buffer[0] = 56 + layer * 36;
        moleculesPerMacroscopicCell.push_back(buffer);
        buffer[0] = 58 + layer * 36;
        moleculesPerMacroscopicCell.push_back(buffer);

        buffer[0] = 49 + layer * 36;
        buffer[1] = 1;
        moleculesPerMacroscopicCell.push_back(buffer);
        buffer[0] = 51 + layer * 36;
        moleculesPerMacroscopicCell.push_back(buffer);
        buffer[0] = 61 + layer * 36;
        moleculesPerMacroscopicCell.push_back(buffer);
        buffer[0] = 63 + layer * 36;
        moleculesPerMacroscopicCell.push_back(buffer);

        buffer[0] = 50 + layer * 36;
        buffer[1] = 3;
        moleculesPerMacroscopicCell.push_back(buffer);
        buffer[0] = 52 + layer * 36;
        moleculesPerMacroscopicCell.push_back(buffer);
        buffer[0] = 62 + layer * 36;
        moleculesPerMacroscopicCell.push_back(buffer);
        buffer[0] = 64 + layer * 36;
        moleculesPerMacroscopicCell.push_back(buffer);
      }
    }

    return moleculesPerMacroscopicCell;
  }

  /** loops over all macroscopic cells on this process and checks for the number of molecules in each cell. Further, we plot the molecule positions of each
   * cell. */
  void testMolecules(coupling::datastructures::MacroscopicCells<LAMMPS_NS::MamicoCell, dim>& macroscopicCells,
                     coupling::cellmappings::MoleculeExtractor<LAMMPS_NS::MamicoCell, dim>& moleculeExtractor,
                     const coupling::IndexConversion<dim>& indexConversion,
                     const std::vector<tarch::la::Vector<2, unsigned int>>& numberMoleculesPerMacroscopicCell) {
    const tarch::la::Vector<dim, unsigned int> localCells = indexConversion.getLocalNumberMacroscopicCells();
    tarch::la::Vector<dim, unsigned int> loop(0);
    if (dim == 2) {
      for (loop[1] = 1; loop[1] < localCells[1] + 1; loop[1]++) {
        for (loop[0] = 1; loop[0] < localCells[0] + 1; loop[0]++) {
          // determine linearised indices
          const unsigned int localIndex = indexConversion.getLocalCellIndex(loop);
          const unsigned int globalIndex = indexConversion.convertLocalToGlobalCellIndex(localIndex);
          // determine number of molecules in this cell, based on cell structure
          (macroscopicCells.getMacroscopicCellsWithLinkedCells())[localIndex].iterateConstCells(moleculeExtractor);
          const unsigned int numberMoleculesFound = (unsigned int)moleculeExtractor.getExtractedMolecules().size();
          // determine the number of molecules as it was expected
          const unsigned int numberMoleculesExpected = findNumberMolecules(globalIndex, numberMoleculesPerMacroscopicCell);
          if (numberMoleculesFound != numberMoleculesExpected) {
            std::cout << "ERROR TestLammpsMoleculeIterator: Found molecules=" << numberMoleculesFound << ", expected=" << numberMoleculesExpected << std::endl;
            std::cout << "This process: " << indexConversion.getThisProcess() << "; considered (global cell): " << globalIndex << std::endl;
            exit(EXIT_FAILURE);
          }
          std::cout << "Global cell index=" << globalIndex << ", molecule coordinates:" << std::endl;
          for (unsigned int i = 0; i < numberMoleculesFound; i++) {
            std::cout << (moleculeExtractor.getExtractedMolecules())[i] << std::endl;
          }
        }
      }
    } else {
      for (loop[2] = 1; loop[2] < localCells[2] + 1; loop[2]++) {
        for (loop[1] = 1; loop[1] < localCells[1] + 1; loop[1]++) {
          for (loop[0] = 1; loop[0] < localCells[0] + 1; loop[0]++) {
            // determine linearised indices
            const unsigned int localIndex = indexConversion.getLocalCellIndex(loop);
            const unsigned int globalIndex = indexConversion.convertLocalToGlobalCellIndex(localIndex);
            // determine number of molecules in this cell, based on cell structure
            (macroscopicCells.getMacroscopicCellsWithLinkedCells())[localIndex].iterateConstCells(moleculeExtractor);
            const unsigned int numberMoleculesFound = (unsigned int)moleculeExtractor.getExtractedMolecules().size();
            // determine the number of molecules as it was expected
            const unsigned int numberMoleculesExpected = findNumberMolecules(globalIndex, numberMoleculesPerMacroscopicCell);
            if (numberMoleculesFound != numberMoleculesExpected) {
              std::cout << "ERROR TestLammpsMoleculeIterator: Found molecules=" << numberMoleculesFound << ", expected=" << numberMoleculesExpected
                        << std::endl;
              std::cout << "This process: " << indexConversion.getThisProcess() << "; considered (global cell): " << globalIndex << std::endl;
              exit(EXIT_FAILURE);
            }
            std::cout << "Global cell index=" << globalIndex << ", molecule coordinates:" << std::endl;
            for (unsigned int i = 0; i < numberMoleculesFound; i++) {
              std::cout << (moleculeExtractor.getExtractedMolecules())[i] << std::endl;
            }
          }
        }
      }
    }
  }

  /** searches for the number of molecules for the macroscopic cell at global cell index "globalIndex" and returns this number from the vector. If no respective
   *  global cell is found, an error is thrown.
   */
  unsigned int findNumberMolecules(const unsigned int& globalIndex,
                                   const std::vector<tarch::la::Vector<2, unsigned int>>& numberMoleculesPerMacroscopicCell) const {
    const unsigned int size = (unsigned int)numberMoleculesPerMacroscopicCell.size();
    for (unsigned int i = 0; i < size; i++) {
      if (globalIndex == numberMoleculesPerMacroscopicCell[i][0]) {
        return numberMoleculesPerMacroscopicCell[i][1];
      }
    }

    std::cout << "ERROR TestLammpsMoleculeIterator::findNumberMolecules(): Could not find molecule number for global cell index " << globalIndex << std::endl;
    exit(EXIT_FAILURE);
    return 2000;
  }
};

#endif // _MOLECULARDYNAMICS_COUPLING_TESTS_TESTLAMMPSMOLECULEITERATOR_H_
