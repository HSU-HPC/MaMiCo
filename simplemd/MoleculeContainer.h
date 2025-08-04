#ifndef _MOLECULARDYNAMICS_MOLECULARCONTAINER_H_
#define _MOLECULARDYNAMICS_MOLECULARCONTAINER_H_

#include <iostream>
#include <vector>
#include <sstream>
#include <algorithm>
#include <numeric>

#include <Kokkos_Core.hpp>

#include <simplemd/MolecularDynamicsDefinitions.h>
#include <simplemd/Molecule.h>
#include <simplemd/LinkedCell.h>

namespace simplemd {
class MoleculeContainer;
}

class simplemd::MoleculeContainer {
public:
  MoleculeContainer(int numCells, int cellSize)
      : _numCells(numCells), _cellSize(cellSize), moleculeData("moleculeData", numCells, cellSize), linkedCellNumMolecules("linkedCellNumMolecules", numCells) {
  }

  void insert(int cellIdx, Molecule& molecule) {
    moleculeData(cellIdx, linkedCellNumMolecules(cellIdx)) = molecule;
    linkedCellNumMolecules(cellIdx) += 1;
  }

  void remove(int cellIdx, int moleculeIdx) {
    moleculeData(cellIdx, moleculeIdx) = moleculeData(cellIdx, linkedCellNumMolecules(cellIdx) - 1);
    linkedCellNumMolecules(cellIdx) -= 1;
  }

  void clearLinkedCell(int cellIdx) { linkedCellNumMolecules(cellIdx) = 0; }

  void sort(int cellIdx, IndexConverter& indexConverter) { // set all outgoing molecules
    for (size_t i = 0; i < linkedCellNumMolecules(cellIdx); i++) {
      int curMolIdx = indexConverter.getIndex(moleculeData(cellIdx, i).pos);
      if (curMolIdx != cellIdx) { // if molecule does not belong to current cell anymore
        // write data to target end
        moleculeData(curMolIdx, linkedCellNumMolecules(curMolIdx)) = moleculeData(cellIdx, i);
        // increment target end
        linkedCellNumMolecules(curMolIdx)++;
        // delete molecule at own position
        remove(cellIdx, i);
        // decrement iterator as the molecule at position i is now new
        i--;
      }
    }
  }

  void sort(IndexConverter& indexConverter) {
    // find red-black cells
    /**
     * The size<> Vector stores the number of cells, plus ghost layer.
     * If there are (1,1,1) ghost cells per dimension, getLocalIndexOfFirstCell will return (1,1,1)
     * Thus this is multiplied by 2 to account for ghost cells in both locations (begin, end) per axis
     * and then added to local number of cells
     */
    const tarch::la::Vector<MD_DIM, unsigned int> size(getLocalNumberOfCells() + 2u * getLocalIndexOfFirstCell());

// iterate over the domain in a red-black manner
#if (MD_DIM > 2)
    for (unsigned int z = 0; z < 2; z++) {
#endif
#if (MD_DIM > 1)
      for (unsigned int y = 0; y < 2; y++) {
#endif
        for (unsigned int x = 0; x < 2; x++) {
          // determine range/ length of blocks for red-black traversal.
          // For odd block sizes, we need to do some more work in the
          // x/y/z==0-traversals. The second x/y/z==1-traversals are reduced by
          // the normal integer-rounding in this case.
          const tarch::la::Vector<MD_DIM, unsigned int> lengthVector((cellRange[0] + (cellRange[0] % 2) * (x == 0)) / 2
#if (MD_DIM > 1)
                                                                     ,
                                                                     (cellRange[1] + (cellRange[1] % 2) * (y == 0)) / 2
#endif
#if (MD_DIM > 2)
                                                                     ,
                                                                     (cellRange[2] + (cellRange[2] % 2) * (z == 0)) / 2
#endif
          );
          const int length = lengthVector[0]
#if (MD_DIM > 1)
                             * lengthVector[1]
#endif
#if (MD_DIM > 2)
                             * lengthVector[2]
#endif
              ;

          // parallelise loop for all cells that are to be traversed in this way
          Kokkos::parallel_for(
              length, KOKKOS_LAMBDA(const unsigned int j) {
                // compute index of the current cell
                unsigned int index = 0;
#if (MD_DIM > 1)
                int helpIndex1 = j;
                int helpIndex2 = 0;
#endif
                unsigned int coordsCell1Buffer;
                unsigned int coordsCell2Buffer;

#if (MD_DIM > 2)
                // determine plane within traversed block
                helpIndex2 = helpIndex1 / (lengthVector[0] * lengthVector[1]);
                // save rest of index in helpIndex1
                helpIndex1 = helpIndex1 - helpIndex2 * (lengthVector[0] * lengthVector[1]);
                // compute contribution to index
                index += (lowerLeftFrontCell[2] + 2 * helpIndex2 + z) * size[0] * size[1];
#endif
#if (MD_DIM > 1)
                // determine plane within traversed block
                helpIndex2 = helpIndex1 / lengthVector[0];
                // save rest of index in helpIndex1
                helpIndex1 = helpIndex1 - helpIndex2 * lengthVector[0];
                // compute contribution to index
                index += (lowerLeftFrontCell[1] + 2 * helpIndex2 + y) * size[0];
                // compute contribution for last dimension
                index += (lowerLeftFrontCell[0] + 2 * helpIndex1 + x);
#else
        index = lowerLeftFrontCell[0] + 2 * j + x;
#endif
#if (MD_DEBUG == MD_YES)
                std::cout << "Handle cell " << index << std::endl;
#endif
                auto linkedCellLocal(linkedCellNumMolecules);
                auto moleculeDataLocal(moleculeData);
                for (size_t i = 0; i < linkedCellLocal(index); i++) {
                  int curMolIdx = indexConverter.getIndex(moleculeDataLocal(index, i).pos);
                  if (curMolIdx != index) { // if molecule does not belong to current cell anymore
                    // write data to target end
                    moleculeDataLocal(curMolIdx, linkedCellLocal(curMolIdx)) = moleculeDataLocal(index, i);
                    // increment target end
                    linkedCellLocal(curMolIdx)++;
                    // delete molecule at own position
                    moleculeDataLocal(index, i) = moleculeDataLocal(index, linkedCellLocal(index) - 1);
                    linkedCellLocal(index) -= 1;
                    // decrement iterator as the molecule at position i is now new
                    i--;
                  }
                }
              }); // j, Kokkos::parallel_for
        } // x
#if (MD_DIM > 1)
      } // y
#endif
#if (MD_DIM > 2)
    } // z
#endif
  }

  simplemd::Molecule& getMoleculeAt(int i, int j) const { return moleculeData(i, j); }

  simplemd::LinkedCell operator[](int idx) {
    Kokkos::View<Molecule*> lcMoleculeSlice(moleculeData, idx, Kokkos::ALL);
    Kokkos::View<int> lcSizeSlice(linkedCellNumMolecules, idx);
    return simplemd::LinkedCell(lcSizeSlice, lcMoleculeSlice);
  }

  int getNumCells() const { return _numCells; }

  Kokkos::View<simplemd::Molecule**> moleculeData;
  Kokkos::View<int*> linkedCellNumMolecules;

private:
  int _numCells;
  int _cellSize;
};
#endif // _MOLECULARDYNAMICS_MOLECULARCONTAINER_H_
