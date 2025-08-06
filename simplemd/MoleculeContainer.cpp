#include "simplemd/MoleculeContainer.h"

simplemd::MoleculeContainer::MoleculeContainer(simplemd::services::ParallelTopologyService parallelTopologyService, int cellCapacity)
    : _numCells(parallelTopologyService.getLocalNumberOfCells(true)), _cellCapacity(cellCapacity),
#if (MD_ERROR == MD_YES)
      _domainSize(parallelTopologyService.getGlobalDomainSize()),
#endif
      _domainOffset(parallelTopologyService.getGlobalDomainOffset()), _meshWidth(parallelTopologyService.getMeshWidth()),
      _globalIndexOfFirstCell(parallelTopologyService.getGlobalIndexOfFirstCell()), _localIndexOfFirstCell(parallelTopologyService.getLocalIndexOfFirstCell()),
      _moleculeData("moleculeData", parallelTopologyService.getLocalNumberOfCellsLinear(true), cellCapacity),
      _linkedCellNumMolecules("linkedCellNumMolecules", parallelTopologyService.getLocalNumberOfCellsLinear(true)) {
}

void simplemd::MoleculeContainer::insert(int cellIdx, simplemd::Molecule& molecule) {
  _moleculeData(cellIdx, _linkedCellNumMolecules(cellIdx)) = molecule;
  _linkedCellNumMolecules(cellIdx) += 1;
}

void simplemd::MoleculeContainer::insert(simplemd::Molecule& molecule) { insert(positionToCellIndex(molecule.getPosition()), molecule); }

void simplemd::MoleculeContainer::remove(int cellIdx, int moleculeIdx) {
  _moleculeData(cellIdx, moleculeIdx) = _moleculeData(cellIdx, _linkedCellNumMolecules(cellIdx) - 1);
  _linkedCellNumMolecules(cellIdx) -= 1;
}

void simplemd::MoleculeContainer::clearLinkedCell(int cellIdx) { _linkedCellNumMolecules(cellIdx) = 0; }

void simplemd::MoleculeContainer::sort(int cellIdx) { // set all outgoing molecules
  for (size_t i = 0; i < _linkedCellNumMolecules(cellIdx); i++) {
    int curMolIdx = positionToCellIndex(_moleculeData(cellIdx, i).getPosition());
    if (curMolIdx != cellIdx) { // if molecule does not belong to current cell anymore
      // write data to target end
      _moleculeData(curMolIdx, _linkedCellNumMolecules(curMolIdx)) = _moleculeData(cellIdx, i);
      // increment target end
      _linkedCellNumMolecules(curMolIdx)++;
      // delete molecule at own position
      remove(cellIdx, i);
      // decrement iterator as the molecule at position i is now new
      i--;
    }
  }
}

void simplemd::MoleculeContainer::sort() {
  // sort all inner cells, exclude all ghost cells
  // find red-black cells

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
        const tarch::la::Vector<MD_DIM, unsigned int> lengthVector((_numCells[0] + ((_numCells[0]) % 2) * (x == 0)) / 2
#if (MD_DIM > 1)
                                                                   ,
                                                                   (_numCells[1] + ((_numCells[1]) % 2) * (y == 0)) / 2
#endif
#if (MD_DIM > 2)
                                                                   ,
                                                                   (_numCells[2] + ((_numCells[2]) % 2) * (z == 0)) / 2
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

#if (MD_DIM > 2)
              // determine plane within traversed block
              helpIndex2 = helpIndex1 / (lengthVector[0] * lengthVector[1]);
              // save rest of index in helpIndex1
              helpIndex1 = helpIndex1 - helpIndex2 * (lengthVector[0] * lengthVector[1]);
              // compute contribution to index (the starting 1 is the z coordinate of the first cell)
              index += (2 * helpIndex2 + z) * _numCells[0] * _numCells[1];
#endif
#if (MD_DIM > 1)
              // determine plane within traversed block
              helpIndex2 = helpIndex1 / lengthVector[0];
              // save rest of index in helpIndex1
              helpIndex1 = helpIndex1 - helpIndex2 * lengthVector[0];
              // compute contribution to index
              index += (2 * helpIndex2 + y) * _numCells[0];
              // compute contribution for last dimension
              index += (2 * helpIndex1 + x);
#else
        index = 2 * j + x;
#endif
#if (MD_DEBUG == MD_YES)
              std::cout << "Handle cell " << index << std::endl;
#endif
              auto linkedCellLocal(_linkedCellNumMolecules);
              auto moleculeDataLocal(_moleculeData);
              for (size_t i = 0; i < linkedCellLocal(index); i++) {
                unsigned int curMolIdx = positionToCellIndex(moleculeDataLocal(index, i).getPosition());
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

simplemd::Molecule& simplemd::MoleculeContainer::getMoleculeAt(int i, int j) const { return _moleculeData(i, j); }

simplemd::LinkedCell simplemd::MoleculeContainer::operator[](unsigned int idx) { return simplemd::LinkedCell(&_moleculeData, &_linkedCellNumMolecules, idx); }

int simplemd::MoleculeContainer::getNumCells() const { return _linkedCellNumMolecules.size(); }

unsigned int simplemd::MoleculeContainer::positionToCellIndex(const tarch::la::Vector<MD_DIM, double>& position) const {
  for (unsigned int d = 0; d < MD_DIM; d++) {
#if (MD_ERROR == MD_YES)
    if ((position[d] < _domainOffset[d] - _meshWidth[d]) || (position[d] > _domainOffset[d] + _domainSize[d] + _meshWidth[d])) {
      std::cout << "ERROR simplemd::MoleculeContainer::positionToCellIndex: Position ";
      std::cout << d << " is out of range!" << std::endl;
      std::cout << "Position: " << position << std::endl;
      exit(EXIT_FAILURE);
    }
#endif
  }
  tarch::la::Vector<MD_DIM, unsigned int> cellVectorIndex(0);

  // determine current cell index (in serial, i.e. 1-D, form)
  for (unsigned int d = 0; d < MD_DIM; d++) {
    // find global cell index

    int index = (int)(floor((position[d] - _domainOffset[d]) / _meshWidth[d]));
    // shift into local cell index
    index += _localIndexOfFirstCell[d];
    index -= _globalIndexOfFirstCell[d];
#if (MD_ERROR == MD_YES)
    if (index < 0) {
      std::cout << "ERROR simplemd::MoleculeContainer::positionToCellIndex: index < 0: index=";
      std::cout << index << std::endl;
      std::cout << "Dimension : " << d << "," << _globalIndexOfFirstCell[d] << "," << _localIndexOfFirstCell[d] << std::endl;
      std::cout << (int)(floor((position[d] - _domainOffset[d]) / _meshWidth[d])) << std::endl;
      for (unsigned int e = 0; e < d; e++) {
        std::cout << cellVectorIndex[e] << std::endl;
      }
      std::cout << "Position: " << position << ", offset: " << _domainOffset << ", meshwidth: " << _meshWidth << std::endl;
      exit(EXIT_FAILURE);
    }
#endif
    cellVectorIndex[d] = (unsigned int)index;
  }
  return vectorIndexToLinear(cellVectorIndex);
}

const unsigned int simplemd::MoleculeContainer::vectorIndexToLinear(const tarch::la::Vector<MD_DIM, unsigned int>& vectorIndex) const {
  unsigned int cellLinearIndex = 0;
  unsigned int stepSize = 1;
  for (int d = 0; d < MD_DIM; d++) {
    cellLinearIndex += vectorIndex[d] * stepSize;
    stepSize *= _numCells[d];
  }
  return cellLinearIndex;
}
