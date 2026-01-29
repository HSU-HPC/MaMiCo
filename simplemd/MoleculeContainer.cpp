#include "simplemd/MoleculeContainer.h"
#include "tarch/utils/RandomNumberService.h"
#include <iostream>

simplemd::MoleculeContainer::MoleculeContainer(simplemd::services::ParallelTopologyService& parallelTopologyService, int cellCapacity)
    : _numCells(parallelTopologyService.getLocalNumberOfCells(true)), _ghostCellLayerThickness(parallelTopologyService.getGhostCellLayerThickness()),
      _numLocalCellsNoGhost(_numCells - 2u * _ghostCellLayerThickness), _cellCapacity(cellCapacity), _domainSize(parallelTopologyService.getGlobalDomainSize()),
      _domainOffset(parallelTopologyService.getGlobalDomainOffset()), _meshWidth(parallelTopologyService.getMeshWidth()),
      _globalIndexOfFirstCell(parallelTopologyService.getGlobalIndexOfFirstCell()), _localIndexOfFirstCell(parallelTopologyService.getLocalIndexOfFirstCell()),
      _moleculeData("moleculeData", parallelTopologyService.getLocalNumberOfCellsLinear(true), cellCapacity),
      _linkedCellNumMolecules("linkedCellNumMolecules", parallelTopologyService.getLocalNumberOfCellsLinear(true)),
      _linkedCellIsGhostCell("linkedCellIsGhostCell", _linkedCellNumMolecules.size()) {
  for (unsigned int i = 0; i < _linkedCellIsGhostCell.size(); i++)
    _linkedCellIsGhostCell(i) = isGhostCell(i);
}

void simplemd::MoleculeContainer::insert(unsigned int cellIdx, const simplemd::Molecule& molecule) {
#if (MD_ERROR == MD_YES)
  if (_linkedCellNumMolecules(cellIdx) + 1 > _cellCapacity) {
    Kokkos::printf("Cell capacity=%d would be exceeded by an operation! Exiting...", _cellCapacity);
    Kokkos::abort("simplemd::MoleculeContainer::insert");
  }
#endif
  _moleculeData(cellIdx, _linkedCellNumMolecules(cellIdx)) = molecule;
  _linkedCellNumMolecules(cellIdx) += 1;
}

void simplemd::MoleculeContainer::insert(const simplemd::Molecule& molecule) { insert(positionToCellIndex(molecule.getConstPosition()), molecule); }

void simplemd::MoleculeContainer::remove(unsigned int cellIdx, unsigned int moleculeIdx) {
#if (MD_ERROR == MD_YES)
  if (moleculeIdx >= _linkedCellNumMolecules(cellIdx)) {
    Kokkos::printf("Deleting particle that does not exist! moleculeIdx: %d, cellIdx: %d, num molecules: %d\n", moleculeIdx, cellIdx,
                   _linkedCellNumMolecules(cellIdx));
    Kokkos::abort("ERROR simplemd::MoleculeContainer::remove");
  }
  if (cellIdx >= _linkedCellNumMolecules(cellIdx)) {
    Kokkos::printf("Deleting particle from cell that does not exist! moleculeIdx: %d, cellIdx: %d, num molecules: %d\n", moleculeIdx, cellIdx,
                   _linkedCellNumMolecules(cellIdx));
    Kokkos::abort("ERROR simplemd::MoleculeContainer::remove");
  }
#endif
  _moleculeData(cellIdx, moleculeIdx) = _moleculeData(cellIdx, _linkedCellNumMolecules(cellIdx) - 1);
  _linkedCellNumMolecules(cellIdx) -= 1;
}

void simplemd::MoleculeContainer::clearLinkedCell(unsigned int cellIdx) { _linkedCellNumMolecules(cellIdx) = 0; }

void simplemd::MoleculeContainer::sort(unsigned int cellIdx) { // set all outgoing molecules
  for (size_t i = 0; i < _linkedCellNumMolecules(cellIdx); i++) {
    unsigned int curMolIdx = positionToCellIndex(_moleculeData(cellIdx, i).getPosition());
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
  constexpr unsigned int stride = 3;
  // sort all cells, including ghost
  // find red-black cells

// iterate over the domain in a red-black manner
#if (MD_DIM > 2)
  for (unsigned int z = 0; z < stride; z++) {
#endif
#if (MD_DIM > 1)
    for (unsigned int y = 0; y < stride; y++) {
#endif
      for (unsigned int x = 0; x < stride; x++) {
        // determine range/ length of blocks for red-black traversal.
        // For odd block sizes, we need to do some more work in the
        // x/y/z==0-traversals. The second x/y/z==1-traversals are reduced by
        // the normal integer-rounding in this case.
        const tarch::la::Vector<MD_DIM, unsigned int> lengthVector((_numCells[0]/stride) + (_numCells[0] % stride > 0) * (x == 0)
#if (MD_DIM > 1)
                                                                   ,
                                                                  (_numCells[1]/stride) + (_numCells[1] % stride > 0) * (y == 0)
#endif
#if (MD_DIM > 2)
                                                                   ,
                                                                   (_numCells[2]/stride) + (_numCells[2] % stride > 0) * (z == 0)
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
        printNonGhostCells(1, true, "host start sort");
        Kokkos::parallel_for(
            Kokkos::RangePolicy<MainExecSpace>(0, length), KOKKOS_CLASS_LAMBDA(const unsigned int j) {
              printNonGhostCells(1, j == 0, "device start sort");
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
              index += (stride * helpIndex2 + z) * _numCells[0] * _numCells[1];
#endif
#if (MD_DIM > 1)
              // determine plane within traversed block
              helpIndex2 = helpIndex1 / lengthVector[0];
              // save rest of index in helpIndex1
              helpIndex1 = helpIndex1 - helpIndex2 * lengthVector[0];
              // compute contribution to index
              index += (stride * helpIndex2 + y) * _numCells[0];
              // compute contribution for last dimension
              index += (stride * helpIndex1 + x);
#else
        index = stride * j + x;
#endif
#if (MD_DEBUG == MD_YES)
              Kokkos::printf("Handle cell %u\n", index);
#endif
#if (MD_ERROR == MD_YES)
              if ( index >= _linkedCellNumMolecules.size()){
                Kokkos::abort("simplemd::MoleculeContainer::sort() out-of-bounds access to linked cell");
              }
#endif
              for (size_t i = 0; i < _linkedCellNumMolecules(index); i++) {
                unsigned int curMolIdx = positionToCellIndex(_moleculeData(index, i).getPosition());
                if (curMolIdx != index) { // if molecule does not belong to current cell anymore
                                          // write data to target end
#if (MD_ERROR == MD_YES)
                  if (_linkedCellNumMolecules(curMolIdx) + 1 > _cellCapacity) {
                    Kokkos::printf("Cell capacity=%d would be exceeded by an operation! Exiting...", _cellCapacity);
                    Kokkos::abort("simplemd::MoleculeContainer::insert");
                  }
#endif
                  _moleculeData(curMolIdx, _linkedCellNumMolecules(curMolIdx)) = _moleculeData(index, i);
                  // increment target end
                  _linkedCellNumMolecules(curMolIdx)++;
                  // delete molecule at own position
                  _moleculeData(index, i) = _moleculeData(index, _linkedCellNumMolecules(index) - 1);
                  _linkedCellNumMolecules(index) -= 1;
                  // decrement iterator as the molecule at position i is now new
                  i--;
                }
              }
              printNonGhostCells(1, j == 0, "device end sorf");
            });          // j, Kokkos::parallel_for
        Kokkos::fence(); // Ensure results are available on the host
        printNonGhostCells(1, true, "host end sort");
      } // x
#if (MD_DIM > 1)
    } // y
#endif
#if (MD_DIM > 2)
  } // z
#endif
}

simplemd::Molecule& simplemd::MoleculeContainer::getMoleculeAt(unsigned int i, unsigned int j) const { return _moleculeData(i, j); }

simplemd::LinkedCell simplemd::MoleculeContainer::operator[](unsigned int idx) const {
  return simplemd::LinkedCell(&_moleculeData, &_linkedCellNumMolecules, idx, isGhostCell(idx));
}

simplemd::LinkedCell simplemd::MoleculeContainer::operator[](tarch::la::Vector<MD_DIM, unsigned int> cellIdx) const {
  return (*this)[vectorIndexToLinear(cellIdx)];
}

unsigned int simplemd::MoleculeContainer::getLocalNumberOfCellsScalarWithGhost() const { return _linkedCellNumMolecules.size(); }

unsigned int simplemd::MoleculeContainer::positionToCellIndex(const tarch::la::Vector<MD_DIM, double>& position) const {
  for (unsigned int d = 0; d < MD_DIM; d++) {
#if (MD_ERROR == MD_YES)
    if ((position[d] < _domainOffset[d] - _meshWidth[d]) || (position[d] > _domainOffset[d] + _domainSize[d] + _meshWidth[d])) {
      const char* highlightError = "(<- !!!)";
      Kokkos::printf("Position: %lf%s %lf%s %lf%s\n", position[0], d == 0 ? highlightError : "", position[1], d == 1 ? highlightError : "",
                     MD_DIM > 2 ? position[2] : 0, d == 2 ? highlightError : "");
      Kokkos::abort("ERROR simplemd::MoleculeContainer::positionToCellIndex: Position is out of range!");
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
      Kokkos::printf("index < 0: index=%d; "
                     "Dimension : dim=%u, GIFC=%u, LIFC=%u, cell=%u; "
                     "cellVectorIndex: %lf %lf %lf; "
                     "Position: %lf %lf %lf; "
                     "offset: %lf %lf %lf; "
                     "meshwidth: %lf %lf %lf"
                     "\n",
                     index, d, _globalIndexOfFirstCell[d], _localIndexOfFirstCell[d], (int)(floor((position[d] - _domainOffset[d]) / _meshWidth[d])),
                     cellVectorIndex[0], cellVectorIndex[1], MD_DIM > 2 ? cellVectorIndex[2] : 0, position[0], position[1], MD_DIM > 2 ? position[2] : 0,
                     _domainOffset[0], _domainOffset[1], MD_DIM > 2 ? _domainOffset[2] : 0, _meshWidth[0], _meshWidth[1], MD_DIM > 2 ? _meshWidth[2] : 0);
      Kokkos::abort("ERROR simplemd::MoleculeContainer::positionToCellIndex");
    }
#endif
    cellVectorIndex[d] = (unsigned int)index;
  }
  return vectorIndexToLinear(cellVectorIndex);
}

unsigned int simplemd::MoleculeContainer::vectorIndexToLinear(const tarch::la::Vector<MD_DIM, unsigned int>& vectorIndex) const {
  unsigned int cellLinearIndex = 0;
  unsigned int stepSize = 1;
  for (int d = 0; d < MD_DIM; d++) {
    cellLinearIndex += vectorIndex[d] * stepSize;
    stepSize *= _numCells[d];
  }
  return cellLinearIndex;
}

tarch::la::Vector<MD_DIM, unsigned int> simplemd::MoleculeContainer::getLocalCellIndexVector(const unsigned int cellIndex) const {
  unsigned int help = cellIndex;
  tarch::la::Vector<MD_DIM, unsigned int> localCellIndexVector(0);

  for (int d = MD_DIM - 1; d > -1; d--) {
    unsigned int div = 1;
    for (int e = 0; e < d; e++) {
      div = div * _numCells[e];
    }
    localCellIndexVector[d] = help / div;

    help = help % div;
  }
  return localCellIndexVector;
}

bool simplemd::MoleculeContainer::isGhostCell(const size_t cellIndex) const {
  auto cellLocalVectorIndex = getLocalCellIndexVector(cellIndex);
  for (int d = 0; d < MD_DIM; d++) {
    // if the coordinate is at the beginning or end within the ghost cell layer, return true; otherwise:
    // consider next coordinate
    auto coord = cellLocalVectorIndex[d];
    if ((coord < _ghostCellLayerThickness[d]) || (coord >= _numCells[d] - _ghostCellLayerThickness[d])) {
      return true;
    }
  }
  // return false if no ghost cell coordinate could be detected
  return false;
}

void simplemd::MoleculeContainer::printNonGhostCells(size_t numCells, bool printContainerContents, const char* const label) const {
#if (MD_DUMP == MD_YES)
  if (!printContainerContents) return;
  Kokkos::printf("=== BEGIN DUMP MOLECULE CONTAINER ===\n");
  Kokkos::printf("Label: %s\n", label);
  Kokkos::printf("cell\tpos_x\tpos_y\tpos_z\n");
  size_t linkedCellCount = _linkedCellNumMolecules.size();
  size_t cellsRemaining = numCells == 0 ? linkedCellCount : min(numCells, linkedCellCount);
  for(size_t i = 0; i < linkedCellCount && cellsRemaining > 0; i++) {
    if (_linkedCellIsGhostCell(i)) continue;
    cellsRemaining--;
    auto cellMoleculeCount = _linkedCellNumMolecules(i);
    for(size_t j = 0; j < cellMoleculeCount; j++) {
        Molecule& molecule = getMoleculeAt(i, j);
        Kokkos::printf("%u\t", i);
        auto position = molecule.getConstPosition();
        Kokkos::printf("%lf\t%lf\t%lf\t", position[0], position[1], MD_DIM > 2 ? position[2] : 0);
        auto velocity = molecule.getConstVelocity();
        Kokkos::printf("%lf\t%lf\t%lf\t", velocity[0], velocity[1], MD_DIM > 2 ? velocity[2] : 0);
        auto force = molecule.getConstForce();
        Kokkos::printf("%lf\t%lf\t%lf\t", force[0], force[1], MD_DIM > 2 ? force[2] : 0);
        Kokkos::printf("\n");
    }
  }
  Kokkos::printf("=== END DUMP MOLECULE CONTAINER ===\n");
#endif
}

size_t simplemd::MoleculeContainer::getLocalNumberOfMoleculesWithGhost() const {
  Kokkos::fence(); // Ensure molecule count per cell is up to date
  size_t moleculeCount = 0;
  for (unsigned int i = 0; i < _linkedCellNumMolecules.size(); i++) {
    moleculeCount += _linkedCellNumMolecules(i);
  }
  return moleculeCount;
}

const tarch::la::Vector<MD_DIM, unsigned int>& simplemd::MoleculeContainer::getLocalIndexOfFirstCell() const { return _ghostCellLayerThickness; }
const tarch::la::Vector<MD_DIM, unsigned int> simplemd::MoleculeContainer::getLocalNumberOfCells() const { return _numLocalCellsNoGhost; }
