#pragma once

#include "coupling/CouplingMDDefinitions.h"
#include "coupling/IndexConversion.h"
#include "coupling/datastructures/MacroscopicCell.h"
#include "coupling/indexing/CellIndex.h"
#include "coupling/interface/MacroscopicSolverInterface.h"
#include "precice/SolverInterface.hpp"
#include "tarch/la/Vector.h"
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
#include <mpi.h>
#endif
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>

namespace coupling {
namespace solvers {
template <unsigned int dim> class PreciceAdapter;
template <unsigned int dim> class PreciceInterface;
} // namespace solvers
} // namespace coupling

/**
 * Adapter for the preCICE library
 * 'Vertices' are used in the preCICE context (elements constituting the preCICE mesh)
 * 'Cells' are used in the MaMiCo context (elements consituting the MaMiCo cartesian grid) 
 */
template <unsigned int dim> class coupling::solvers::PreciceAdapter {
public:
  PreciceAdapter(const std::string M2mMeshName, const std::string m2MMeshName, const std::string M2mVelocityName, const std::string m2MVelocityName)
      : _M2mMeshName(M2mMeshName), _m2MMeshName(m2MMeshName), _M2mVelocityName(M2mVelocityName), _m2MVelocityName(m2MVelocityName), _rank(0) {
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
    MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
#endif
    int size = 1;
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
    MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
    _solverInterface = new precice::SolverInterface("mamico", "../precice-config.xml", _rank, size);
  }

  virtual ~PreciceAdapter() {
    if (_solverInterface != NULL) {
      _solverInterface->finalize();
      delete _solverInterface;
      _solverInterface = NULL;
    }
    if (_M2mVertexCoords != NULL) {
      delete[] _M2mVertexCoords;
    }
    if (_M2mVertexVelocities != NULL) {
      delete[] _M2mVertexVelocities;
    }
    if (_M2mVertexIndices != NULL) {
      delete[] _M2mVertexIndices;
    }
    if (_m2MVertexCoords != NULL) {
      delete[] _m2MVertexCoords;
    }
    if (_m2MVertexVelocities != NULL) {
      delete[] _m2MVertexVelocities;
    }
    if (_m2MVertexIndices != NULL) {
      delete[] _m2MVertexIndices;
    }
  }

  void setMeshes(coupling::solvers::PreciceInterface<dim>* _macroscopicSolverInterface, const tarch::la::Vector<3, double> mdDomainOffset, const tarch::la::Vector<3, double> macroscopicCellSize) {
    if (_rank == 0) {
      if (_solverInterface->hasMesh(_M2mMeshName)) {
        std::vector<double> M2mVertexCoords;
        using namespace coupling::indexing;
        for (CellIndex<dim, IndexTrait::vector> cellIndex_v : CellIndex<dim, IndexTrait::vector>()) {
          if (_macroscopicSolverInterface->sendMacroscopicQuantityToMDSolver(static_cast<tarch::la::Vector<dim, unsigned int>>(cellIndex_v.get()))) {
            for (unsigned int currentDim = 0; currentDim < dim; currentDim++) {
              M2mVertexCoords.push_back(mdDomainOffset[currentDim] + cellIndex_v.get()[currentDim] * macroscopicCellSize[currentDim] -
                                      macroscopicCellSize[currentDim] + 0.5 * macroscopicCellSize[currentDim]);
            }
          }
        }
        _M2mVertexNumbers = M2mVertexCoords.size() / dim;
        _M2mVertexCoords = new double[M2mVertexCoords.size()];
        std::copy(M2mVertexCoords.begin(), M2mVertexCoords.end(), _M2mVertexCoords);
        _M2mVertexIndices = new int[_M2mVertexNumbers];
        _solverInterface->setMeshVertices(_solverInterface->getMeshID(_M2mMeshName), _M2mVertexNumbers, _M2mVertexCoords, _M2mVertexIndices);
        _M2mVertexVelocities = new double[_M2mVertexNumbers * dim];
      }
      if (_solverInterface->hasMesh(_m2MMeshName)) {
        std::vector<double> m2MVertexCoords;
        std::vector<unsigned int> m2MCellIndices;
        using namespace coupling::indexing;
        for (CellIndex<dim, IndexTrait::vector> cellIndex_v : CellIndex<dim, IndexTrait::vector>()) {
          if (_macroscopicSolverInterface->receiveMacroscopicQuantityFromMDSolver(static_cast<tarch::la::Vector<dim, unsigned int>>(cellIndex_v.get()))) {
            coupling::indexing::CellIndex<dim> cellIndex_s = cellIndex_v;
            m2MCellIndices.push_back(cellIndex_s.get());
            for (unsigned int currentDim = 0; currentDim < dim; currentDim++) {
              m2MVertexCoords.push_back(mdDomainOffset[currentDim] + cellIndex_v.get()[currentDim] * macroscopicCellSize[currentDim] -
                                      macroscopicCellSize[currentDim] + 0.5 * macroscopicCellSize[currentDim]);
            }
          }
        }
        _m2MVertexNumbers = m2MVertexCoords.size() / dim;
        _m2MVertexCoords = new double[m2MVertexCoords.size()];
        std::copy(m2MVertexCoords.begin(), m2MVertexCoords.end(), _m2MVertexCoords);
        _m2MVertexIndices = new int[_m2MVertexNumbers];
        _solverInterface->setMeshVertices(_solverInterface->getMeshID(_m2MMeshName), _m2MVertexNumbers, _m2MVertexCoords, _m2MVertexIndices);
        _m2MVertexVelocities = new double[_m2MVertexNumbers * dim];

        for (size_t i = 0; i < _m2MVertexNumbers; i++) {
          std::vector<int> tetrahedronsVertexIndices = getTetrahedronsVertexIndices(i, _m2MVertexIndices, m2MCellIndices);
          for (size_t tetrahedronIndex = 0; tetrahedronIndex < tetrahedronsVertexIndices.size() / 4; tetrahedronIndex++) {
            _solverInterface->setMeshTetrahedron(_solverInterface->getMeshID(_m2MMeshName), tetrahedronsVertexIndices[tetrahedronIndex * 4],
                                                tetrahedronsVertexIndices[tetrahedronIndex * 4 + 1], tetrahedronsVertexIndices[tetrahedronIndex * 4 + 2],
                                                tetrahedronsVertexIndices[tetrahedronIndex * 4 + 3]);
          }
        }
      }
    }
  }

  double initialize() { return _solverInterface->initialize(); }

  bool isCouplingOngoing() { return _solverInterface->isCouplingOngoing(); }

  double advance(const double dt) { return _solverInterface->advance(dt); }

  void readData() {
    if (_solverInterface->hasMesh(_M2mMeshName)) {
      _solverInterface->readBlockVectorData(_solverInterface->getDataID(_M2mVelocityName, _solverInterface->getMeshID(_M2mMeshName)), _M2mVertexNumbers,
                                            _M2mVertexIndices, _M2mVertexVelocities);
    }
  }

  void writeData() {
    if (_solverInterface->hasMesh(_m2MMeshName)) {
      _solverInterface->writeBlockVectorData(_solverInterface->getDataID(_m2MVelocityName, _solverInterface->getMeshID(_m2MMeshName)), _m2MVertexNumbers,
                                             _m2MVertexIndices, _m2MVertexVelocities);
    }
  }

  // Couette solver methods, to be moved.

  tarch::la::Vector<3, double> getVelocity(const tarch::la::Vector<3, double> &pos) const {
    tarch::la::Vector<3, double> vel(0.0);
    if (_solverInterface->hasMesh(_M2mMeshName)) {
      unsigned int vertexIndex = 0;
      while (vertexIndex < _M2mVertexNumbers &&
             !(_M2mVertexCoords[dim * vertexIndex] == pos[0] && _M2mVertexCoords[dim * vertexIndex + 1] == pos[1] && _M2mVertexCoords[dim * vertexIndex + 2] == pos[2]))
        vertexIndex++;
      if (vertexIndex < _M2mVertexNumbers) {
        for (unsigned int currentDim = 0; currentDim < dim; currentDim++) {
          vel[currentDim] = _M2mVertexVelocities[dim * vertexIndex + currentDim];
        }
      }
    }
    return vel;
  }

  void setVelocity(const tarch::la::Vector<3, double> &pos, const tarch::la::Vector<3, double> &velocity) const {
    if (_solverInterface->hasMesh(_m2MMeshName)) {
      unsigned int vertexIndex = 0;
      while (vertexIndex < _m2MVertexNumbers &&
             !(_m2MVertexCoords[dim * vertexIndex] == pos[0] && _m2MVertexCoords[dim * vertexIndex + 1] == pos[1] && _m2MVertexCoords[dim * vertexIndex + 2] == pos[2]))
        vertexIndex++;
      if (vertexIndex < _m2MVertexNumbers) {
        for (unsigned int currentDim = 0; currentDim < dim; currentDim++) {
          _m2MVertexVelocities[dim * vertexIndex + currentDim] = velocity[currentDim];
        }
      }
    }
  }

  // void setMDBoundaryValues(std::vector<coupling::datastructures::MacroscopicCell<3>*>& m2MBuffer, const unsigned int* const m2MCellIndices) {
  //   if (_solverInterface->hasMesh(_m2MMeshName)) {
  //     for (size_t i = 0; i < _m2MVertexNumbers; i++) {
  //       if (m2MBuffer[i]->getMacroscopicMass() != 0.0) {
  //         tarch::la::Vector<3, double> vel((1.0 / m2MBuffer[i]->getMacroscopicMass()) * m2MBuffer[i]->getMacroscopicMomentum());
  //         for (unsigned int currentDim = 0; currentDim < dim; currentDim++)
  //           _m2MVertexVelocities[dim * i + currentDim] = vel[currentDim];
  //       }
  //     }
  //   }
  // }

private:
  std::vector<int> getTetrahedronsVertexIndices(int i, int* m2MVertexIndices, std::vector<unsigned int>& m2MCellIndices) {
    int m2MVertexIndex = m2MVertexIndices[i];
    using CellIndex_v = coupling::indexing::CellIndex<dim, coupling::indexing::IndexTrait::vector>;
    using CellIndex_s = coupling::indexing::CellIndex<dim>;
    CellIndex_v cellIndex_v = CellIndex_s{m2MCellIndices[i]};
    std::vector<int> tetrahedronsVertexIndices;
    const int numberOfNeighbors = 7;
    tarch::la::Vector<3, int> directions[numberOfNeighbors] = {{1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}, {0, 0, 1}};
    int vertexIndices[numberOfNeighbors]; 
    for (size_t i = 0; i < numberOfNeighbors; i++) {
      CellIndex_v neighborCellIndex_v{cellIndex_v.get() + directions[i]};
      CellIndex_s neighborCellIndex_s = neighborCellIndex_v;
      std::vector<unsigned int>::iterator it = std::find(m2MCellIndices.begin(), m2MCellIndices.end(), neighborCellIndex_s.get());
      if (it != m2MCellIndices.end()) {
        vertexIndices[i] = m2MVertexIndices[it - m2MCellIndices.begin()];
      } else {
        vertexIndices[i] = -1;
      }
    }
    if (vertexIndices[0] > -1 && vertexIndices[5] > -1 && vertexIndices[6] > -1)
      tetrahedronsVertexIndices.insert(tetrahedronsVertexIndices.end(), {m2MVertexIndex, vertexIndices[0], vertexIndices[5], vertexIndices[6]});
    if (vertexIndices[0] > -1 && vertexIndices[2] > -1 && vertexIndices[5] > -1)
      tetrahedronsVertexIndices.insert(tetrahedronsVertexIndices.end(), {m2MVertexIndex, vertexIndices[0], vertexIndices[2], vertexIndices[5]});
    if (vertexIndices[0] > -1 && vertexIndices[5] > -1 && vertexIndices[6] > -1 && vertexIndices[3] > -1)
      tetrahedronsVertexIndices.insert(tetrahedronsVertexIndices.end(), {vertexIndices[0], vertexIndices[5], vertexIndices[6], vertexIndices[3]});
    if (vertexIndices[0] > -1 && vertexIndices[1] > -1 && vertexIndices[2] > -1 && vertexIndices[5] > -1)
      tetrahedronsVertexIndices.insert(tetrahedronsVertexIndices.end(), {vertexIndices[0], vertexIndices[1], vertexIndices[2], vertexIndices[5]});
    if (vertexIndices[0] > -1 && vertexIndices[1] > -1 && vertexIndices[4] > -1 && vertexIndices[5] > -1)
      tetrahedronsVertexIndices.insert(tetrahedronsVertexIndices.end(), {vertexIndices[0], vertexIndices[1], vertexIndices[4], vertexIndices[5]});
    if (vertexIndices[0] > -1 && vertexIndices[3] > -1 && vertexIndices[4] > -1 && vertexIndices[5] > -1)
      tetrahedronsVertexIndices.insert(tetrahedronsVertexIndices.end(), {vertexIndices[0], vertexIndices[3], vertexIndices[4], vertexIndices[5]});
    return tetrahedronsVertexIndices;
  }

  const std::string _M2mMeshName;
  const std::string _m2MMeshName;
  const std::string _M2mVelocityName;
  const std::string _m2MVelocityName;

  int _rank;

  precice::SolverInterface* _solverInterface = nullptr;

  int* _M2mVertexIndices = nullptr;
  double* _M2mVertexCoords = nullptr;
  unsigned int _M2mVertexNumbers = 0;
  double* _M2mVertexVelocities = nullptr;

  int* _m2MVertexIndices = nullptr;
  double* _m2MVertexCoords = nullptr;
  unsigned int _m2MVertexNumbers = 0;
  double* _m2MVertexVelocities = nullptr;
};

template <unsigned int dim> class coupling::solvers::PreciceInterface : public coupling::interface::MacroscopicSolverInterface<dim> {
public:
  PreciceInterface(const tarch::la::Vector<dim, int> globalNumberMacroscopicCells, const unsigned int overlap, const unsigned int rank,
  const coupling::IndexConversion<dim>* indexConversion)
      : _globalNumberMacroscopicCells(globalNumberMacroscopicCells), _overlap(overlap), _rank(rank), _indexConversion(indexConversion) {}

  bool receiveMacroscopicQuantityFromMDSolver(tarch::la::Vector<dim, unsigned int> globalCellIndex) override {
    // if (_indexConversion->getUniqueRankForMacroscopicCell(globalCellIndex) != _rank) return false;
    bool rcv = true;
    for (unsigned int currentDim = 0; currentDim < dim; currentDim++) {
      rcv &= globalCellIndex[currentDim] >= 1 + (_overlap - 1);
      rcv &= globalCellIndex[currentDim] < _globalNumberMacroscopicCells[currentDim] + 1 - (_overlap - 1);
    }
    return rcv;
  }

  bool sendMacroscopicQuantityToMDSolver(tarch::la::Vector<3, unsigned int> globalCellIndex) override {
    // std::vector<unsigned int> ranks = _indexConversion->getRanksForMacroscopicCell(globalCellIndex);
    // if (std::find(ranks.begin(), ranks.end(), _rank) == ranks.end()) return false;
    bool isGhostCell = false;
    bool isInner = true;
    for (unsigned int currentDim = 0; currentDim < dim; currentDim++) {
      isGhostCell |= globalCellIndex[currentDim] > _globalNumberMacroscopicCells[currentDim];
      isGhostCell |= globalCellIndex[currentDim] < 1;
      isInner &= globalCellIndex[currentDim] >= 1 + _overlap;
      isInner &= globalCellIndex[currentDim] < _globalNumberMacroscopicCells[currentDim] + 1 - _overlap;
    }
    return (!isGhostCell) && (!isInner);
  }

  std::vector<unsigned int> getRanks(tarch::la::Vector<dim, unsigned int> globalCellIndex) override { return {0}; }

  std::vector<unsigned int> getSourceRanks(tarch::la::Vector<dim, unsigned int> globalCellIndex) override { return {0}; }

  std::vector<unsigned int> getTargetRanks(tarch::la::Vector<dim, unsigned int> globalCellIndex) override { return {0}; }

private:
  const tarch::la::Vector<3, unsigned int> _globalNumberMacroscopicCells;
  const unsigned int _overlap;
  const unsigned int _rank;
  const coupling::IndexConversion<dim>* _indexConversion;
};

// template <unsigned int dim> class coupling::solvers::PreciceInterface : public coupling::interface::MacroscopicSolverInterface<dim> {
// public:
//   PreciceInterface(const tarch::la::Vector<dim, int> globalNumberMacroscopicCells, const unsigned int overlap, const unsigned int rank,
//   const coupling::IndexConversion<dim>* indexConversion)
//       : _globalNumberMacroscopicCells(globalNumberMacroscopicCells), _overlap(overlap), _rank(rank), _indexConversion(indexConversion) {}

//   bool receiveMacroscopicQuantityFromMDSolver(tarch::la::Vector<dim, unsigned int> globalCellIndex) override {
//     bool rcv = true;
//     for (unsigned int currentDim = 0; currentDim < dim; currentDim++) {
//       rcv &= globalCellIndex[currentDim] >= 1 + (_overlap - 1);
//       rcv &= globalCellIndex[currentDim] < _globalNumberMacroscopicCells[currentDim] + 1 - (_overlap - 1);
//     }
//     return rcv;
//   }

//   bool sendMacroscopicQuantityToMDSolver(tarch::la::Vector<3, unsigned int> globalCellIndex) override {
//     bool isGhostCell = false;
//     bool isInner = true;
//     for (unsigned int currentDim = 0; currentDim < dim; currentDim++) {
//       isGhostCell |= globalCellIndex[currentDim] > _globalNumberMacroscopicCells[currentDim];
//       isGhostCell |= globalCellIndex[currentDim] < 1;
//       isInner &= globalCellIndex[currentDim] >= 1 + _overlap;
//       isInner &= globalCellIndex[currentDim] < _globalNumberMacroscopicCells[currentDim] + 1 - _overlap;
//     }
//     return (!isGhostCell) && (!isInner);
//   }

//   std::vector<unsigned int> getRanks(tarch::la::Vector<dim, unsigned int> globalCellIndex) override { return {0}; }

//   std::vector<unsigned int> getSourceRanks(tarch::la::Vector<dim, unsigned int> globalCellIndex) override { return {0}; }

//   std::vector<unsigned int> getTargetRanks(tarch::la::Vector<dim, unsigned int> globalCellIndex) override { return {0}; }

// private:
//   const tarch::la::Vector<3, unsigned int> _globalNumberMacroscopicCells;
//   const unsigned int _overlap;
//   const unsigned int _rank;
//   const coupling::IndexConversion<dim>* _indexConversion;
// };