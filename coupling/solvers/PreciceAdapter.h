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
    if (_coordsM2mCells != NULL) {
      delete[] _coordsM2mCells;
    }
    if (_velocityM2mCells != NULL) {
      delete[] _velocityM2mCells;
    }
    if (_vertexM2mCellIDs != NULL) {
      delete[] _vertexM2mCellIDs;
    }
    if (_coordsm2MCells != NULL) {
      delete[] _coordsm2MCells;
    }
    if (_velocitym2MCells != NULL) {
      delete[] _velocitym2MCells;
    }
    if (_vertexm2MCellIDs != NULL) {
      delete[] _vertexm2MCellIDs;
    }
  }

  void setMeshes(const unsigned int* const M2mCellGlobalIndices, size_t numberOfM2mCells, const unsigned int* const m2MCellGlobalIndices,
                 size_t numberOfm2MCells, const tarch::la::Vector<3, double> mdDomainOffset, const tarch::la::Vector<3, double> macroscopicCellSize) {
    if (_solverInterface->hasMesh(_M2mMeshName)) {
      std::vector<double> coordsM2mCells;
      _numberOfM2mCells = numberOfM2mCells;
      for (size_t i = 0; i < _numberOfM2mCells; i++) {
        const unsigned int M2mCellGlobalIndex = M2mCellGlobalIndices[i];
        tarch::la::Vector<dim, int> M2mCellGlobalVectorIndex = coupling::indexing::convertToVector<dim>({M2mCellGlobalIndex});
        for (unsigned int currentDim = 0; currentDim < dim; currentDim++) {
          coordsM2mCells.push_back(mdDomainOffset[currentDim] + M2mCellGlobalVectorIndex[currentDim] * macroscopicCellSize[currentDim] -
                                   macroscopicCellSize[currentDim] + 0.5 * macroscopicCellSize[currentDim]);
        }
      }
      _coordsM2mCells = new double[coordsM2mCells.size()];
      std::copy(coordsM2mCells.begin(), coordsM2mCells.end(), _coordsM2mCells);
      _vertexM2mCellIDs = new int[_numberOfM2mCells];
      _solverInterface->setMeshVertices(_solverInterface->getMeshID(_M2mMeshName), _numberOfM2mCells, _coordsM2mCells, _vertexM2mCellIDs);
      _velocityM2mCells = new double[_numberOfM2mCells * dim];
    }
    if (_solverInterface->hasMesh(_m2MMeshName)) {
      std::vector<double> coordsm2MCells;
      _numberOfm2MCells = numberOfm2MCells;
      for (size_t i = 0; i < _numberOfm2MCells; i++) {
        const unsigned int m2MCellGlobalIndex = m2MCellGlobalIndices[i];
        tarch::la::Vector<dim, int> m2MCellGlobalVectorIndex = coupling::indexing::convertToVector<dim>({m2MCellGlobalIndex});
        for (unsigned int currentDim = 0; currentDim < dim; currentDim++) {
          coordsm2MCells.push_back(mdDomainOffset[currentDim] + m2MCellGlobalVectorIndex[currentDim] * macroscopicCellSize[currentDim] -
                                   macroscopicCellSize[currentDim] + 0.5 * macroscopicCellSize[currentDim]);
        }
      }
      _coordsm2MCells = new double[coordsm2MCells.size()];
      std::copy(coordsm2MCells.begin(), coordsm2MCells.end(), _coordsm2MCells);
      _vertexm2MCellIDs = new int[_numberOfm2MCells];
      _solverInterface->setMeshVertices(_solverInterface->getMeshID(_m2MMeshName), _numberOfm2MCells, _coordsm2MCells, _vertexm2MCellIDs);
      _velocitym2MCells = new double[_numberOfm2MCells * dim];
      for (size_t i = 0; i < _numberOfm2MCells * dim; i++)
        _velocitym2MCells[i] = 0.0;

      for (size_t i = 0; i < _numberOfm2MCells; i++) {
        const unsigned int m2MCellGlobalIndex = m2MCellGlobalIndices[i];
        tarch::la::Vector<dim, int> m2MCellGlobalVectorIndex = coupling::indexing::convertToVector<dim>({m2MCellGlobalIndex});
        std::vector<int> tetrahedronsNodeIDs = getTetrahedronsNodeIDs(i, m2MCellGlobalVectorIndex, m2MCellGlobalIndices, numberOfm2MCells);
        for (size_t tetrahedronIndex = 0; tetrahedronIndex < tetrahedronsNodeIDs.size() / 4; tetrahedronIndex++) {
          _solverInterface->setMeshTetrahedron(_solverInterface->getMeshID(_m2MMeshName), tetrahedronsNodeIDs[tetrahedronIndex * 4],
                                               tetrahedronsNodeIDs[tetrahedronIndex * 4 + 1], tetrahedronsNodeIDs[tetrahedronIndex * 4 + 2],
                                               tetrahedronsNodeIDs[tetrahedronIndex * 4 + 3]);
        }
      }
    }
  }

  double initialize() { return _solverInterface->initialize(); }

  bool isCouplingOngoing() { return _solverInterface->isCouplingOngoing(); }

  double advance(const double dt) { return _solverInterface->advance(dt); }

  void readData() {
    if (_solverInterface->hasMesh(_M2mMeshName)) {
      _solverInterface->readBlockVectorData(_solverInterface->getDataID(_M2mVelocityName, _solverInterface->getMeshID(_M2mMeshName)), _numberOfM2mCells,
                                            _vertexM2mCellIDs, _velocityM2mCells);
    }
  }

  void writeData() {
    if (_solverInterface->hasMesh(_m2MMeshName)) {
      _solverInterface->writeBlockVectorData(_solverInterface->getDataID(_m2MVelocityName, _solverInterface->getMeshID(_m2MMeshName)), _numberOfm2MCells,
                                             _vertexm2MCellIDs, _velocitym2MCells);
    }
  }

  // Couette solver methods, to be moved.

  tarch::la::Vector<3, double> getVelocity(tarch::la::Vector<3, double> pos) const {
    tarch::la::Vector<3, double> vel(0.0);
    if (_solverInterface->hasMesh(_M2mMeshName)) {
      unsigned int cellIndex = 0;
      while (cellIndex < _numberOfM2mCells &&
             !(_coordsM2mCells[dim * cellIndex] == pos[0] && _coordsM2mCells[dim * cellIndex + 1] == pos[1] && _coordsM2mCells[dim * cellIndex + 2] == pos[2]))
        cellIndex++;
      if (cellIndex < _numberOfM2mCells) {
        for (unsigned int currentDim = 0; currentDim < dim; currentDim++) {
          vel[currentDim] = _velocityM2mCells[dim * cellIndex + currentDim];
        }
      }
    }
    return vel;
  }

  void setMDBoundaryValues(std::vector<coupling::datastructures::MacroscopicCell<3>*>& m2MBuffer, const unsigned int* const m2MCellGlobalIndices) {
    if (_solverInterface->hasMesh(_m2MMeshName)) {
      for (size_t i = 0; i < _numberOfm2MCells; i++) {
        if (m2MBuffer[i]->getMacroscopicMass() != 0.0) {
          tarch::la::Vector<3, double> vel((1.0 / m2MBuffer[i]->getMacroscopicMass()) * m2MBuffer[i]->getMacroscopicMomentum());
          for (unsigned int currentDim = 0; currentDim < dim; currentDim++)
            _velocitym2MCells[dim * i + currentDim] = vel[currentDim];
        }
      }
    }
  }

private:
  std::vector<int> getTetrahedronsNodeIDs(const int nodePreciceID, const tarch::la::Vector<3, int> nodeMamicoVIndex,
                                          const unsigned int* const m2MCellGlobalIndices, size_t numberOfm2MCells) {
    std::vector<int> tetrahedronsNodeIDs;
    const int numberOfNeighbors = 7;
    tarch::la::Vector<3, int> directionNeighbors[numberOfNeighbors] = {{1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}, {0, 0, 1}};
    int neighborVertexIDSs[numberOfNeighbors];
    for (size_t iDirectionNeighbors = 0; iDirectionNeighbors < numberOfNeighbors; iDirectionNeighbors++) {
      tarch::la::Vector<3, int> neighborMamicoVIndex = nodeMamicoVIndex + directionNeighbors[iDirectionNeighbors];
      using coupling::indexing::BaseIndex;
      using coupling::indexing::convertToScalar;
      using coupling::indexing::IndexTrait;
      int neighborMamicoIndex = convertToScalar<dim>(BaseIndex<dim>(neighborMamicoVIndex));
      int neighborPreciceID = 0;
      while (neighborPreciceID < (int)numberOfm2MCells && (int)m2MCellGlobalIndices[neighborPreciceID] != neighborMamicoIndex)
        neighborPreciceID++;
      if (neighborPreciceID < (int)numberOfm2MCells) {
        neighborVertexIDSs[iDirectionNeighbors] = neighborPreciceID;
      } else {
        neighborVertexIDSs[iDirectionNeighbors] = -1;
      }
    }
    if (neighborVertexIDSs[0] > -1 && neighborVertexIDSs[5] > -1 && neighborVertexIDSs[6] > -1)
      tetrahedronsNodeIDs.insert(tetrahedronsNodeIDs.end(), {nodePreciceID, neighborVertexIDSs[0], neighborVertexIDSs[5], neighborVertexIDSs[6]});
    if (neighborVertexIDSs[0] > -1 && neighborVertexIDSs[2] > -1 && neighborVertexIDSs[5] > -1)
      tetrahedronsNodeIDs.insert(tetrahedronsNodeIDs.end(), {nodePreciceID, neighborVertexIDSs[0], neighborVertexIDSs[2], neighborVertexIDSs[5]});
    if (neighborVertexIDSs[0] > -1 && neighborVertexIDSs[5] > -1 && neighborVertexIDSs[6] > -1 && neighborVertexIDSs[3] > -1)
      tetrahedronsNodeIDs.insert(tetrahedronsNodeIDs.end(), {neighborVertexIDSs[0], neighborVertexIDSs[5], neighborVertexIDSs[6], neighborVertexIDSs[3]});
    if (neighborVertexIDSs[0] > -1 && neighborVertexIDSs[1] > -1 && neighborVertexIDSs[2] > -1 && neighborVertexIDSs[5] > -1)
      tetrahedronsNodeIDs.insert(tetrahedronsNodeIDs.end(), {neighborVertexIDSs[0], neighborVertexIDSs[1], neighborVertexIDSs[2], neighborVertexIDSs[5]});
    if (neighborVertexIDSs[0] > -1 && neighborVertexIDSs[1] > -1 && neighborVertexIDSs[4] > -1 && neighborVertexIDSs[5] > -1)
      tetrahedronsNodeIDs.insert(tetrahedronsNodeIDs.end(), {neighborVertexIDSs[0], neighborVertexIDSs[1], neighborVertexIDSs[4], neighborVertexIDSs[5]});
    if (neighborVertexIDSs[0] > -1 && neighborVertexIDSs[3] > -1 && neighborVertexIDSs[4] > -1 && neighborVertexIDSs[5] > -1)
      tetrahedronsNodeIDs.insert(tetrahedronsNodeIDs.end(), {neighborVertexIDSs[0], neighborVertexIDSs[3], neighborVertexIDSs[4], neighborVertexIDSs[5]});
    return tetrahedronsNodeIDs;
  }
  const std::string _M2mMeshName;
  const std::string _m2MMeshName;
  const std::string _M2mVelocityName;
  const std::string _m2MVelocityName;

  int _rank;

  precice::SolverInterface* _solverInterface = nullptr;

  int* _vertexM2mCellIDs = nullptr;
  double* _coordsM2mCells = nullptr;
  unsigned int _numberOfM2mCells = 0;
  double* _velocityM2mCells = nullptr;

  int* _vertexm2MCellIDs = nullptr;
  double* _coordsm2MCells = nullptr;
  unsigned int _numberOfm2MCells = 0;
  double* _velocitym2MCells = nullptr;
};

template <unsigned int dim> class coupling::solvers::PreciceInterface : public coupling::interface::MacroscopicSolverInterface<dim> {
public:
  PreciceInterface(const tarch::la::Vector<3, int> globalNumberMacroscopicCells, const unsigned int overlap, const unsigned int rank)
      : _globalNumberMacroscopicCells(globalNumberMacroscopicCells), _overlap(overlap), _rank(rank) {}

  bool receiveMacroscopicQuantityFromMDSolver(tarch::la::Vector<dim, unsigned int> globalCellIndex) override {
    bool rcv = true;
    for (unsigned int currentDim = 0; currentDim < dim; currentDim++) {
      rcv &= globalCellIndex[currentDim] >= 1 + (_overlap - 1);
      rcv &= globalCellIndex[currentDim] < _globalNumberMacroscopicCells[currentDim] + 1 - (_overlap - 1);
    }
    return rcv;
  }

  bool sendMacroscopicQuantityToMDSolver(tarch::la::Vector<3, unsigned int> globalCellIndex) override {
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
};