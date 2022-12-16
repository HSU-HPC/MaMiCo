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

using namespace coupling::indexing;

template <unsigned int dim> class PreciceAdapter {
public:
  PreciceAdapter(const tarch::la::Vector<3, double> mdDomainOffset, const tarch::la::Vector<3, double> macroscopicCellSize)
      : _rank(0), _mdDomainOffset(mdDomainOffset), _macroscopicCellSize(macroscopicCellSize) {
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
    if (_coordsm2MCells != NULL) {
      delete[] _coordsm2MCells;
    }
    if (_velocitym2MCells != NULL) {
      delete[] _velocitym2MCells;
    }
  }

  void setMeshes(const unsigned int* const M2mCellGlobalIndices, size_t numberOfM2mCells, const unsigned int* const m2MCellGlobalIndices,
                 size_t numberOfm2MCells) {
    if (_solverInterface->hasMesh(M2m_MESH_NAME)) {
      std::vector<double> coordsM2mCells;
      _numberOfM2mCells = numberOfM2mCells;
      for (size_t i = 0; i < _numberOfM2mCells; i++) {
        const unsigned int M2mCellGlobalIndex = M2mCellGlobalIndices[i];
        tarch::la::Vector<dim, int> M2mCellGlobalVectorIndex = coupling::indexing::convertToVector<dim>({M2mCellGlobalIndex});
        for (unsigned int currentDim = 0; currentDim < dim; currentDim++) {
          coordsM2mCells.push_back(_mdDomainOffset[currentDim] + M2mCellGlobalVectorIndex[currentDim] * _macroscopicCellSize[currentDim] -
                                   _macroscopicCellSize[currentDim] + 0.5 * _macroscopicCellSize[currentDim]);
        }
      }
      _coordsM2mCells = new double[coordsM2mCells.size()];
      std::copy(coordsM2mCells.begin(), coordsM2mCells.end(), _coordsM2mCells);
      _vertexM2mCellIDs = new int[_numberOfM2mCells];
      _solverInterface->setMeshVertices(_solverInterface->getMeshID(M2m_MESH_NAME), _numberOfM2mCells, _coordsM2mCells, _vertexM2mCellIDs);
      _velocityM2mCells = new double[_numberOfM2mCells * dim];
    }
    if (_solverInterface->hasMesh(m2M_MESH_NAME)) {
      std::vector<double> coordsm2MCells;
      _numberOfm2MCells = numberOfm2MCells;
      for (size_t i = 0; i < _numberOfm2MCells; i++) {
        const unsigned int m2MCellGlobalIndex = m2MCellGlobalIndices[i];
        tarch::la::Vector<dim, int> m2MCellGlobalVectorIndex = coupling::indexing::convertToVector<dim>({m2MCellGlobalIndex});
        for (unsigned int currentDim = 0; currentDim < dim; currentDim++) {
          coordsm2MCells.push_back(_mdDomainOffset[currentDim] + m2MCellGlobalVectorIndex[currentDim] * _macroscopicCellSize[currentDim] -
                                   _macroscopicCellSize[currentDim] + 0.5 * _macroscopicCellSize[currentDim]);
        }
      }
      _coordsm2MCells = new double[coordsm2MCells.size()];
      std::copy(coordsm2MCells.begin(), coordsm2MCells.end(), _coordsm2MCells);
      _vertexm2MCellIDs = new int[_numberOfm2MCells];
      _solverInterface->setMeshVertices(_solverInterface->getMeshID(m2M_MESH_NAME), _numberOfm2MCells, _coordsm2MCells, _vertexm2MCellIDs);
      _velocitym2MCells = new double[_numberOfm2MCells * dim];
      for (size_t i = 0; i < _numberOfm2MCells * dim; i++)
        _velocitym2MCells[i] = 0.0;

      for (size_t i = 0; i < _numberOfm2MCells; i++) {
        const unsigned int m2MCellGlobalIndex = m2MCellGlobalIndices[i];
        tarch::la::Vector<dim, int> m2MCellGlobalVectorIndex = coupling::indexing::convertToVector<dim>({m2MCellGlobalIndex});
        std::vector<int> tetrahedronsNodeIDs = getTetrahedronsNodeIDs(i, m2MCellGlobalVectorIndex, m2MCellGlobalIndices, numberOfm2MCells);
        for (size_t tetrahedronIndex = 0; tetrahedronIndex < tetrahedronsNodeIDs.size() / 4; tetrahedronIndex++) {
          _solverInterface->setMeshTetrahedron(_solverInterface->getMeshID(m2M_MESH_NAME), tetrahedronsNodeIDs[tetrahedronIndex * 4],
                                               tetrahedronsNodeIDs[tetrahedronIndex * 4 + 1], tetrahedronsNodeIDs[tetrahedronIndex * 4 + 2],
                                               tetrahedronsNodeIDs[tetrahedronIndex * 4 + 3]);
        }
      }
    }
  }

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

  double initialize() {
    return _solverInterface->initialize();
  }

  bool isCouplingOngoing() { return _solverInterface->isCouplingOngoing(); }

  double advance(const double dt) { return _solverInterface->advance(dt); }

  void readData() {
    if (_solverInterface->hasMesh(M2m_MESH_NAME)) {
      _solverInterface->readBlockVectorData(_solverInterface->getDataID(M2m_VELOCITY_NAME, _solverInterface->getMeshID(M2m_MESH_NAME)), _numberOfM2mCells,
                                            _vertexM2mCellIDs, _velocityM2mCells);
    }
  }

  tarch::la::Vector<3, double> getVelocity(tarch::la::Vector<3, double> pos) const {
    tarch::la::Vector<3, double> vel(0.0);
    if (_solverInterface->hasMesh(M2m_MESH_NAME)) {
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

  void writeData() {
    if (_solverInterface->hasMesh(m2M_MESH_NAME)) {
      _solverInterface->writeBlockVectorData(_solverInterface->getDataID(m2M_VELOCITY_NAME, _solverInterface->getMeshID(m2M_MESH_NAME)), _numberOfm2MCells,
                                             _vertexm2MCellIDs, _velocitym2MCells);
    }
  }

  void setMDBoundaryValues(std::vector<coupling::datastructures::MacroscopicCell<3>*>& m2MBuffer, const unsigned int* const m2MCellGlobalIndices) {
    if (_solverInterface->hasMesh(m2M_MESH_NAME)) {
      for (size_t i = 0; i < _numberOfm2MCells; i++) {
        tarch::la::Vector<3, double> vel((1.0 / m2MBuffer[i]->getMacroscopicMass()) * m2MBuffer[i]->getMacroscopicMomentum());
        for (unsigned int currentDim = 0; currentDim < dim; currentDim++)
          _velocitym2MCells[dim * i + currentDim] = vel[currentDim];
      }
    }
  }

private:
  const static std::string M2m_MESH_NAME;
  const static std::string M2m_VELOCITY_NAME;
  const static std::string m2M_MESH_NAME;
  const static std::string m2M_VELOCITY_NAME;

  int _rank;
  const tarch::la::Vector<3, double> _mdDomainOffset;
  const tarch::la::Vector<3, double> _macroscopicCellSize;

  precice::SolverInterface* _solverInterface = NULL;

  int* _vertexM2mCellIDs;
  double* _coordsM2mCells;
  unsigned int _numberOfM2mCells;
  double* _velocityM2mCells;

  int* _vertexm2MCellIDs;
  double* _coordsm2MCells;
  unsigned int _numberOfm2MCells;
  double* _velocitym2MCells;
};

template<> const std::string PreciceAdapter<3>::M2m_MESH_NAME = "mamico-M2m-mesh";
template<> const std::string PreciceAdapter<3>::m2M_MESH_NAME = "mamico-m2M-mesh";
template<> const std::string PreciceAdapter<3>::M2m_VELOCITY_NAME = "VelocityMacro";
template<> const std::string PreciceAdapter<3>::m2M_VELOCITY_NAME = "VelocityMicro";