#pragma once

#include "coupling/CouplingMDDefinitions.h"
#include "coupling/IndexConversion.h"
#include "coupling/datastructures/MacroscopicCell.h"
#include "coupling/indexing/IndexingService.h"
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

namespace precice_scenario{
  template <unsigned int dim> class PreciceAdapter;
}

/**
 * Adapter for the preCICE library
 * 'Vertices' are used in the preCICE context (elements constituting the preCICE mesh)
 * 'Cells' are used in the MaMiCo context (elements consituting the MaMiCo cartesian grid)
 * CellIndexing system should be initialized prior to using this class
 */
template <unsigned int dim> class precice_scenario::PreciceAdapter {
public:
  PreciceAdapter(): _M2mMeshName("mamico-M2m-mesh"), _m2MMeshName("mamico-m2M-mesh"), _M2mVelocityName("VelocityMacro"), _m2MVelocityName("VelocityMicro"), _rank(0) {
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
    if (_M2mCellIndices != NULL) {
      delete[] _M2mCellIndices;
    }
    if (_m2MCellIndices != NULL) {
      delete[] _m2MCellIndices;
    }
    deleteBuffer(_M2mCells);
    deleteBuffer(_m2MCells);
  }

  std::vector<coupling::datastructures::MacroscopicCell<3>*> getM2mCells() {
    return _M2mCells;
  }

  std::vector<coupling::datastructures::MacroscopicCell<3>*> getm2MCells() {
    return _m2MCells;
  }

  unsigned int* getM2mCellIndices() {
    return _M2mCellIndices;
  }

  unsigned int* getm2MCellIndices() {
    return _m2MCellIndices;
  }

  void setMeshes(coupling::interface::MacroscopicSolverInterface<dim>* _macroscopicSolverInterface, const tarch::la::Vector<3, double> mdDomainOffset, const tarch::la::Vector<3, double> macroscopicCellSize) {
    std::vector<double> M2mVertexCoords;
    std::vector<unsigned int> M2mCellIndices;
    using namespace coupling::indexing;
    for (auto cellIndex_v : CellIndex<dim, IndexTrait::vector>()) {
      tarch::la::Vector<3, unsigned int> cellIndex_V = static_cast<tarch::la::Vector<dim, unsigned int>>(cellIndex_v.get());
      if (_macroscopicSolverInterface->sendMacroscopicQuantityToMDSolver(cellIndex_V)) {
        std::vector<unsigned int> ranks = _macroscopicSolverInterface->getSourceRanks(cellIndex_V);
        if (std::find(ranks.begin(), ranks.end(), _rank) != ranks.end()) {
          for (unsigned int currentDim = 0; currentDim < dim; currentDim++) {
            M2mVertexCoords.push_back(mdDomainOffset[currentDim] + cellIndex_v.get()[currentDim] * macroscopicCellSize[currentDim] -
                                    macroscopicCellSize[currentDim] + 0.5 * macroscopicCellSize[currentDim]);
          }
          _M2mCells.push_back(new coupling::datastructures::MacroscopicCell<3>());
          M2mCellIndices.push_back(convertToScalar<dim>(cellIndex_v));
        }
      }
    }
    _M2mVertexNumbers = M2mVertexCoords.size() / dim;
    _M2mVertexCoords = new double[M2mVertexCoords.size()];
    std::copy(M2mVertexCoords.begin(), M2mVertexCoords.end(), _M2mVertexCoords);
    _M2mVertexIndices = new int[_M2mVertexNumbers];
    if (_solverInterface->hasMesh(_M2mMeshName)) {
      _solverInterface->setMeshVertices(_solverInterface->getMeshID(_M2mMeshName), _M2mVertexNumbers, _M2mVertexCoords, _M2mVertexIndices);
    }
    _M2mVertexVelocities = new double[_M2mVertexNumbers * dim];

    _M2mCellIndices = new unsigned int[_M2mCells.size()];
    std::copy(M2mCellIndices.begin(), M2mCellIndices.end(), _M2mCellIndices);
  
    std::vector<double> m2MVertexCoords;
    double offset[dim] = {0, /*MD size*/-120/*CFD size*/-20/*half a mamico cell*/+1.25, 0}; 
    std::vector<unsigned int> m2MCellIndices;
    using namespace coupling::indexing;
    for (auto cellIndex_v : CellIndex<dim, IndexTrait::vector>()) {
      tarch::la::Vector<3, unsigned int> cellIndex_V = static_cast<tarch::la::Vector<dim, unsigned int>>(cellIndex_v.get());
      if (_macroscopicSolverInterface->receiveMacroscopicQuantityFromMDSolver(cellIndex_V)) {
        std::vector<unsigned int> ranks = _macroscopicSolverInterface->getTargetRanks(cellIndex_V);
        if (std::find(ranks.begin(), ranks.end(), _rank) != ranks.end()) {
          for (unsigned int currentDim = 0; currentDim < dim; currentDim++) {
            m2MVertexCoords.push_back(mdDomainOffset[currentDim] + cellIndex_v.get()[currentDim] * macroscopicCellSize[currentDim] -
                                    macroscopicCellSize[currentDim] + 0.5 * macroscopicCellSize[currentDim] + offset[currentDim]);
          }
          _m2MCells.push_back(new coupling::datastructures::MacroscopicCell<3>());
          m2MCellIndices.push_back(convertToScalar<dim>(cellIndex_v));
        }
      }
    }
    _m2MVertexNumbers = m2MVertexCoords.size() / dim;
    _m2MVertexCoords = new double[m2MVertexCoords.size()];
    std::copy(m2MVertexCoords.begin(), m2MVertexCoords.end(), _m2MVertexCoords);
    _m2MVertexIndices = new int[_m2MVertexNumbers];
    if (_solverInterface->hasMesh(_m2MMeshName)) {
      _solverInterface->setMeshVertices(_solverInterface->getMeshID(_m2MMeshName), _m2MVertexNumbers, _m2MVertexCoords, _m2MVertexIndices);
    }
    _m2MVertexVelocities = new double[_m2MVertexNumbers * dim];

    _m2MCellIndices = new unsigned int[_m2MCells.size()];
    std::copy(m2MCellIndices.begin(), m2MCellIndices.end(), _m2MCellIndices);

    if (_solverInterface->hasMesh(_m2MMeshName)) {
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

  double initialize() { return _solverInterface->initialize(); }

  bool isCouplingOngoing() { return _solverInterface->isCouplingOngoing(); }

  double advance(const double dt) { return _solverInterface->advance(dt); }

  void readData(double massCell) {
    if (_solverInterface->hasMesh(_M2mMeshName)) {
      _solverInterface->readBlockVectorData(_solverInterface->getDataID(_M2mVelocityName, _solverInterface->getMeshID(_M2mMeshName)), _M2mVertexNumbers,
                                            _M2mVertexIndices, _M2mVertexVelocities);
      for (size_t i = 0; i < _M2mCells.size(); i++) {
        tarch::la::Vector<3, double> velocity{0.0};
        for (unsigned int currentDim = 0; currentDim < dim; currentDim++) {
          velocity[currentDim] = _M2mVertexVelocities[dim * i + currentDim];
        }
        tarch::la::Vector<3, double> momentum(massCell * velocity);
        _M2mCells[i]->setMicroscopicMass(massCell);
        _M2mCells[i]->setMicroscopicMomentum(momentum);
      }
    }
  }

  void writeData() {
    if (_solverInterface->hasMesh(_m2MMeshName)) {
      for (size_t i = 0; i < _m2MCells.size(); i++) {
        tarch::la::Vector<3, double> velocity{0.0};
        if (_m2MCells[i]->getMacroscopicMass() != 0.0) {
          velocity = (1.0 / _m2MCells[i]->getMacroscopicMass()) * _m2MCells[i]->getMacroscopicMomentum();
        }
        for (unsigned int currentDim = 0; currentDim < dim; currentDim++) {
          _m2MVertexVelocities[dim * i + currentDim] = velocity[currentDim];
        }
      }
      _solverInterface->writeBlockVectorData(_solverInterface->getDataID(_m2MVelocityName, _solverInterface->getMeshID(_m2MMeshName)), _m2MVertexNumbers,
                                             _m2MVertexIndices, _m2MVertexVelocities);
    }
  }

  void writeData(double value) {
    if (_solverInterface->hasMesh(_m2MMeshName)) {
      for (size_t i = 0; i < _m2MCells.size(); i++) {
        _m2MVertexVelocities[dim * i] = 0.0;
        _m2MVertexVelocities[dim * i + 1] = value;
        _m2MVertexVelocities[dim * i + 2] = 0.0;
      }
      _solverInterface->writeBlockVectorData(_solverInterface->getDataID(_m2MVelocityName, _solverInterface->getMeshID(_m2MMeshName)), _m2MVertexNumbers,
                                             _m2MVertexIndices, _m2MVertexVelocities);
    }
  }

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

  void deleteBuffer(std::vector<coupling::datastructures::MacroscopicCell<dim>*>& buffer) const {
    for (unsigned int i = 0; i < buffer.size(); i++) {
      if (buffer[i] != NULL) {
        delete buffer[i];
        buffer[i] = NULL;
      }
    }
    buffer.clear();
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

  std::vector<coupling::datastructures::MacroscopicCell<dim>*> _M2mCells;
  unsigned int* _M2mCellIndices;
  std::vector<coupling::datastructures::MacroscopicCell<dim>*> _m2MCells;
  unsigned int* _m2MCellIndices;
};