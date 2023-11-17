#pragma once

#include "coupling/CouplingMDDefinitions.h"
#include "coupling/IndexConversion.h"
#include "coupling/datastructures/MacroscopicCell.h"
#include "coupling/indexing/IndexingService.h"
#include "coupling/indexing/CellIndex.h"
#include "coupling/scenario/precice/PreciceInterface.h"
#include "precice/precice.hpp"
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
namespace preciceadapter {
  template <unsigned int dim> class PreciceAdapter;
} // namespace preciceadapter
} // namespace coupling

/**
 * Adapter for the preCICE library
 * 'Vertices' are used in the preCICE context (elements constituting the preCICE mesh)
 * 'Cells' are used in the MaMiCo context (elements consituting the MaMiCo cartesian grid)
 * CellIndexing system should be initialized prior to using this class
 */
template <unsigned int dim> class coupling::preciceadapter::PreciceAdapter {
public:
  PreciceAdapter(): _rank(0) {
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
    MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
#endif
    int size = 1;
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
    MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
    _participant = new precice::Participant("mamico", "../precice-config.xml", _rank, size);
  }

  virtual ~PreciceAdapter() {
    if (_participant != NULL) {
      _participant->finalize();
      delete _participant;
      _participant = NULL;
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

  std::vector<coupling::datastructures::MacroscopicCell<3>*> getM2mCells() const {
    return _M2mCells;
  }

  std::vector<coupling::datastructures::MacroscopicCell<3>*> getm2MCells() const {
    return _m2MCells;
  }

  unsigned int* getM2mCellIndices() const {
    return _M2mCellIndices;
  }

  unsigned int* getm2MCellIndices() const {
    return _m2MCellIndices;
  }

  void setMeshes(coupling::preciceadapter::PreciceInterface<dim>* preciceInterface, const tarch::la::Vector<3, double> mdDomainOffset, 
    const tarch::la::Vector<3, double> macroscopicCellSize) {
    std::vector<unsigned int> M2mCellIndices;
    std::vector<unsigned int> m2MCellIndices; 
    using namespace coupling::indexing;
    for (auto cellIndex_v : CellIndex<dim, IndexTrait::vector>()) {
      tarch::la::Vector<dim, unsigned int> cellIndex_V = static_cast<tarch::la::Vector<dim, unsigned int>>(cellIndex_v.get());
      if (preciceInterface->sendMacroscopicQuantityToMDSolver(cellIndex_V)) {
        std::vector<unsigned int> ranks = preciceInterface->getSourceRanks(cellIndex_V);
        if (std::find(ranks.begin(), ranks.end(), _rank) != ranks.end()) {
          addCell(preciceInterface->getMacroscopicToMDSolverMeshName(cellIndex_V), cellIndex_V, _M2mVertexCoordinates, _M2mCells, M2mCellIndices, _M2mVertexToCell, _M2mCellToVertex, mdDomainOffset, macroscopicCellSize, preciceInterface->getMacroscopicToMDSolverMeshOffset(cellIndex_V));
        }
      }
      if (preciceInterface->receiveMacroscopicQuantityFromMDSolver(cellIndex_V)) {
        std::vector<unsigned int> ranks = preciceInterface->getTargetRanks(cellIndex_V);
        if (std::find(ranks.begin(), ranks.end(), _rank) != ranks.end()) {
          addCell(preciceInterface->getMDToMacroscopicSolverMeshName(cellIndex_V), cellIndex_V, _m2MVertexCoordinates, _m2MCells, m2MCellIndices, _m2MVertexToCell, _m2MCellToVertex, mdDomainOffset, macroscopicCellSize, preciceInterface->getMDToMacroscopicSolverMeshOffset(cellIndex_V));
        }
      }
    }
    _M2mCellIndices = new unsigned int[_M2mCells.size()];
    std::copy(M2mCellIndices.begin(), M2mCellIndices.end(), _M2mCellIndices);
    _m2MCellIndices = new unsigned int[_m2MCells.size()];
    std::copy(m2MCellIndices.begin(), m2MCellIndices.end(), _m2MCellIndices);
    initializeVectors(_M2mVertexCoordinates, _M2mVertexIndices, _M2mVertexData, preciceInterface);
    setMeshTetrahedra(_M2mVertexIndices, M2mCellIndices, _M2mVertexToCell, _M2mCellToVertex);
    if (preciceInterface->twoWayCoupling()) {
      initializeVectors(_m2MVertexCoordinates, _m2MVertexIndices, _m2MVertexData, preciceInterface);
      setMeshTetrahedra(_m2MVertexIndices, m2MCellIndices, _m2MVertexToCell, _m2MCellToVertex);
    }
  }

  void initialize() { _participant->initialize(); }

  bool isCouplingOngoing() const { return _participant->isCouplingOngoing(); }

  bool isTimeWindowComplete() const { return _participant->isTimeWindowComplete(); }

  double getMaxTimeStepSize() const { return _participant->getMaxTimeStepSize(); }

  void advance(const double dt) { _participant->advance(dt); }

  void readData(coupling::preciceadapter::PreciceInterface<dim>* preciceInterface) {
    std::map<std::string, std::map<std::string, std::vector<double>>>::iterator it1VertexData;
    for (it1VertexData = _M2mVertexData.begin(); it1VertexData != _M2mVertexData.end(); ++it1VertexData) {
      std::map<std::string, std::vector<double>>::iterator it2VertexData;
      for (it2VertexData = it1VertexData->second.begin(); it2VertexData != it1VertexData->second.end(); ++it2VertexData) {
        _participant->readData(it1VertexData->first, it2VertexData->first, _M2mVertexIndices[it1VertexData->first], 0, it2VertexData->second);
        Data data = preciceInterface->getData(it1VertexData->first, it2VertexData->first);
        std::vector<int> vertexIndices = _M2mVertexIndices[it1VertexData->first];
        for (size_t i = 0; i < vertexIndices.size(); ++i) {
          coupling::datastructures::MacroscopicCell<dim>* cell = _M2mCells[_M2mVertexToCell[it1VertexData->first][i]];
          switch (data.type) {
            case DataType::scalar:
              preciceInterface->readScalarData(it1VertexData->first, it2VertexData->first, cell, it2VertexData->second[i]);
              break;
            case DataType::vector:
              preciceInterface->readVectorData(it1VertexData->first, it2VertexData->first, cell, it2VertexData->second[i*dim], it2VertexData->second[i*dim+1], it2VertexData->second[i*dim+2]);
              break; 
          }
        }
      }
    }
  }

  void writeData(coupling::preciceadapter::PreciceInterface<dim>* preciceInterface) {
    std::map<std::string, std::map<std::string, std::vector<double>>>::iterator it1VertexData;
    for (it1VertexData = _m2MVertexData.begin(); it1VertexData != _m2MVertexData.end(); ++it1VertexData) {
      std::map<std::string, std::vector<double>>::iterator it2VertexData;
      for (it2VertexData = it1VertexData->second.begin(); it2VertexData != it1VertexData->second.end(); ++it2VertexData) {
        Data data = preciceInterface->getData(it1VertexData->first, it2VertexData->first);
        std::vector<int> vertexIndices = _m2MVertexIndices[it1VertexData->first];
        for (size_t i = 0; i < vertexIndices.size(); ++i) {
          coupling::datastructures::MacroscopicCell<dim>* cell = _m2MCells[_m2MVertexToCell[it1VertexData->first][i]];
          switch (data.type) {
            case DataType::scalar:
              preciceInterface->writeScalarData(it1VertexData->first, it2VertexData->first, cell, it2VertexData->second[i]);
              break;
            case DataType::vector:
              preciceInterface->writeVectorData(it1VertexData->first, it2VertexData->first, cell, it2VertexData->second[i*dim], it2VertexData->second[i*dim+1], it2VertexData->second[i*dim+2]);
              break; 
          }
        }
        _participant->writeData(it1VertexData->first, it2VertexData->first, _m2MVertexIndices[it1VertexData->first], it2VertexData->second);
      }
    }
  }

private:
  /**
   * add a vertex and a cell in vertexCoordinates map and cells map for this mesh name 
   * and compute the vertex/cell mapping
   */ 
  void addCell(const std::string& meshName, tarch::la::Vector<dim, unsigned int>& cellIndex, std::map<std::string, std::vector<double>>& vertexCoordinates,
    std::vector<coupling::datastructures::MacroscopicCell<dim>*>& cells, std::vector<unsigned int>& cellIndices,
    std::map<std::string, std::map<int, unsigned int>>& vertexToCell, std::map<std::string, std::map<unsigned int, int>>& cellToVertex,
    const tarch::la::Vector<dim, double>& mdDomainOffset, const tarch::la::Vector<3, double>& macroscopicCellSize, 
    const tarch::la::Vector<dim, double> offset) {
    for (unsigned int currentDim = 0; currentDim < dim; currentDim++) {
      vertexCoordinates[meshName].push_back(mdDomainOffset[currentDim] + cellIndex[currentDim] * macroscopicCellSize[currentDim] -
                              macroscopicCellSize[currentDim] + 0.5 * macroscopicCellSize[currentDim] + offset[currentDim]);
    }
    cells.push_back(new coupling::datastructures::MacroscopicCell<dim>());
    coupling::indexing::CellIndex<dim, coupling::indexing::IndexTrait::vector> cellIndex_v{static_cast<tarch::la::Vector<dim,int>>(cellIndex)};
    cellIndices.push_back(coupling::indexing::convertToScalar<dim>(cellIndex_v));
    vertexToCell[meshName][vertexCoordinates[meshName].size()/dim - 1] = cells.size()-1;
    cellToVertex[meshName][cells.size()-1] = vertexCoordinates[meshName].size()/dim-1;
  }

  /**
   * initialize the vertex indices and vertex data maps based on the vertex coordinates map
   * set the participant precice mesh vertices 
   */
  void initializeVectors(std::map<std::string, std::vector<double>>& vertexCoordinates, 
    std::map<std::string, std::vector<int>>& vertexIndices, std::map<std::string, std::map<std::string, std::vector<double>>>& vertexData,
    coupling::preciceadapter::PreciceInterface<dim>* preciceInterface) {
    std::map<std::string, std::vector<double>>::iterator itVertexCoordinates;
    for (itVertexCoordinates = vertexCoordinates.begin(); itVertexCoordinates != vertexCoordinates.end(); ++itVertexCoordinates) {
      size_t numberOfCoordinates = itVertexCoordinates->second.size();
      vertexIndices[itVertexCoordinates->first] = std::vector<int>(numberOfCoordinates/dim);
      _participant->setMeshVertices(itVertexCoordinates->first, itVertexCoordinates->second, vertexIndices[itVertexCoordinates->first]);
      for (const Data& data : preciceInterface->getData(itVertexCoordinates->first)) {
        size_t dataSize = numberOfCoordinates;
        if (data.type == DataType::scalar) dataSize/=dim;
        vertexData[itVertexCoordinates->first][data.name] = std::vector<double>(dataSize);
      }
    }
  }

  /**
   * set the participant mesh tetrahedra
   */
  void setMeshTetrahedra(const std::map<std::string, std::vector<int>>& vertexIndices, const std::vector<unsigned int>& cellIndices, 
    const std::map<std::string, std::map<int, unsigned int>>& vertexToCell, const std::map<std::string, std::map<unsigned int, int>>& cellToVertex) {
    std::map<std::string, std::vector<int>>::const_iterator itVertexIndices;
    for (itVertexIndices = vertexIndices.begin(); itVertexIndices != vertexIndices.end(); ++itVertexIndices) {
      for (size_t i = 0; i < itVertexIndices->second.size(); ++i) {
        using CellIndex_v = coupling::indexing::CellIndex<dim, coupling::indexing::IndexTrait::vector>;
        using CellIndex_s = coupling::indexing::CellIndex<dim>;
        CellIndex_v cellIndex_v = CellIndex_s{cellIndices[vertexToCell.at(itVertexIndices->first).at(i)]};
        const int numberOfNeighbors = 7;
        tarch::la::Vector<3, int> directions[numberOfNeighbors] = {{1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}, {0, 0, 1}};
        int neighborVertexIndices[numberOfNeighbors];
        for (size_t j = 0; j < numberOfNeighbors; j++) {
          CellIndex_s neighborCellIndex_s = CellIndex_v{cellIndex_v.get() + directions[j]};
          std::vector<unsigned int>::const_iterator itCells = std::find(cellIndices.begin(), cellIndices.end(), neighborCellIndex_s.get());
          if (itCells != cellIndices.end()) {
            neighborVertexIndices[j] = itVertexIndices->second[cellToVertex.at(itVertexIndices->first).at(itCells - cellIndices.begin())];
          } else {
            neighborVertexIndices[j] = -1;
          }
        }
        int vertexIndex = itVertexIndices->second[i];
        if (neighborVertexIndices[0] > -1 && neighborVertexIndices[5] > -1 && neighborVertexIndices[6] > -1)
          _participant->setMeshTetrahedron(itVertexIndices->first, vertexIndex, neighborVertexIndices[0], neighborVertexIndices[5], neighborVertexIndices[6]);
        if (neighborVertexIndices[0] > -1 && neighborVertexIndices[2] > -1 && neighborVertexIndices[5] > -1)
          _participant->setMeshTetrahedron(itVertexIndices->first, vertexIndex, neighborVertexIndices[0], neighborVertexIndices[2], neighborVertexIndices[5]);
        if (neighborVertexIndices[0] > -1 && neighborVertexIndices[5] > -1 && neighborVertexIndices[6] > -1 && neighborVertexIndices[3] > -1)
          _participant->setMeshTetrahedron(itVertexIndices->first, neighborVertexIndices[0], neighborVertexIndices[5], neighborVertexIndices[6], neighborVertexIndices[3]);
        if (neighborVertexIndices[0] > -1 && neighborVertexIndices[1] > -1 && neighborVertexIndices[2] > -1 && neighborVertexIndices[5] > -1)
          _participant->setMeshTetrahedron(itVertexIndices->first, neighborVertexIndices[0], neighborVertexIndices[1], neighborVertexIndices[2], neighborVertexIndices[5]);
        if (neighborVertexIndices[0] > -1 && neighborVertexIndices[1] > -1 && neighborVertexIndices[4] > -1 && neighborVertexIndices[5] > -1)
          _participant->setMeshTetrahedron(itVertexIndices->first, neighborVertexIndices[0], neighborVertexIndices[1], neighborVertexIndices[4], neighborVertexIndices[5]);
        if (neighborVertexIndices[0] > -1 && neighborVertexIndices[3] > -1 && neighborVertexIndices[4] > -1 && neighborVertexIndices[5] > -1)
          _participant->setMeshTetrahedron(itVertexIndices->first, neighborVertexIndices[0], neighborVertexIndices[3], neighborVertexIndices[4], neighborVertexIndices[5]);
      }
    }
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

  // rank of this adapter
  int _rank;
  // pointer to the preCICE participant
  precice::Participant* _participant = nullptr;
  // macro to micro preCICE data containers
  std::map<std::string, std::vector<int>> _M2mVertexIndices;
  std::map<std::string, std::vector<double>> _M2mVertexCoordinates;
  std::map<std::string, std::map<std::string, std::vector<double>>> _M2mVertexData;
  // micro to macro preCICE data containers
  std::map<std::string, std::vector<int>> _m2MVertexIndices;
  std::map<std::string, std::vector<double>> _m2MVertexCoordinates;
  std::map<std::string, std::map<std::string, std::vector<double>>> _m2MVertexData;
  // macro to micro mapping between vertex arrays and cell array
  std::map<std::string, std::map<int, unsigned int>> _M2mVertexToCell;
  std::map<std::string, std::map<unsigned int, int>> _M2mCellToVertex;
  // micro to macro mapping between vertex arrays and cell array
  std::map<std::string, std::map<int, unsigned int>> _m2MVertexToCell;
  std::map<std::string, std::map<unsigned int, int>> _m2MCellToVertex;
  // macro to micro MaMiCo data containers
  std::vector<coupling::datastructures::MacroscopicCell<dim>*> _M2mCells;
  unsigned int* _M2mCellIndices;
  // micro to macro MaMiCo data containers
  std::vector<coupling::datastructures::MacroscopicCell<dim>*> _m2MCells;
  unsigned int* _m2MCellIndices;
};