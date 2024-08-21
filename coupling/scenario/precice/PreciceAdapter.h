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
namespace precice {
  template <unsigned int dim> class PreciceAdapter;
} // namespace preciceadapter
} // namespace coupling

/**
 * MaMiCo adapter for the preCICE library
 * 'Vertices' are used in the preCICE context (elements constituting the preCICE mesh)
 * 'Cells' are used in the MaMiCo context (elements consituting the MaMiCo cartesian grid)
 * Indexing system should be initialized prior to using this class
 */
template <unsigned int dim> class coupling::precice::PreciceAdapter {
public:

  /**
   * Create a new preCICE participant in the coupling scheme and read the preCICE config file
   */
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

  /**
   * Call finalize on the participant and delete it
   */
  virtual ~PreciceAdapter() {
    if (_participant != NULL) {
      _participant->finalize();
      delete _participant;
      _participant = NULL;
    }
  }

  void setMeshes(coupling::preciceadapter::PreciceInterface<dim>* preciceInterface) {
    using namespace coupling::indexing;
    for (auto idx : I08()) {
      if (!I12::contains(idx)) {
        if (tarch::utils::contains(preciceInterface->getSourceRanks(idx), (unsigned int)_rank)) {
          addCell(idx, _Macro2MDCouplingCells, 
            preciceInterface->getMacro2MDSolverMeshName(idx), _Macro2MDMeshes, 
            _Macro2MDVertexToCouplingCell, _Macro2MDCouplingCellToVertex, 
            preciceInterface->getMacro2MDSolverMeshOffset(idx));
        }
      }
    }
    for (auto idx : I12()) {
      if (tarch::utils::contains(preciceInterface->getTargetRanks(idx), (unsigned int)_rank)) {
        addCell(idx, _MD2MacroCouplingCells, 
          preciceInterface->getMD2MacroSolverMeshName(idx), _MD2MacroMeshes, 
          _MD2MacroVertexToCouplingCell, _MD2MacroCouplingCell,
          preciceInterface->getMD2MacroSolverMeshOffset(idx));
      }
    }
    setMeshesVertices(_Macro2MDMeshes, _Macro2MDIndices, _Macro2MDData, preciceInterface);
    setMeshTetrahedra(_M2mVertexIndices, M2mCellIndices, _M2mVertexToCell, _M2mCellToVertex);
    if (preciceInterface->twoWayCoupling()) {
      initializeVectors(_MD2MacroMeshes, _MD2MacroIndices, _MD2MacroData, preciceInterface);
      setMeshTetrahedra(_m2MVertexIndices, m2MCellIndices, _m2MVertexToCell, _m2MCellToVertex);
    }
  }

  void initialize() { _participant->initialize(); }

  bool isCouplingOngoing() const { return _participant->isCouplingOngoing(); }

  bool isTimeWindowComplete() const { return _participant->isTimeWindowComplete(); }

  double getMaxTimeStepSize() const { return _participant->getMaxTimeStepSize(); }

  bool requiresWritingCheckpoint() const { return _participant->requiresWritingCheckpoint(); }

  bool requiresReadingCheckpoint() const { return _participant->requiresReadingCheckpoint(); }

  void advance(const double dt) { _participant->advance(dt); }

  void readData() {
    for (auto const& [meshName, data] : _Macro2MDData) {
      for (auto const& [dataName, dataValues] : data) {
        _participant->readData(meshName, dataName, meshIndices, 0, dataValues);
      }
    }
  }

  template <class MY_LINKEDCELL, unsigned int dim>
  void sendFromMacro2MD(coupling::preciceadapter::PreciceInterface<dim>* preciceInterface, coupling::services::MultiMDCellService<MY_LINKEDCELL, dim>* multiMDCellService) {
    for (auto pair : _Macro2MDCouplingCells) { // loop over all the macro to micro coupling cells
      I01 idx:
      coupling::datastructures::CouplingCell<3>* couplingCell;
      std::tie(couplingCell, idx) = pair;
      // this cell belongs to one macro to micro mesh
      std::string meshName;
      unsigned int vertexIndex;
      std::tie(meshName, vertexIndex) = couplingCellMapping[I00{idx}];
      for (auto const& [dataName, dataValues] : _Macro2MDData[meshName]) {
        auto dataDescription = preciceInterface->getDataDescription(meshName, dataName);
        switch (dataDescription.type) {
          case DataType::scalar:
            preciceInterface->readScalarData(meshName, dataName, couplingCell, dataValues[vertexIndex]);
            break;
          case DataType::vector:
            preciceInterface->readVectorData(meshName, dataName, couplingCell, dataValues[vertexIndex*dim], dataValues[vertexIndex*dim+1], dataValues[vertexIndex*dim+2]);
            break; 
        }
      }
    }
    multiMDCellService->sendFromMacro2MD(_Macro2MDCouplingCells);
  }

  template <class MY_LINKEDCELL, unsigned int dim>
  void sendFromMD2Macro(coupling::preciceadapter::PreciceInterface<dim>* preciceInterface, coupling::services::MultiMDCellService<MY_LINKEDCELL, dim>* multiMDCellService) {
    multiMDCellService->sendFromMacro2MD(_MD2MacroCouplingCells);
    for (auto pair : _MD2MacroCouplingCells) { // loop over all the macro to micro coupling cells
      I01 idx:
      coupling::datastructures::CouplingCell<3>* couplingCell;
      std::tie(couplingCell, idx) = pair;
      // this cell belongs to one macro to micro mesh
      std::string meshName;
      unsigned int vertexIndex;
      std::tie(meshName, vertexIndex) = couplingCellMapping[I00{idx}];
      for (auto const& [dataName, dataValues] : _MD2MacroData[meshName]) {
        auto dataDescription = preciceInterface->getDataDescription(meshName, dataName);
        switch (dataDescription.type) {
          case DataType::scalar:
            preciceInterface->readScalarData(meshName, dataName, couplingCell, dataValues[vertexIndex]);
            break;
          case DataType::vector:
            preciceInterface->readVectorData(meshName, dataName, couplingCell, dataValues[vertexIndex*dim], dataValues[vertexIndex*dim+1], dataValues[vertexIndex*dim+2]);
            break; 
        }
      }
    }
  }


  void writeData() {
    for (auto const& [meshName, data] : _MD2MacroData) {
      for (auto const& [dataName, dataValues] : data) {
        _participant->writeData(meshName, dataName, meshIndices, 0, dataValues);
      }
    }
  }

private:
  /**
   * Add a cell with indice idx to the couplingCells container
   * Compute the coordinates of this cell and add it to the mesh named meshName
   * There mignt be an offset for this cell, only taken into account for the
   * calculation of the coordinates
   */
  void addCell(I01 idx, coupling::datastructures::FlexibleCellContainer<3>& couplingCells,
    const std::string& meshName, std::map<std::string, std::vector<double>>& meshes,
    std::map<std::string, std::map<int, unsigned int>>& vertexToCouplingCell, std::map<std::string, std::map<unsigned int, int>>& couplingCellToVertex,
    tarch::la::Vector<dim, double> offset) {
    
    auto mdDomainOffset = IndexingService::getGlobalMDDomainOffset();
    auto couplingCellSize = IndexingService::getCouplingCellSize();
    auto cellMidPoint = idx.getCellMidPoint();

    std::vector<double> coordinates;
    for (unsigned int currentDim = 0; currentDim < dim; currentDim++) {
      // meshes[meshName].push_back(cellMidPoint[currentDim] + offset[currentDim]);
      coordinates.push_back(cellMidPoint[currentDim] + offset[currentDim]);
    }
    int index = _participant->setMeshVertex(meshName, coordinates);
    indices[meshName].push_back(index);

    coupling::datastructures::CouplingCell<3>* couplingCell = new coupling::datastructures::CouplingCell<3>();
    if (couplingCell == nullptr)
      throw std::runtime_error(std::string("ERROR PreciceAdater::addCell: couplingCell==NULL!"));
    couplingCells << std::make_pair(couplingCell, idx);

    vertexToCouplingCell[meshName][meshes[meshName].size()/dim - 1] = cells.size()-1;
    couplingCellToVertex[meshName][cells.size()-1] = meshes[meshName].size()/dim-1;
  }

  /**
   * For all the given meshes:
   * - call preCICE method setMeshVertices to calculate the indices of the coordinates of each mesh
   * - initialize all the vector which will contain all data for each
   */
  void setMeshesVertices(std::map<std::string, std::vector<double>>& meshes, std::map<std::string, std::map<std::string, std::vector<double>>>& data,
    coupling::preciceadapter::PreciceInterface<dim>* preciceInterface) {
    for (auto const& [meshName, meshCoordinates] : meshes) {
      size_t dataSize = coordinates.size();
      for (const DataDescription& dataDescription : preciceInterface->getDataDescription(name)) {
        if (dataDescription.type == DataType::scalar) dataSize/=dim;
        data[meshName][data.name] = std::vector<double>(dataSize);
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
      if (_participant->requiresMeshConnectivityFor(itVertexIndices->first)) {
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
  }

  // rank of this adapter
  int _rank;
  // pointer to the preCICE participant
  precice::Participant* _participant = nullptr;
  // macro to micro preCICE data containers
  std::map<std::string, std::vector<double>> _Macro2MDMeshes;
  std::map<std::string, std::vector<int>> _Macro2MDIndices;
  std::map<std::string, std::map<std::string, std::vector<double>>> _Macro2MDData;
  // micro to macro preCICE data containers
  std::map<std::string, std::vector<double>> _MD2MacroMeshes;
  std::map<std::string, std::vector<int>> _MD2MacroIndices;
  std::map<std::string, std::map<std::string, std::vector<double>>> _MD2MacroData;
  // macro to micro mapping between vertex arrays and cell array
  std::map<std::string, std::map<int, unsigned int>> _Macro2MDVertexToCell;
  std::map<std::string, std::map<unsigned int, int>> _Macro2MDCouplingCellToVertex;
  // micro to macro mapping between vertex arrays and cell array
  std::map<std::string, std::map<int, unsigned int>> _MD2MacroVertexToCouplingCell;
  std::map<std::string, std::map<unsigned int, int>> _MD2MacroCouplingCellToVertex;
  // macro to micro MaMiCo data buffer
  coupling::datastructures::FlexibleCellContainer<3> _Macro2MDCouplingCells;
  // micro to macro MaMiCo data buffer
  coupling::datastructures::FlexibleCellContainer<3> _MD2MacroCouplingCells;
};