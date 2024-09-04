#pragma once

#include "coupling/CouplingMDDefinitions.h"
#include "coupling/datastructures/CouplingCell.h"
#include "coupling/indexing/IndexingService.h"
#include "coupling/indexing/CellIndex.h"
#include "coupling/scenario/precice/PreciceInterface.h"
#include "precice/precice.hpp"
#include "tarch/la/Vector.h"
#include "tarch/utils/Utils.h"
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
  template <class LinkedCell, unsigned int dim> class PreciceAdapter;
} // namespace precice
} // namespace coupling

/**
 * MaMiCo adapter for the preCICE library
 * 'Vertices' are used in the preCICE context (elements constituting the preCICE mesh)
 * 'Cells' are used in the MaMiCo context (elements consituting the MaMiCo cartesian grid)
 * Indexing system should be initialized prior to using this class
 */
template <class LinkedCell, unsigned int dim> class coupling::preciceadapter::PreciceAdapter {
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

  /**
   * Construct the coupling cell container buffers used for
   * communication between preCICE and MaMiCo
   * Construct the MaMiCo cartesian meshes used to communicate
   * between the participants of the coupling and send it to 
   * preCICE
   * Initialize the data buffers
   * Finally, call the preCICE initialize method
   */
  void initialize(coupling::preciceadapter::PreciceInterface<dim>* preciceInterface) { 
    initializeMeshes(preciceInterface);
    _participant->initialize(); 
  }

  /**
   * Calls preCICE isCouplingOnGoing() method of this participant
   */
  bool isCouplingOngoing() const { return _participant->isCouplingOngoing(); }

  /**
   * Calls preCICE isTimeWindowComplete() method of this participant
   */
  bool isTimeWindowComplete() const { return _participant->isTimeWindowComplete(); }

  /**
   * Calls preCICE getMaxTimeStepSize() method of this participant
   */
  double getMaxTimeStepSize() const { return _participant->getMaxTimeStepSize(); }

  /**
   * Calls preCICE requiresWritingCheckpoint() method of this participant
   */
  bool requiresWritingCheckpoint() const { return _participant->requiresWritingCheckpoint(); }

  /**
   * Calls preCICE requiresReadingCheckpoint() method of this participant
   */
  bool requiresReadingCheckpoint() const { return _participant->requiresReadingCheckpoint(); }

  /**
   * Calls preCICE advance() method of this participant
   */
  void advance(const double dt) { _participant->advance(dt); }

  /**
   * Read ´the data received from other participants of the coupling
   * and store it in the data buffers
   * Send the data stored in the data buffers to
   * MaMiCo coupling cells
   */
  void sendFromMacro2MD(coupling::preciceadapter::PreciceInterface<dim>* preciceInterface, coupling::services::MultiMDCellService<LinkedCell, dim>* multiMDCellService) {
    for (auto const& [meshName, data] : _macro2MDData) {
      for (auto& [dataName, value] : data) {
        // need this cast from Tuple_element to std::vector, don't get why yet
        std::vector<double> dataValues = value;
        _participant->readData(meshName, dataName, _macro2MDIndices[meshName], 0, dataValues);
      }
    }
    for (auto pair : _macro2MDCouplingCells) { // loop over all the macro to micro coupling cells
      I01 idx;
      coupling::datastructures::CouplingCell<3>* couplingCell;
      std::tie(couplingCell, idx) = pair;
      // this cell belongs to one macro to micro mesh
      // and is mapped to one vertex of this mesh
      std::string meshName;
      unsigned int vertexIndex;
      std::tie(meshName, vertexIndex) = _macro2MDCellMapping[I00{idx}.get()];
      for (auto const& [dataName, dataValues] : _macro2MDData[meshName]) {
        auto dataDescription = preciceInterface->getDataDescription(meshName, dataName);
        switch (dataDescription.type) {
          case DataType::scalar:
            preciceInterface->readScalarData(meshName, dataName, couplingCell, idx, dataValues[vertexIndex]);
            break;
          case DataType::vector:
            preciceInterface->readVectorData(meshName, dataName, couplingCell, idx, dataValues[vertexIndex*dim], dataValues[vertexIndex*dim+1], dataValues[vertexIndex*dim+2]);
            break; 
        }
      }
    }
    multiMDCellService->sendFromMacro2MD(_macro2MDCouplingCells);
  }

  void sendFromMD2Macro(coupling::preciceadapter::PreciceInterface<dim>* preciceInterface, coupling::services::MultiMDCellService<LinkedCell, dim>* multiMDCellService) {
    multiMDCellService->sendFromMD2Macro(_MD2MacroCouplingCells);
    for (auto pair : _MD2MacroCouplingCells) { // loop over all the macro to micro coupling cells
      I01 idx;
      coupling::datastructures::CouplingCell<3>* couplingCell;
      std::tie(couplingCell, idx) = pair;
      // this cell belongs to one macro to micro mesh
      std::string meshName;
      unsigned int vertexIndex;
      std::tie(meshName, vertexIndex) = _MD2MacroCellMapping[I00{idx}.get()];
      for (auto const& [dataName, dataValues] : _MD2MacroData[meshName]) {
        auto dataDescription = preciceInterface->getDataDescription(meshName, dataName);
        switch (dataDescription.type) {
          case DataType::scalar:
            preciceInterface->writeScalarData(meshName, dataName, couplingCell, idx, dataValues[vertexIndex]);
            break;
          case DataType::vector:
            preciceInterface->writeVectorData(meshName, dataName, couplingCell, idx, dataValues[vertexIndex*dim], dataValues[vertexIndex*dim+1], dataValues[vertexIndex*dim+2]);
            break; 
        }
      }
    }
    for (auto const& [meshName, data] : _MD2MacroData) {
      for (auto const& [dataName, value] : data) {
        std::vector<double> dataValues = value;
        _participant->writeData(meshName, dataName, _MD2MacroIndices[meshName], dataValues);
      }
    }
  }

private:
  /**
   * Construct the coupling cell container buffers used for
   * communication (receive and send according to the given preCICE interface) 
   * between this preCICE adapter and MaMiCo
   * Construct the MaMiCo cartesian meshes used to communicate
   * between this preCICE adapter and the other participant preCICE
   * adapter (adapter for a CFD solver)
   */
  void initializeMeshes(coupling::preciceadapter::PreciceInterface<dim>* preciceInterface) {
    using namespace coupling::indexing;
    for (auto idx : I08()) {
      if (!I12::contains(idx)) {
        if (tarch::utils::contains(preciceInterface->getSourceRanks(idx), (unsigned int)_rank)) {
          addCell(idx, _macro2MDCouplingCells);
          addMeshVertex(preciceInterface->getMacro2MDSolverMeshName(idx), _macro2MDMeshes, _macro2MDIndices, _macro2MDCellMapping, preciceInterface->getMacro2MDSolverMeshOffset(idx));
        }
      }
    }
    for (auto idx : I12()) {
      if (tarch::utils::contains(preciceInterface->getTargetRanks(idx), (unsigned int)_rank)) {
        addCell(idx, _MD2MacroCouplingCells);
        addMeshVertex(preciceInterface->getMD2MacroSolverMeshName(idx), _MD2MacroMeshes, _MD2MacroIndices, _MD2MacroCellMapping, preciceInterface->getMD2MacroSolverMeshOffset(idx)); 
      }
    }
    initializeData(_macro2MDMeshes, _macro2MDData, preciceInterface);
    // setMeshTetrahedra(_M2mVertexIndices, M2mCellIndices, _M2mVertexToCell, _M2mCellToVertex);
    if (preciceInterface->twoWayCoupling()) {
      initializeData(_MD2MacroMeshes, _MD2MacroData, preciceInterface);
      // setMeshTetrahedra(_m2MVertexIndices, m2MCellIndices, _m2MVertexToCell, _m2MCellToVertex);
    }
  }

  /**
   * Add a cell with indice idx to the couplingCells container
   */
  void addCell(I01 idx, coupling::datastructures::FlexibleCellContainer<3>& couplingCells) {
    coupling::datastructures::CouplingCell<3>* couplingCell = new coupling::datastructures::CouplingCell<3>();
    if (couplingCell == nullptr)
      throw std::runtime_error(std::string("ERROR PreciceAdater::addCell: couplingCell==NULL!"));
    couplingCells << std::make_pair(couplingCell, idx);
  }

  /**
   * Add a mesh vertex, calls preCICE setMeshVertex method
   */
  void addMeshVertex(I01 idx, const std::string& meshName, std::map<std::string, std::vector<double>>& meshes, std::map<std::string, std::vector<int>> indices,
    std::map<unsigned int, std::pair<std::string, int>> cellMapping, tarch::la::Vector<dim, double> offset) {
    auto cellMidPoint = idx.getCellMidPoint();
    std::vector<double> coordinates;
    for (unsigned int currentDim = 0; currentDim < dim; currentDim++) {
      auto coordinate = cellMidPoint[currentDim] + offset[currentDim];
      meshes[meshName].push_back(coordinate);
      coordinates.push_back(coordinate);
    }
    int index = _participant->setMeshVertex(meshName, coordinates);
    indices[meshName].push_back(index);
    cellMapping[I00{idx}.get()] = std::make_pair(meshName, index);
  }

  /**
   * For all the given meshes, initialize data buffers used to communicate data with 
   * the other participant preCICE adapter
   */
  void initializeData(std::map<std::string, std::vector<double>>& meshes, std::map<std::string, std::map<std::string, std::vector<double>>>& data,
    coupling::preciceadapter::PreciceInterface<dim>* preciceInterface) {
    for (auto const& [meshName, meshCoordinates] : meshes) {
      size_t dataSize = meshCoordinates.size();
      for (const DataDescription& dataDescription : preciceInterface->getDataDescription(meshName)) {
        if (dataDescription.type == DataType::scalar) dataSize/=dim;
        data[meshName][dataDescription.name] = std::vector<double>(dataSize);
      }
    }
  }

  // /**
  //  * set the participant mesh tetrahedra
  //  */
  // void setMeshTetrahedra(const std::map<std::string, std::vector<int>>& vertexIndices, const std::vector<unsigned int>& cellIndices, 
  //   const std::map<std::string, std::map<int, unsigned int>>& vertexToCell, const std::map<std::string, std::map<unsigned int, int>>& cellToVertex) {
  //   std::map<std::string, std::vector<int>>::const_iterator itVertexIndices;
  //   for (itVertexIndices = vertexIndices.begin(); itVertexIndices != vertexIndices.end(); ++itVertexIndices) {
  //     if (_participant->requiresMeshConnectivityFor(itVertexIndices->first)) {
  //       for (size_t i = 0; i < itVertexIndices->second.size(); ++i) {
  //         using CellIndex_v = coupling::indexing::CellIndex<dim, coupling::indexing::IndexTrait::vector>;
  //         using CellIndex_s = coupling::indexing::CellIndex<dim>;
  //         CellIndex_v cellIndex_v = CellIndex_s{cellIndices[vertexToCell.at(itVertexIndices->first).at(i)]};
  //         const int numberOfNeighbors = 7;
  //         tarch::la::Vector<3, int> directions[numberOfNeighbors] = {{1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}, {0, 0, 1}};
  //         int neighborVertexIndices[numberOfNeighbors];
  //         for (size_t j = 0; j < numberOfNeighbors; j++) {
  //           CellIndex_s neighborCellIndex_s = CellIndex_v{cellIndex_v.get() + directions[j]};
  //           std::vector<unsigned int>::const_iterator itCells = std::find(cellIndices.begin(), cellIndices.end(), neighborCellIndex_s.get());
  //           if (itCells != cellIndices.end()) {
  //             neighborVertexIndices[j] = itVertexIndices->second[cellToVertex.at(itVertexIndices->first).at(itCells - cellIndices.begin())];
  //           } else {
  //             neighborVertexIndices[j] = -1;
  //           }
  //         }
  //         int vertexIndex = itVertexIndices->second[i];
  //         if (neighborVertexIndices[0] > -1 && neighborVertexIndices[5] > -1 && neighborVertexIndices[6] > -1)
  //           _participant->setMeshTetrahedron(itVertexIndices->first, vertexIndex, neighborVertexIndices[0], neighborVertexIndices[5], neighborVertexIndices[6]);
  //         if (neighborVertexIndices[0] > -1 && neighborVertexIndices[2] > -1 && neighborVertexIndices[5] > -1)
  //           _participant->setMeshTetrahedron(itVertexIndices->first, vertexIndex, neighborVertexIndices[0], neighborVertexIndices[2], neighborVertexIndices[5]);
  //         if (neighborVertexIndices[0] > -1 && neighborVertexIndices[5] > -1 && neighborVertexIndices[6] > -1 && neighborVertexIndices[3] > -1)
  //           _participant->setMeshTetrahedron(itVertexIndices->first, neighborVertexIndices[0], neighborVertexIndices[5], neighborVertexIndices[6], neighborVertexIndices[3]);
  //         if (neighborVertexIndices[0] > -1 && neighborVertexIndices[1] > -1 && neighborVertexIndices[2] > -1 && neighborVertexIndices[5] > -1)
  //           _participant->setMeshTetrahedron(itVertexIndices->first, neighborVertexIndices[0], neighborVertexIndices[1], neighborVertexIndices[2], neighborVertexIndices[5]);
  //         if (neighborVertexIndices[0] > -1 && neighborVertexIndices[1] > -1 && neighborVertexIndices[4] > -1 && neighborVertexIndices[5] > -1)
  //           _participant->setMeshTetrahedron(itVertexIndices->first, neighborVertexIndices[0], neighborVertexIndices[1], neighborVertexIndices[4], neighborVertexIndices[5]);
  //         if (neighborVertexIndices[0] > -1 && neighborVertexIndices[3] > -1 && neighborVertexIndices[4] > -1 && neighborVertexIndices[5] > -1)
  //           _participant->setMeshTetrahedron(itVertexIndices->first, neighborVertexIndices[0], neighborVertexIndices[3], neighborVertexIndices[4], neighborVertexIndices[5]);
  //       }
  //     }
  //   }
  // }

  // rank of this adapter
  int _rank;
  // pointer to the preCICE participant
  precice::Participant* _participant = nullptr;
  // macro to micro preCICE data containers
  std::map<std::string, std::vector<double>> _macro2MDMeshes;
  std::map<std::string, std::vector<int>> _macro2MDIndices;
  std::map<std::string, std::map<std::string, std::vector<double>>> _macro2MDData;
  // micro to macro preCICE data containers
  std::map<std::string, std::vector<double>> _MD2MacroMeshes;
  std::map<std::string, std::vector<int>> _MD2MacroIndices;
  std::map<std::string, std::map<std::string, std::vector<double>>> _MD2MacroData;
  // macro to micro mapping from a cell to a mesh
  std::map<unsigned int, std::pair<std::string, int>> _macro2MDCellMapping;
  // micro to macro mapping from a cell to a mesh
  std::map<unsigned int, std::pair<std::string, int>> _MD2MacroCellMapping;
  // macro to micro MaMiCo data buffer
  coupling::datastructures::FlexibleCellContainer<3> _macro2MDCouplingCells;
  // micro to macro MaMiCo data buffer
  coupling::datastructures::FlexibleCellContainer<3> _MD2MacroCouplingCells;
};