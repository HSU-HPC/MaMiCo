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
  PreciceAdapter(coupling::preciceadapter::PreciceInterface<dim>* _preciceInterface): _preciceInterface(_preciceInterface), _rank(0) {
    if (_preciceInterface == nullptr)
      throw new std::runtime_error("PreciceAdapter: preCICE interface == nullptr");
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
  void initialize() { 
    init();
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
   * Call the participant readData, send the values stored in the data buffer to the participant
   */
  void readData() {
    for (auto& meshData : _macro2MDData) {
      for (auto& dataValues : meshData.second) {
        _participant->readData(meshData.first, dataValues.first, _macro2MDIndices[meshData.first], 0, dataValues.second);
      }
    }
  }

  /**
   * Call the participant writeData, store the values in the data buffer
   */
  void writeData() {
    for (auto const& meshData : _MD2MacroData) {
      for (auto const& dataValues : meshData.second) {
        const std::string dataName = dataValues.first;
        _participant->writeData(meshData.first, dataValues.first, _MD2MacroIndices[meshData.first], dataValues.second);
      }
    }
  }

  /**
   * Read ´the data received from other participants of the coupling
   * and store it in the data buffers
   * Send the data stored in the data buffers to
   * MaMiCo coupling cells
   */
  void sendFromMacro2MD(coupling::services::MultiMDCellService<LinkedCell, dim>* multiMDCellService) {
    for (auto& pair : _macro2MDCouplingCells) { // loop over all the macro to micro coupling cells
      I01 idx;
      coupling::datastructures::CouplingCell<3>* couplingCell;
      std::tie(couplingCell, idx) = pair;
      // this cell belongs to one macro to micro mesh
      // and is mapped to one vertex of this mesh
      std::string meshName;
      unsigned int vertexIndex;
      std::tie(meshName, vertexIndex) = _macro2MDCellMapping[I00{idx}.get()];
      for (auto const& dataValues : _macro2MDData[meshName]) {
        auto dataDescription = _preciceInterface->getDataDescription(meshName, dataValues.first);
        switch (dataDescription.type) {
          case DataType::scalar:
            _preciceInterface->readScalarData(meshName, dataValues.first, couplingCell, idx, dataValues.second[vertexIndex]);
            break;
          case DataType::vector:
            _preciceInterface->readVectorData(meshName, dataValues.first, couplingCell, idx, 
              dataValues.second[vertexIndex*dim], dataValues.second[vertexIndex*dim+1], dataValues.second[vertexIndex*dim+2]);
            break; 
        }
      }
    }
    multiMDCellService->sendFromMacro2MD(_macro2MDCouplingCells);
  }

  void sendFromMD2Macro(coupling::services::MultiMDCellService<LinkedCell, dim>* multiMDCellService) {
    multiMDCellService->sendFromMD2Macro(_MD2MacroCouplingCells);
    for (auto pair : _MD2MacroCouplingCells) { // loop over all the macro to micro coupling cells
      I01 idx;
      coupling::datastructures::CouplingCell<3>* couplingCell;
      std::tie(couplingCell, idx) = pair;
      // this cell belongs to one macro to micro mesh
      std::string meshName;
      unsigned int vertexIndex;
      std::tie(meshName, vertexIndex) = _MD2MacroCellMapping[I00{idx}.get()];
      for (auto& [dataName, dataValues] : _MD2MacroData[meshName]) {
        auto dataDescription = _preciceInterface->getDataDescription(meshName, dataName);
        switch (dataDescription.type) {
          case DataType::scalar:
            _preciceInterface->writeScalarData(meshName, dataName, couplingCell, idx, dataValues[vertexIndex]);
            break;
          case DataType::vector:
            _preciceInterface->writeVectorData(meshName, dataName, couplingCell, idx, dataValues[vertexIndex*dim], dataValues[vertexIndex*dim+1], dataValues[vertexIndex*dim+2]);
            break; 
        }
      }
    }
  }

  /**
   * Write the MD to macro coupling cell values to a csv file
   */
  void write2csv(int couplingCycle) {
    if (_MD2MacroCouplingCells.size() > 0) {
      std::stringstream ss;
      ss << "CouetteAvgMultiMDCells_0_" << _rank << "_" << couplingCycle << ".csv";
      std::ofstream file(ss.str().c_str());
      if (!file.is_open()) {
        exit(EXIT_FAILURE);
      }
      file << "I01_x,I01_y,I01_z,vel_x,vel_y,vel_z,mass" << std::endl;
      for (auto pair : _MD2MacroCouplingCells) {
        I01 idx;
        coupling::datastructures::CouplingCell<3>* couplingCell;
        std::tie(couplingCell, idx) = pair;
        tarch::la::Vector<3, double> vel(couplingCell->getMacroscopicMomentum());
        if (couplingCell->getMacroscopicMass() != 0.0)
          vel = (1.0 / couplingCell->getMacroscopicMass()) * vel;
        file   << idx.get()[0] << "," << idx.get()[1] << "," << idx.get()[2] << ","
               << vel[0] << "," << vel[1] << "," << vel[2] << "," 
               << couplingCell->getMacroscopicMass();
        file   << std::endl;
      }
      file.close();
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
  void init() {
    using namespace coupling::indexing;
    for (auto idx : I01()) {
      if (_preciceInterface->isMacro2MD(idx)) {
        if (tarch::utils::contains(_preciceInterface->getSourceRanks(idx), (unsigned int)_rank)) {
          addCell(idx, _macro2MDCouplingCells);
          std::string meshName = _preciceInterface->getMacro2MDSolverMeshName(idx);
          std::vector<double> coordinates = addVertex(idx, meshName, _macro2MDMeshes,
            _preciceInterface->getMacro2MDSolverMeshOffset(idx));
          int index = setMeshVertex(meshName, coordinates, _macro2MDIndices);
          setCellMapping(idx, index, meshName, _macro2MDCellMapping);
        }
      } 
      if (_preciceInterface->isMD2Macro(idx)) {
        if (tarch::utils::contains(_preciceInterface->getTargetRanks(idx), (unsigned int)_rank)) {
          addCell(idx, _MD2MacroCouplingCells);
          std::string meshName = _preciceInterface->getMD2MacroSolverMeshName(idx);
          std::vector<double> coordinates = addVertex(idx, meshName, _MD2MacroMeshes,
            _preciceInterface->getMacro2MDSolverMeshOffset(idx));
          if (_preciceInterface->twoWayCoupling()) {
            int index = setMeshVertex(meshName, coordinates, _MD2MacroIndices);
            setCellMapping(idx, index, meshName, _MD2MacroCellMapping); 
          }
        }    
      }
    }
    initializeData(_macro2MDMeshes, _macro2MDData);
    setMeshTetrahedra(_macro2MDCouplingCells, _macro2MDCellMapping);
    initializeData(_MD2MacroMeshes, _MD2MacroData);
    if (_preciceInterface->twoWayCoupling()) {
      setMeshTetrahedra(_MD2MacroCouplingCells, _MD2MacroCellMapping);
    }
  }

  /**
   * Add a cell with indice idx to the couplingCells container
   */
  void addCell(const I01& idx, coupling::datastructures::FlexibleCellContainer<3>& couplingCells) {
    coupling::datastructures::CouplingCell<3>* couplingCell = new coupling::datastructures::CouplingCell<3>();
    if (couplingCell == nullptr)
      throw std::runtime_error(std::string("ERROR PreciceAdater::addCell: couplingCell==NULL!"));
    couplingCells << std::make_pair(couplingCell, idx);
  }

  /**
   * Add a mesh vertex, calls preCICE setMeshVertex method
   */
  std::vector<double> addVertex(const I01& idx, const std::string& meshName, std::map<std::string, std::vector<double>>& meshes, const tarch::la::Vector<dim, double>& offset) {
    auto cellMidPoint = idx.getCellMidPoint();
    std::vector<double> coordinates;
    for (unsigned int currentDim = 0; currentDim < dim; currentDim++) {
      auto coordinate = cellMidPoint[currentDim] + offset[currentDim];
      meshes[meshName].push_back(coordinate);
      coordinates.push_back(coordinate);
    }
    return coordinates;
  }

  /**
   * Call the participant setMeshVertex method and add the returned index in the indices map
   */
  int setMeshVertex(std::string meshName, std::vector<double> coordinates, std::map<std::string, std::vector<int>>& indices) {
    int index = -1;
    index = _participant->setMeshVertex(meshName, coordinates);
    indices[meshName].push_back(index);
    return index;
  }

  /**
   * Compute the cell mapping between the MaMico index and the participant index for this coupling cell
   */
  void setCellMapping(const I01& idx, int index, std::string meshName, std::map<unsigned int, std::pair<std::string, int>>& cellMapping) {
    cellMapping[I00{idx}.get()] = std::make_pair(meshName, index);
  }

  /**
   * For all the given meshes, initialize data buffers used to communicate data with 
   * the other participant preCICE adapter
   */
  void initializeData(const std::map<std::string, std::vector<double>>& meshes, std::map<std::string, std::map<std::string, std::vector<double>>>& data) {
    for (auto& meshCoordinates : meshes) {
      const std::string meshName = meshCoordinates.first;
      const std::vector<double> coordinates = meshCoordinates.second;
      size_t dataSize = coordinates.size();
      for (const DataDescription& dataDescription : _preciceInterface->getDataDescriptions(meshName)) {
        if (dataDescription.type == DataType::scalar) dataSize/=dim;
        data[meshName][dataDescription.name] = std::vector<double>(dataSize);
      }
    }
  }

  /**
   * Set the MaMiCo mesh tetrahedra to preCICE
   */
  void setMeshTetrahedra(const coupling::datastructures::FlexibleCellContainer<3>& couplingCells, const std::map<unsigned int, std::pair<std::string, int>>& cellMapping) {
    for (auto pair : couplingCells) {
      I01 idx;
      coupling::datastructures::CouplingCell<3>* couplingCell;
      std::tie(couplingCell, idx) = pair;
      std::string meshName;
      unsigned int vertexIndex;
      std::tie(meshName, vertexIndex) = cellMapping.at(I00{idx}.get());
      if (_participant->requiresMeshConnectivityFor(meshName)) {
          const int numberOfNeighbors = 7;
          tarch::la::Vector<3, int> directions[numberOfNeighbors] = {{1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}, {0, 0, 1}};
          int neighborVertexIndices[numberOfNeighbors];
          for (size_t j = 0; j < numberOfNeighbors; j++) {
            I01 neighborCellIdx{idx.get() + directions[j]};
            neighborVertexIndices[j] = -1;
            if (_preciceInterface->contains(meshName, neighborCellIdx)) {
              int neighborVertexIndex;
              std::tie(std::ignore, neighborVertexIndex) = cellMapping.at(I00{neighborCellIdx}.get());
              neighborVertexIndices[j] = neighborVertexIndex;
            }
          }
          if (neighborVertexIndices[0] > -1 && neighborVertexIndices[5] > -1 && neighborVertexIndices[6] > -1)
            _participant->setMeshTetrahedron(meshName, vertexIndex, neighborVertexIndices[0], neighborVertexIndices[5], neighborVertexIndices[6]);
          if (neighborVertexIndices[0] > -1 && neighborVertexIndices[2] > -1 && neighborVertexIndices[5] > -1)
            _participant->setMeshTetrahedron(meshName, vertexIndex, neighborVertexIndices[0], neighborVertexIndices[2], neighborVertexIndices[5]);
          if (neighborVertexIndices[0] > -1 && neighborVertexIndices[5] > -1 && neighborVertexIndices[6] > -1 && neighborVertexIndices[3] > -1)
            _participant->setMeshTetrahedron(meshName, neighborVertexIndices[0], neighborVertexIndices[5], neighborVertexIndices[6], neighborVertexIndices[3]);
          if (neighborVertexIndices[0] > -1 && neighborVertexIndices[1] > -1 && neighborVertexIndices[2] > -1 && neighborVertexIndices[5] > -1)
            _participant->setMeshTetrahedron(meshName, neighborVertexIndices[0], neighborVertexIndices[1], neighborVertexIndices[2], neighborVertexIndices[5]);
          if (neighborVertexIndices[0] > -1 && neighborVertexIndices[1] > -1 && neighborVertexIndices[4] > -1 && neighborVertexIndices[5] > -1)
            _participant->setMeshTetrahedron(meshName, neighborVertexIndices[0], neighborVertexIndices[1], neighborVertexIndices[4], neighborVertexIndices[5]);
          if (neighborVertexIndices[0] > -1 && neighborVertexIndices[3] > -1 && neighborVertexIndices[4] > -1 && neighborVertexIndices[5] > -1)
            _participant->setMeshTetrahedron(meshName, neighborVertexIndices[0], neighborVertexIndices[3], neighborVertexIndices[4], neighborVertexIndices[5]);
      }
    }
  }

  // Interface between MaMiCo and the participant
  coupling::preciceadapter::PreciceInterface<dim>* _preciceInterface = nullptr;
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