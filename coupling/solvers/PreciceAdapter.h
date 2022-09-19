#ifndef _MOLECULARDYNAMICS_COUPLING_SOLVERS_PRECICESOLVER_H_
#define _MOLECULARDYNAMICS_COUPLING_SOLVERS_PRECICESOLVER_H_

#include "coupling/CouplingMDDefinitions.h"
#include "coupling/IndexConversion.h"
#include "coupling/datastructures/MacroscopicCell.h"
#include "coupling/indexing/CellIndex.h"
#include "coupling/interface/MacroscopicSolverInterface.h"
#include "precice/SolverInterface.hpp"
#include "tarch/la/Vector.h"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>

namespace coupling {
namespace solvers {

using namespace indexing;

template <unsigned int dim> class PreciceAdapter : public coupling::interface::MacroscopicSolverInterface<dim> {
public:
  PreciceAdapter() {}

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
    if (_velocitym2MCells != NULL) {
      delete[] _velocitym2MCells;
    }
  }

  double initialize(const tarch::la::Vector<3, double> mdDomainOffset, const tarch::la::Vector<3, double> macroscopicCellSize, const unsigned int overlap,
    const unsigned int* const mamicoRecvdIndices, size_t mamicoRecvdIndicesSize) {
    _mdDomainOffset = mdDomainOffset;
    _macroscopicCellSize = macroscopicCellSize;
    _overlap = overlap;
    int rank = 0;
    int size = 1;
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
    _solverInterface = new precice::SolverInterface("mamico", "../precice-config.xml", rank, size);

    _numberOfM2mCells = 0;
    for (CellIndex<dim, IndexTrait::vector> cellIndex : CellIndex<dim, IndexTrait::vector>()) {
      if (sendMacroscopicQuantityToMDSolver(static_cast<tarch::la::Vector<dim, unsigned int>>(cellIndex.get())))
        _numberOfM2mCells++;
    }
    std::cout << "MaMiCo: mamico-M2m-mesh size=" << _numberOfM2mCells << std::endl;
    _coordsM2mCells = new double[dim * _numberOfM2mCells];
    unsigned int vertexPreciceID = 0;
    for (CellIndex<dim, IndexTrait::vector> cellIndex : CellIndex<dim, IndexTrait::vector>()) {
      tarch::la::Vector<3, unsigned int> cellVectorIndex = static_cast<tarch::la::Vector<dim, unsigned int>>(cellIndex.get());
      if (sendMacroscopicQuantityToMDSolver(cellVectorIndex)) {
        for (unsigned int currentDim = 0; currentDim < dim; currentDim++) {
          _coordsM2mCells[dim * vertexPreciceID + currentDim] =
              mdDomainOffset[currentDim] + cellVectorIndex[currentDim] * macroscopicCellSize[0] - macroscopicCellSize[0] + 0.5 * macroscopicCellSize[currentDim];
        }
        vertexPreciceID++;
      }
    }
    _vertexM2mCellIDs = new int[_numberOfM2mCells];
    _solverInterface->setMeshVertices(_solverInterface->getMeshID("mamico-M2m-mesh"), _numberOfM2mCells, _coordsM2mCells, _vertexM2mCellIDs);
    _velocityM2mCells = new double[_numberOfM2mCells * dim];

    _coordsm2MCells = new double[dim * _numberOfm2MCells];
    vertexPreciceID = 0;
    for (size_t iMamicoRecvdIndices = 0; iMamicoRecvdIndices < mamicoRecvdIndicesSize; iMamicoRecvdIndices++) {
      const unsigned int mamicoRecvdIndex = mamicoRecvdIndices[iMamicoRecvdIndices];
      tarch::la::Vector<3, unsigned int> mamicoRecvdVIndex = static_cast<tarch::la::Vector<dim, unsigned int>>(coupling::indexing::convertToVector<dim>({mamicoRecvdIndex}));
      if (receiveMacroscopicQuantityFromMDSolver(mamicoRecvdVIndex)) {
        for (unsigned int currentDim = 0; currentDim < dim; currentDim++) {
          _coordsm2MCells[dim * vertexPreciceID + currentDim] =
              _mdDomainOffset[currentDim] + mamicoRecvdVIndex[currentDim] * macroscopicCellSize[0] - macroscopicCellSize[0] + 0.5 * _macroscopicCellSize[currentDim];
        }
        vertexPreciceID++;
      }
    }
    
    _vertexm2MCellIDs = new int[_numberOfm2MCells];
    _solverInterface->setMeshVertices(_solverInterface->getMeshID("mamico-m2M-mesh"), _numberOfm2MCells, _coordsm2MCells, _vertexm2MCellIDs);
    _velocitym2MCells = new double[_numberOfm2MCells * dim];

    /* vertexPreciceID = 0;
    for (size_t iMamicoRecvdIndices = 0; iMamicoRecvdIndices < mamicoRecvdIndicesSize; iMamicoRecvdIndices++) {
      const unsigned int mamicoRecvdIndex = mamicoRecvdIndices[iMamicoRecvdIndices];
      tarch::la::Vector<3, int> mamicoRecvdVIndex = coupling::indexing::convertToVector<dim>({mamicoRecvdIndex});
      if (receiveMacroscopicQuantityFromMDSolver(static_cast<tarch::la::Vector<dim, unsigned int>>(mamicoRecvdVIndex))) {
        std::vector<int> tetrahedronsNodeIDs = getTetrahedronsNodeIDs(vertexPreciceID, mamicoRecvdVIndex, mamicoRecvdIndices, mamicoRecvdIndicesSize);
        for (size_t tetrahedronIndex = 0; tetrahedronIndex < tetrahedronsNodeIDs.size()/4; tetrahedronIndex++) {
          _solverInterface->setMeshTetrahedron(_solverInterface->getMeshID("mamico-m2M-mesh"),
          tetrahedronsNodeIDs[tetrahedronIndex*4], tetrahedronsNodeIDs[tetrahedronIndex*4+1], tetrahedronsNodeIDs[tetrahedronIndex*4+2], tetrahedronsNodeIDs[tetrahedronIndex*4+3]);
        }
        vertexPreciceID++;
      }
    } */

    return _solverInterface->initialize();
  }

/*   std::vector<int> getTetrahedronsNodeIDs(const int nodePreciceID, const tarch::la::Vector<3, int> nodeMamicoVIndex, const unsigned int* const mamicoRecvdIndices, size_t mamicoRecvdIndicesSize) {
    std::vector<int> tetrahedronsNodeIDs;
    const int numberOfNeighbors = 6;
    tarch::la::Vector<3, int> directionNeighbors[numberOfNeighbors] = {{0,0,-1}, {0,0,1}, {1,0,0}, {0,1,0}, {-1,0,0}, {0,-1,0}};
    int neighborVertexIDSs[numberOfNeighbors];
    for (size_t iDirectionNeighbors= 0; iDirectionNeighbors < numberOfNeighbors; iDirectionNeighbors++) {
      tarch::la::Vector<3, int> neighborMamicoVIndex = nodeMamicoVIndex + directionNeighbors[iDirectionNeighbors];
      if (receiveMacroscopicQuantityFromMDSolver(static_cast<tarch::la::Vector<dim, unsigned int>>(neighborMamicoVIndex))) {
        using coupling::indexing::IndexTrait;
        using coupling::indexing::CellIndex;
        using coupling::indexing::convertToScalar;
        int neighborMamicoIndex = convertToScalar<dim>(CellIndex<dim,IndexTrait::vector>(neighborMamicoVIndex));
        int neighborPreciceID = 0;
        while (neighborPreciceID < (int)mamicoRecvdIndicesSize && (int)mamicoRecvdIndices[neighborPreciceID] != neighborMamicoIndex) neighborPreciceID++;
        neighborVertexIDSs[iDirectionNeighbors]=neighborPreciceID;
      } else {
        neighborVertexIDSs[iDirectionNeighbors]=-1;
      }
    }
    for (size_t zNeighborIndex = 0; zNeighborIndex < 2; zNeighborIndex++) {
      if (neighborVertexIDSs[zNeighborIndex] > 0) {
        for (size_t xyNeighborIndex = 0; xyNeighborIndex < 4; xyNeighborIndex++) {
          if (neighborVertexIDSs[xyNeighborIndex+2] > 0 && neighborVertexIDSs[(xyNeighborIndex+1)%4 + 2] > 0) {
            tetrahedronsNodeIDs.insert(tetrahedronsNodeIDs.end(), {nodePreciceID, neighborVertexIDSs[zNeighborIndex],
                        neighborVertexIDSs[xyNeighborIndex+2], neighborVertexIDSs[(xyNeighborIndex+1)%4 + 2]});
          }
        }
      }
    }
    return tetrahedronsNodeIDs;
  } */

  tarch::la::Vector<3, double> getVelocity(tarch::la::Vector<3, double> pos) const {
    tarch::la::Vector<3, double> vel(0.0);
    unsigned int cellIndex = 0;
    while (cellIndex < _numberOfM2mCells &&
           !(_coordsM2mCells[dim * cellIndex] == pos[0] && _coordsM2mCells[dim * cellIndex + 1] == pos[1] && _coordsM2mCells[dim * cellIndex + 2] == pos[2]))
      cellIndex++;
    if (cellIndex < _numberOfM2mCells) {
      for (unsigned int currentDim = 0; currentDim < dim; currentDim++) {
        vel[currentDim] = _velocityM2mCells[dim * cellIndex + currentDim];
      }
    }
    return vel;
  }

  bool isCouplingOngoing() {
    return _solverInterface->isCouplingOngoing();
  }

  double advance(const double dt) {
    return _solverInterface->advance(dt);
  }

  bool isReadDataAvailable() {
    return (_solverInterface->isCouplingOngoing() && _solverInterface->isReadDataAvailable());
  }

  void readData() {
      _solverInterface->readBlockVectorData(_solverInterface->getDataID("VelocityMacro", _solverInterface->getMeshID("mamico-M2m-mesh")), _numberOfM2mCells, _vertexM2mCellIDs,
                                      _velocityM2mCells);
  }

  bool isWriteDataRequired(const double dt) {
    return (_solverInterface->isCouplingOngoing() && _solverInterface->isWriteDataRequired(dt));
  }

  void writeData() {
    _solverInterface->writeBlockVectorData(_solverInterface->getDataID("VelocityMicro", _solverInterface->getMeshID("mamico-m2M-mesh")), _numberOfm2MCells, _vertexm2MCellIDs, 
                                    _velocitym2MCells);      
  }

  void setMDBoundaryValues(std::vector<coupling::datastructures::MacroscopicCell<3>*>& recvBuffer, const unsigned int* const mamicoRecvIndices) {
    const size_t mamicoRecvIndicesSize = recvBuffer.size();
    int cellIndex = 0;
    for (size_t indexRecvIndices = 0; indexRecvIndices < mamicoRecvIndicesSize; indexRecvIndices++) {
      const unsigned int recvIndex = mamicoRecvIndices[indexRecvIndices];
      tarch::la::Vector<3, unsigned int> recvVectorIndex = static_cast<tarch::la::Vector<dim, unsigned int>>(coupling::indexing::convertToVector<dim>({recvIndex}));
      if (receiveMacroscopicQuantityFromMDSolver(recvVectorIndex)) {
        tarch::la::Vector<3, double> vel((1.0 / recvBuffer[indexRecvIndices]->getMacroscopicMass()) * recvBuffer[indexRecvIndices]->getMacroscopicMomentum());
        for (unsigned int currentDim = 0; currentDim < dim; currentDim++) _velocitym2MCells[dim * cellIndex + currentDim] = vel[currentDim];
        cellIndex++;
      }
    }
  }

  bool receiveMacroscopicQuantityFromMDSolver(tarch::la::Vector<dim, unsigned int> globalCellIndex) override {
		tarch::la::Vector<3, int> lowerBoundary = CellIndex<3>::lowerBoundary.get();
		tarch::la::Vector<3, int> upperBoundary = CellIndex<3>::upperBoundary.get();
		bool rcv = true;
		for (unsigned int currentDim = 0; currentDim < dim; currentDim++) {
			rcv &= (int)globalCellIndex[currentDim] >= lowerBoundary[currentDim] + (int)_overlap;
      rcv &= (int)globalCellIndex[currentDim] <= upperBoundary[currentDim] - (int)_overlap;
		}
    return rcv;
  }

  bool sendMacroscopicQuantityToMDSolver(tarch::la::Vector<3, unsigned int> globalCellIndex) override {
		tarch::la::Vector<3, int> lowerBoundary = CellIndex<3>::lowerBoundary.get();
    tarch::la::Vector<3, int> upperBoundary = CellIndex<3>::upperBoundary.get();
    bool isInnerCell = true;
		bool isGhostCell = false;
    for (unsigned int currentDim = 0; currentDim < dim; currentDim++) {
      isInnerCell &= (int)globalCellIndex[currentDim] > lowerBoundary[currentDim] + (int)_overlap;
      isInnerCell &= (int)globalCellIndex[currentDim] < upperBoundary[currentDim] - (int)_overlap;
			isGhostCell |= globalCellIndex[currentDim] == 0;
      isGhostCell |= (int)globalCellIndex[currentDim] == upperBoundary[currentDim];
    }
    return (!isGhostCell) && (!isInnerCell);
  }

  std::vector<unsigned int> getRanks(tarch::la::Vector<dim, unsigned int> globalCellIndex) override {
    std::vector<unsigned int> ranks;
    ranks.push_back(0);
    return ranks;
  }


private:
  unsigned int _overlap;
  tarch::la::Vector<3, double> _mdDomainOffset;
  tarch::la::Vector<3, double> _macroscopicCellSize;

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
} // namespace solvers
} // namespace coupling

#endif
