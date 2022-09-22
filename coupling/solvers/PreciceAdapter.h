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
  PreciceAdapter(const tarch::la::Vector<3, double> mdDomainOffset, const tarch::la::Vector<3, double> macroscopicCellSize, const unsigned int overlap, const int rank) : _rank(rank), _overlap(overlap), _mdDomainOffset(mdDomainOffset), 
  _macroscopicCellSize(macroscopicCellSize) {}

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

  double initialize(const unsigned int* const M2mCellGlobalIndices, size_t numberOfM2mCells,
    const unsigned int* const m2MCellGlobalIndices, size_t numberOfm2MCells) {
    int size = 1;
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
    MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
    _solverInterface = new precice::SolverInterface("mamico", "../precice-config.xml", _rank, size);

    std::vector<double> coordsM2mCells;
    _numberOfM2mCells = numberOfM2mCells;
    for (size_t i = 0; i < _numberOfM2mCells; i++) {
      const unsigned int M2mCellGlobalIndex = M2mCellGlobalIndices[i];
      tarch::la::Vector<dim, int> M2mCellGlobalVectorIndex = coupling::indexing::convertToVector<dim>({M2mCellGlobalIndex});
      for (unsigned int currentDim = 0; currentDim < dim; currentDim++) {
        coordsM2mCells.push_back(_mdDomainOffset[currentDim] + M2mCellGlobalVectorIndex[currentDim] * _macroscopicCellSize[currentDim] 
        - _macroscopicCellSize[currentDim] + 0.5 * _macroscopicCellSize[currentDim]);
      }
    }
    _coordsM2mCells = new double[coordsM2mCells.size()];
    std::copy(coordsM2mCells.begin(), coordsM2mCells.end(), _coordsM2mCells);
    _vertexM2mCellIDs = new int[_numberOfM2mCells];
    _solverInterface->setMeshVertices(_solverInterface->getMeshID("mamico-M2m-mesh"), _numberOfM2mCells, _coordsM2mCells, _vertexM2mCellIDs);
    _velocityM2mCells = new double[_numberOfM2mCells * dim];

    std::vector<double> coordsm2MCells;
    _numberOfm2MCells = numberOfm2MCells;
    for (size_t i = 0; i < _numberOfm2MCells; i++) {
      const unsigned int m2MCellGlobalIndex = m2MCellGlobalIndices[i];
      tarch::la::Vector<dim, int> m2MCellGlobalVectorIndex = coupling::indexing::convertToVector<dim>({m2MCellGlobalIndex});
      for (unsigned int currentDim = 0; currentDim < dim; currentDim++) {
        coordsm2MCells.push_back(_mdDomainOffset[currentDim] + m2MCellGlobalVectorIndex[currentDim] * _macroscopicCellSize[currentDim] 
        - _macroscopicCellSize[currentDim] + 0.5 * _macroscopicCellSize[currentDim]);
      }
    }
    _coordsm2MCells = new double[coordsm2MCells.size()];
    std::copy(coordsm2MCells.begin(), coordsm2MCells.end(), _coordsm2MCells);
    _vertexm2MCellIDs = new int[_numberOfm2MCells];
    _solverInterface->setMeshVertices(_solverInterface->getMeshID("mamico-m2M-mesh"), _numberOfm2MCells, _coordsm2MCells, _vertexm2MCellIDs);
    _velocitym2MCells = new double[_numberOfm2MCells * dim];
    for (size_t i = 0; i < _numberOfm2MCells * dim; i++) _velocitym2MCells[i] = 0.0;
    

    /* vertexPreciceID = 0;
    for (size_t im2MCellGlobalIndices = 0; im2MCellGlobalIndices < numberOfm2MCells; im2MCellGlobalIndices++) {
      const unsigned int mamicoRecvdIndex = m2MCellGlobalIndices[im2MCellGlobalIndices];
      tarch::la::Vector<3, int> mamicoRecvdVIndex = coupling::indexing::convertToVector<dim>({mamicoRecvdIndex});
      if (receiveMacroscopicQuantityFromMDSolver(static_cast<tarch::la::Vector<dim, unsigned int>>(mamicoRecvdVIndex))) {
        std::vector<int> tetrahedronsNodeIDs = getTetrahedronsNodeIDs(vertexPreciceID, mamicoRecvdVIndex, m2MCellGlobalIndices, numberOfm2MCells);
        for (size_t tetrahedronIndex = 0; tetrahedronIndex < tetrahedronsNodeIDs.size()/4; tetrahedronIndex++) {
          _solverInterface->setMeshTetrahedron(_solverInterface->getMeshID("mamico-m2M-mesh"),
          tetrahedronsNodeIDs[tetrahedronIndex*4], tetrahedronsNodeIDs[tetrahedronIndex*4+1], tetrahedronsNodeIDs[tetrahedronIndex*4+2], tetrahedronsNodeIDs[tetrahedronIndex*4+3]);
        }
        vertexPreciceID++;
      }
    } */
    return _solverInterface->initialize();
  }

/*   std::vector<int> getTetrahedronsNodeIDs(const int nodePreciceID, const tarch::la::Vector<3, int> nodeMamicoVIndex, const unsigned int* const m2MCellGlobalIndices, size_t numberOfm2MCells) {
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
        while (neighborPreciceID < (int)numberOfm2MCells && (int)m2MCellGlobalIndices[neighborPreciceID] != neighborMamicoIndex) neighborPreciceID++;
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
                std::cout << "Reading data" << std::endl;
  }

  bool isWriteDataRequired(const double dt) {
    return (_solverInterface->isCouplingOngoing() && _solverInterface->isWriteDataRequired(dt));
  }

  void writeData() {
    _solverInterface->writeBlockVectorData(_solverInterface->getDataID("VelocityMicro", _solverInterface->getMeshID("mamico-m2M-mesh")), _numberOfm2MCells, _vertexm2MCellIDs, 
                                    _velocitym2MCells);      
  }

  void setMDBoundaryValues(std::vector<coupling::datastructures::MacroscopicCell<3>*>& m2MBuffer, const unsigned int* const m2MCellGlobalIndices) {
    for (size_t i = 0; i < _numberOfm2MCells; i++) {
      tarch::la::Vector<3, double> vel((1.0 / m2MBuffer[i]->getMacroscopicMass()) * m2MBuffer[i]->getMacroscopicMomentum());
      for (unsigned int currentDim = 0; currentDim < dim; currentDim++) _velocitym2MCells[dim * i + currentDim] = vel[currentDim];
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
    ranks.push_back(_rank);
    return ranks;
  }


private:
  int _rank;
  const unsigned int _overlap;
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
} // namespace solvers
} // namespace coupling

#endif
