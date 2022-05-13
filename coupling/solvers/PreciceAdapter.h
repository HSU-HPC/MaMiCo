#ifndef _MOLECULARDYNAMICS_COUPLING_SOLVERS_PRECICESOLVER_H_
#define _MOLECULARDYNAMICS_COUPLING_SOLVERS_PRECICESOLVER_H_

#include "coupling/CouplingMDDefinitions.h"
#include "coupling/IndexConversion.h"
#include "coupling/datastructures/MacroscopicCell.h"
#include "coupling/indexing/CellIndex.h"
#include "coupling/interface/MacroscopicSolverInterface.h"
#include "coupling/solvers/CouetteSolver.h"
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

template <unsigned int dim> class PreciceAdapter : public AbstractCouetteSolver<3>, public interface::MacroscopicSolverInterface<dim> {
public:
  PreciceAdapter(const double channelHeight, const double dx, const double dt, const unsigned int plotEveryTimestep, const std::string filestem,
                 const unsigned int overlap)
      : coupling::solvers::AbstractCouetteSolver<3>(), coupling::interface::MacroscopicSolverInterface<dim>(), _channelHeight(channelHeight), _dx(dx), _dt(dt),
        _plotEveryTimestep(plotEveryTimestep), _filestem(filestem) {
  }

  virtual ~PreciceAdapter() {
    if (_interface != NULL) {
      _interface->finalize();
      delete _interface;
      _interface = NULL;
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

  void setInterface(const unsigned int overlap) { _overlap = overlap; }

  void setCouplingMesh(const tarch::la::Vector<3, double> mdDomainOffset, const tarch::la::Vector<3, double> macroscopicCellSize,
    const unsigned int* const recvIndices, size_t recvIndicesSize) {
    _mdDomainOffset = mdDomainOffset;
    _macroscopicCellSize = macroscopicCellSize;
    int rank = 0;
    int size = 1;
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
    _interface = new precice::SolverInterface("mamico", "../precice-config.xml", rank, size);

    _numberOfM2mCells = 0;
    for (CellIndex<dim, IndexTrait::vector> cellIndex : CellIndex<dim, IndexTrait::vector>()) {
      if (sendMacroscopicQuantityToMDSolver(static_cast<tarch::la::Vector<dim, unsigned int>>(cellIndex.get())))
        _numberOfM2mCells++;
    }
    std::cout << "MaMiCo >> mamico-M2m-mesh size: " << _numberOfM2mCells << std::endl;
    _coordsM2mCells = new double[dim * _numberOfM2mCells];
    unsigned int index = 0;
    for (CellIndex<dim, IndexTrait::vector> cellIndex : CellIndex<dim, IndexTrait::vector>()) {
      tarch::la::Vector<3, unsigned int> cellVectorIndex = static_cast<tarch::la::Vector<dim, unsigned int>>(cellIndex.get());
      if (sendMacroscopicQuantityToMDSolver(cellVectorIndex)) {
        for (unsigned int currentDim = 0; currentDim < dim; currentDim++) {
          _coordsM2mCells[dim * index + currentDim] =
              mdDomainOffset[currentDim] + cellVectorIndex[currentDim] * _dx - _dx + 0.5 * macroscopicCellSize[currentDim];
        }
        index++;
      }
    }
    _vertexM2mCellIDs = new int[_numberOfM2mCells];
    _interface->setMeshVertices(_interface->getMeshID("mamico-M2m-mesh"), _numberOfM2mCells, _coordsM2mCells, _vertexM2mCellIDs);
    _velocityM2mCells = new double[_numberOfM2mCells * dim];

		_numberOfm2MCells = 0;
    for (CellIndex<dim, IndexTrait::vector> cellIndex : CellIndex<dim, IndexTrait::vector>()) {
      if (sendMacroscopicQuantityToPreCICE(static_cast<tarch::la::Vector<dim, unsigned int>>(cellIndex.get())))
        _numberOfm2MCells++;
    }
    std::cout << "MaMiCo >> mamico-m2M-mesh size: " << _numberOfm2MCells << std::endl;
    _coordsm2MCells = new double[dim * _numberOfm2MCells];
    index = 0;
    for (CellIndex<dim, IndexTrait::vector> cellIndex : CellIndex<dim, IndexTrait::vector>()) {
      tarch::la::Vector<3, unsigned int> cellVectorIndex = static_cast<tarch::la::Vector<dim, unsigned int>>(cellIndex.get());
      if (sendMacroscopicQuantityToPreCICE(cellVectorIndex)) {
        for (unsigned int currentDim = 0; currentDim < dim; currentDim++) {
          _coordsm2MCells[dim * index + currentDim] =
              mdDomainOffset[currentDim] + cellVectorIndex[currentDim] * _dx - _dx + 0.5 * macroscopicCellSize[currentDim];
        }
        index++;
      }
    }
    _vertexm2MCellIDs = new int[_numberOfm2MCells];
    _interface->setMeshVertices(_interface->getMeshID("mamico-m2M-mesh"), _numberOfm2MCells, _coordsm2MCells, _vertexm2MCellIDs);
    /*for (size_t i = 0; i < _numberOfm2MCells; i++) {
      std::cout << "vertex id:" << _vertexm2MCellIDs[i]
                << ", coords:[" << _coordsm2MCells[dim* i]
                << "," << _coordsm2MCells[dim* i+1]
                << "," << _coordsm2MCells[dim* i +2]
                << "]" << std::endl;
    }*/
    _velocitym2MCells = new double[_numberOfm2MCells * dim];

    _precice_dt = _interface->initialize();

    for (size_t index = 0; index < recvIndicesSize; index++) {
      const unsigned int recvIndex = recvIndices[index];
      tarch::la::Vector<3, unsigned int> cellVectorIndex = static_cast<tarch::la::Vector<dim, unsigned int>>(coupling::indexing::convertToVector<dim>({recvIndex}));
      if (sendMacroscopicQuantityToPreCICE(cellVectorIndex)) {
        double* coords = new double[dim];
        for (unsigned int currentDim = 0; currentDim < dim; currentDim++) {
          coords[currentDim] =
              _mdDomainOffset[currentDim] + cellVectorIndex[currentDim] * _dx - _dx + 0.5 * _macroscopicCellSize[currentDim];
        }
        unsigned int m2MCellID = 0;
        while(m2MCellID < _numberOfm2MCells && (coords[0]!=_coordsm2MCells[m2MCellID*dim] || coords[1]!=_coordsm2MCells[m2MCellID*dim+1] || coords[2]!=_coordsm2MCells[m2MCellID*dim+2])) m2MCellID++;
        if (m2MCellID >= _numberOfm2MCells) {
          std::cout << "ERROR Could not find received macroscopic cell among precice cells " << std::endl;
          exit(EXIT_FAILURE);
        }
        _vertexRecvIDsTom2MCellIDs.insert({recvIndex,m2MCellID});
      }
    }
    std::cout << "HEYYYYYYYYYYY" << std::endl;
    for (const auto& [key, value] : _vertexRecvIDsTom2MCellIDs) {
        std::cout << '[' << key << "] = " << value << "; ";
    }
    std::cout << std::endl;
  }

  void setWallVelocity(const tarch::la::Vector<3, double> wallVelocity) override { _wallVelocity = wallVelocity; }

  tarch::la::Vector<3, double> getVelocity(tarch::la::Vector<3, double> pos) const override {
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

  double getDensity(tarch::la::Vector<3, double> pos) const { return 1.0; }

  void advance(double dt) override {
    if (_interface->isCouplingOngoing()) {
      if (_interface->isReadDataAvailable()) {
        std::cout << "MaMiCo >> Reading macro velocities from preCICE !" << std::endl;
        // velocity from the conitnuum solver
        _interface->readBlockVectorData(_interface->getDataID("VelocityMacro", _interface->getMeshID("mamico-M2m-mesh")), _numberOfM2mCells, _vertexM2mCellIDs,
                                        _velocityM2mCells);
      }

      // Solving the time step
      // Normally does nothing, everything is done on the MD side
      if (_interface->isWriteDataRequired(dt)) {
        int meshID = _interface->getMeshID("mamico-m2M-mesh");
        std::cout << "MaMiCo >> Writing micro velocities to preCICE !" << std::endl;
        // Velocit from the md solver
        int dataID = _interface->getDataID("VelocityMicro", meshID);
        _interface->writeBlockVectorData(dataID, _numberOfm2MCells, _vertexm2MCellIDs, _velocitym2MCells);
      }

      double computed_dt = std::min(_precice_dt, dt);
      _precice_dt = _interface->advance(computed_dt);
    }
  }

  bool sendMacroscopicQuantityToPreCICE(tarch::la::Vector<dim, unsigned int> globalCellIndex) {
    tarch::la::Vector<3, int> lowerBoundary = CellIndex<3>::lowerBoundary.get();
    tarch::la::Vector<3, int> upperBoundary = CellIndex<3>::upperBoundary.get();
    bool isInnerCell = true;
    bool isOuterCell = false;
    for (unsigned int currentDim = 0; currentDim < dim; currentDim++) {
      isInnerCell &= (int)globalCellIndex[currentDim] > lowerBoundary[currentDim] + (int)_overlap + 1;
      isInnerCell &= (int)globalCellIndex[currentDim] < upperBoundary[currentDim] - (int)_overlap - 1;
      isOuterCell |= (int)globalCellIndex[currentDim] <= lowerBoundary[currentDim] + (int)_overlap - 1;
      isOuterCell |= (int)globalCellIndex[currentDim] >= upperBoundary[currentDim] - (int)_overlap + 1;
    }
    return !isInnerCell && !isOuterCell;
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

  bool sendMacroscopicQuantityToMDSolver(tarch::la::Vector<3, unsigned int> globalCellIndex) {
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


  void setMDBoundaryValues(std::vector<coupling::datastructures::MacroscopicCell<3>*>& recvBuffer, const unsigned int* const recvIndices) {
    std::cout << "setMDBoundaryValues" << std::endl;
    const size_t size = recvBuffer.size();
    for (size_t index = 0; index < size; index++) {
      const unsigned int recvIndex = recvIndices[index];
      if (_vertexRecvIDsTom2MCellIDs.count(recvIndex) > 0) {
        const unsigned int m2MCellID = _vertexRecvIDsTom2MCellIDs[recvIndex];
        tarch::la::Vector<3, double> vel((1.0 / recvBuffer[index]->getMacroscopicMass()) * recvBuffer[index]->getMacroscopicMomentum());
        /*std::cout << "recvIndex: " << recvIndex << std::endl;
        std::cout << "preciceIndex: " << i << std::endl;
        std::cout << "vel: " << vel << std::endl;*/
        for (unsigned int currentDim = 0; currentDim < dim; currentDim++) _velocitym2MCells[dim * m2MCellID + currentDim] = vel[currentDim];
      }
    }
  }


private:
  const double _channelHeight;
  const double _dx;
  const double _dt;
  const unsigned int _plotEveryTimestep;
  const std::string _filestem;
  unsigned int _counter{0};
  tarch::la::Vector<3, double> _wallVelocity;
  unsigned int _overlap;
  tarch::la::Vector<3, double> _mdDomainOffset;
  tarch::la::Vector<3, double> _macroscopicCellSize;

  precice::SolverInterface* _interface = NULL;
  double _precice_dt;

  int* _vertexM2mCellIDs;
  double* _coordsM2mCells;
  unsigned int _numberOfM2mCells;
  double* _velocityM2mCells;

	int* _vertexm2MCellIDs;
  double* _coordsm2MCells;
  unsigned int _numberOfm2MCells;
  double* _velocitym2MCells;

  std::map<unsigned int, unsigned int> _vertexRecvIDsTom2MCellIDs;
};
} // namespace solvers
} // namespace coupling

#endif
