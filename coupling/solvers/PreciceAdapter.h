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
  }

  void setInterface(const unsigned int overlap) { _overlap = overlap; }

  void setCouplingMesh(const tarch::la::Vector<3, double> mdDomainOffset, const tarch::la::Vector<3, double> macroscopicCellSize) {
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
    unsigned int cellIndex = 0;
    for (CellIndex<dim, IndexTrait::vector> cellIndex : CellIndex<dim, IndexTrait::vector>()) {
      tarch::la::Vector<3, unsigned int> cellVectorIndex = static_cast<tarch::la::Vector<dim, unsigned int>>(cellIndex.get());
      if (sendMacroscopicQuantityToMDSolver(cellVectorIndex)) {
        for (unsigned int currentDim = 0; currentDim < dim; currentDim++) {
          _coordsM2mCells[dim * cellIndex + currentDim] =
              mdDomainOffset[currentDim] + cellVectorIndex[currentDim] * _dx - _dx + 0.5 * macroscopicCellSize[currentDim];
        }
        cellIndex++;
      }
    }
    _vertexM2mCellIDs = new int[_numberOfM2mCells];
    _interface->setMeshVertices(_interface->getMeshID("mamico-M2m-mesh"), _numberOfM2mCells, _coordsM2mCells, _vertexM2mCellIDs);
    _velocityM2mCells = new double[_numberOfM2mCells * dim];

		_numberOfm2MCells = 0;
    for (CellIndex<dim, IndexTrait::vector> cellIndex : CellIndex<dim, IndexTrait::vector>()) {
      if (receiveMacroscopicQuantityFromMDSolver(static_cast<tarch::la::Vector<dim, unsigned int>>(cellIndex.get())))
        _numberOfm2MCells++;
    }
    std::cout << "MaMiCo >> mamico-m2M-mesh size: " << _numberOfM2mCells << std::endl;
    _coordsm2MCells = new double[dim * _numberOfM2mCells];
    cellIndex = 0;
    for (CellIndex<dim, IndexTrait::vector> cellIndex : CellIndex<dim, IndexTrait::vector>()) {
      tarch::la::Vector<3, unsigned int> cellVectorIndex = static_cast<tarch::la::Vector<dim, unsigned int>>(cellIndex.get());
      if (receiveMacroscopicQuantityFromMDSolver(cellVectorIndex)) {
        for (unsigned int currentDim = 0; currentDim < dim; currentDim++) {
          _coordsm2MCells[dim * cellIndex + currentDim] =
              mdDomainOffset[currentDim] + cellVectorIndex[currentDim] * _dx - _dx + 0.5 * macroscopicCellSize[currentDim];
        }
        cellIndex++;
      }
    }
    _vertexm2MCellIDs = new int[_numberOfm2MCells];
    _interface->setMeshVertices(_interface->getMeshID("mamico-m2M-mesh"), _numberOfm2MCells, _coordsm2MCells, _vertexm2MCellIDs);
    _precice_dt = _interface->initialize();
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
      double* mdVelocities = new double[_numberOfm2MCells * dim];
      for (unsigned int vertexIndex = 0; vertexIndex < _numberOfm2MCells; vertexIndex++) {
        mdVelocities[dim * vertexIndex] = 1.;
        mdVelocities[dim * vertexIndex + 1] = 1.;
        mdVelocities[dim * vertexIndex + 2] = 1.;
      }

      if (_interface->isWriteDataRequired(dt)) {
        int meshID = _interface->getMeshID("mamico-m2M-mesh");
        std::cout << "MaMiCo >> Writing micro velocities to preCICE !" << std::endl;
        // Velocit from the md solver
        int dataID = _interface->getDataID("VelocityMicro", meshID);
        _interface->writeBlockVectorData(dataID, _numberOfm2MCells, _vertexm2MCellIDs, mdVelocities);
      }
      delete[] mdVelocities;
      double computed_dt = std::min(_precice_dt, dt);
      _precice_dt = _interface->advance(computed_dt);
    }
  }

  bool receiveMacroscopicQuantityFromMDSolver(tarch::la::Vector<dim, unsigned int> globalCellIndex) override {
		tarch::la::Vector<3, int> lowerBoundary = CellIndex<3>::lowerBoundary.get();
		tarch::la::Vector<3, int> upperBoundary = CellIndex<3>::upperBoundary.get();
		bool isInnerCell = true;
		bool isOuterCell = false;
		for (unsigned int currentDim = 0; currentDim < dim; currentDim++) {
			isInnerCell &= (int)globalCellIndex[currentDim] > lowerBoundary[currentDim] + (int)_overlap + 1;
      isInnerCell &= (int)globalCellIndex[currentDim] < upperBoundary[currentDim] - (int)_overlap - 1;
			isOuterCell |= (int)globalCellIndex[currentDim] < lowerBoundary[currentDim] + (int)_overlap - 1;
      isOuterCell |= (int)globalCellIndex[currentDim] > upperBoundary[currentDim] - (int)_overlap + 1;
		}
    return !isInnerCell && !isOuterCell;
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
};
} // namespace solvers
} // namespace coupling

#endif
