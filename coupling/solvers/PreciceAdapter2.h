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

template <unsigned int dim> class PreciceAdapter {
public:
  PreciceAdapter()
  {
  }
  virtual ~PreciceAdapter()
  {
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

  void setCouplingMesh() {
    int rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    _interface = new precice::SolverInterface("mamico", "../precice-config.xml", rank, size);

    _numberOfM2mCells = 3;
    std::cout << "MaMiCo >> mamico-M2m-mesh size: " << _numberOfM2mCells << std::endl;
    _coordsM2mCells = new double[] {100, 100, 6.25, 100, 100, 8.75, 100, 100, 11.25};
    _vertexM2mCellIDs = new int[_numberOfM2mCells];
    _interface->setMeshVertices(_interface->getMeshID("mamico-M2m-mesh"), _numberOfM2mCells, _coordsM2mCells, _vertexM2mCellIDs);
    _velocityM2mCells = new double[_numberOfM2mCells * dim];
    _precice_dt = _interface->initialize();
  }

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
    std::cout << "pos:" << pos << std::endl;
    std::cout << "vel:" << vel << std::endl;
    return vel;
  }

  double getDensity(tarch::la::Vector<3, double> pos) const { return 1.0; }

  void advance(double dt) {
    if (_interface->isCouplingOngoing()) {
      if (_interface->isReadDataAvailable()) {
        std::cout << "MaMiCo >> Reading macro velocities from preCICE !" << std::endl;
        // velocity from the conitnuum solver
        _interface->readBlockVectorData(_interface->getDataID("VelocityMacro", _interface->getMeshID("mamico-M2m-mesh")), _numberOfM2mCells, _vertexM2mCellIDs,
                                        _velocityM2mCells);
      }

      // Solving the time step
      // Normally does nothing, everything is done on the MD side

      double computed_dt = std::min(_precice_dt, dt);
      _precice_dt = _interface->advance(computed_dt);
    }
  }


private:
  precice::SolverInterface* _interface = NULL;
  double _precice_dt;

  int* _vertexM2mCellIDs;
  double* _coordsM2mCells;
  unsigned int _numberOfM2mCells;
  double* _velocityM2mCells;

};
} // namespace solvers
} // namespace coupling

#endif
