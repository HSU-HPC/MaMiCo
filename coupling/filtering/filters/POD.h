// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#pragma once
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <string>
#include <vector>

//#define DEBUG_POD
#include "coupling/filtering/interfaces/FilterInterface.h"

namespace coupling {
namespace filtering {
template <unsigned int dim> class POD;
}
} // namespace coupling

/** Noise reduction algorithm using Proper orthogonal decomposition
 *
 *  @author Piet Jarmatz, Felix Maurer
 */
template <unsigned int dim> class coupling::filtering::POD : public coupling::filtering::FilterInterface<dim> {
public:
  POD(const std::vector<coupling::datastructures::MacroscopicCell<dim>*>& inputCellVector,
      const std::vector<coupling::datastructures::MacroscopicCell<dim>*>& outputCellVector,
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
      MPI_Comm comm,
#endif
      const std::array<bool, 7> filteredValues, int tws, int kmax)
      : coupling::filtering::FilterInterface<dim>(inputCellVector, outputCellVector, filteredValues, "POD"), _timeWindowSize(tws), _kMax(kmax),
        _cycleCounter(0), _spatialIndex(0), _t(0), _data(NULL), _C(NULL), _A(NULL), _A_T(NULL) {
    int spatialSize = inputCellVector.size();
    _data = new Eigen::MatrixXd[dim + 1]; // separate data matrices for: mass,
                                          // momentum0, momentum1, momentum2
    _C = new Eigen::MatrixXd[dim + 1];
    _A = new Eigen::MatrixXd[dim + 1];
    _A_T = new Eigen::MatrixXd[dim + 1];
    for (unsigned int i = 0; i < dim + 1; i++) {
      _data[i] = Eigen::MatrixXd::Constant(_timeWindowSize, spatialSize, (i == 0) ? 1 : 0);
    }

#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
    // get MPI parameters
    _comm = comm;
    MPI_Comm_rank(comm, (int*)&_rank);
    MPI_Comm_size(comm, (int*)&_commSize);
#endif

#ifdef DEBUG_POD
    std::cout << "		POD: Created Proper Orthogonal Decomposition instance." << std::endl;
    // TODO selection of filtered properties
    std::cout << "			WARNING: Regardless of configuration, "
                 "POD always filters macroscopic mass and momentum."
              << std::endl;
#endif
  }

  ~POD() {
    if (_data != NULL) {
      delete[] _data;
      _data = NULL;
    }
    if (_C != NULL) {
      delete[] _C;
      _C = NULL;
    }
    if (_A != NULL) {
      delete[] _A;
      _A = NULL;
    }
    if (_A_T != NULL) {
      delete[] _A_T;
      _A_T = NULL;
    }

#ifdef DEBUG_POD
    std::cout << "		POD: Deleted Proper Orthogonal Decomposition instance." << std::endl;
#endif
  }

  void operator()();

private:
  unsigned int _timeWindowSize; // number of snapshots / coupling cycles taken
                                // into consideration for noise reduction
  const unsigned int _kMax;     // number of dominant eigenvalues
  unsigned int _cycleCounter;   // coupling cycle counter, indicates how many data
                                // snapshots are available already
  unsigned int _spatialIndex;   // cell counter, should run from zero to
                                // getLocalNumberMacroscopicCells()-1 within an
                                // iteration of ProcessInnerMacroscopicCell
  unsigned int _t;              // active temporal index, iterates cyclic between zero and
                                // _timeWindowSize
  Eigen::MatrixXd* _data;       // set of snapshots (sampled by transferStrategy)
  Eigen::MatrixXd* _C;          // temporal auto-correlation covariance matrix of _data
  Eigen::MatrixXd* _A;          // POD temporal modes / eigenvectors of C
  Eigen::MatrixXd* _A_T;        // Transpose of A

#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
  MPI_Comm _comm;
  unsigned int _rank;
  unsigned int _commSize;
#endif
};

// include implementation of header
#include "POD.cpph"
