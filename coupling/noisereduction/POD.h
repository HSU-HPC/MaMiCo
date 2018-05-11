#ifndef _MOLECULARDYNAMICS_COUPLING_NOISEREDUCTION_POD_H_
#define _MOLECULARDYNAMICS_COUPLING_NOISEREDUCTION_POD_H_

#include "coupling/noisereduction/NoiseReduction.h"
//#define EIGEN_NO_DEBUG
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

namespace coupling {
  namespace noisereduction {
    template<class LinkedCell, unsigned int dim>
    class POD;
  }
}

/** Noise reduction algorithm using Proper orthogonal decomposition
 *
 *  @author Piet Jarmatz
 */
template<class LinkedCell,unsigned int dim>
class coupling::noisereduction::POD:
public coupling::noisereduction::NoiseReduction<LinkedCell,dim> {
  public:
    POD(const coupling::IndexConversion<dim> &indexConversion); // @todo get POD parameters from configuration
    virtual ~POD();

    virtual void processInnerMacroscopicCell(
      coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell,dim> &cell, const unsigned int &index
    );
    virtual void beginProcessInnerMacroscopicCells();
    virtual void endProcessInnerMacroscopicCells();
  private:
  	/** returns the local number of macroscopic cells EXCL. ghost layers */
    unsigned int getLocalNumberMacroscopicCells(const coupling::IndexConversion<dim> &indexConversion) const;

  	const unsigned int _timeWindowSize; // number of snapshots / coupling cycles taken into consideration for noise reduction
  	const unsigned int _kMax; // number of dominant eigenvalues
  	unsigned int _cycleCounter; // coupling cycle counter, indicates how many data snapshots are available already
  	unsigned int _spatialIndex; // cell counter, should run from zero to getLocalNumberMacroscopicCells()-1 within an iteration of ProcessInnerMacroscopicCells
  	unsigned int _t; // active temporal index, iterates cyclic between zero and _timeWindowSize
  	Eigen::MatrixXd *_data;  // set of snapshots (sampled by transferStrategy)
  	Eigen::MatrixXd *_C;  // temporal auto-correlation covariance matrix of _data
  	Eigen::MatrixXd *_A;  // POD temporal modes / eigenvectors of C
  	Eigen::MatrixXd *_A_T;  // Transpose of A
    bool _firstTraversal;   // distinguishes between reading and writing cell traversals
};
#include "coupling/noisereduction/POD.cpph"

#endif // _MOLECULARDYNAMICS_COUPLING_NOISEREDUCTION_POD_H_
