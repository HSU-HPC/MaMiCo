#ifndef _MOLECULARDYNAMICS_COUPLING_NOISEREDUCTION_NLM_H_
#define _MOLECULARDYNAMICS_COUPLING_NOISEREDUCTION_NLM_H_

#include "coupling/noisereduction/NoiseReduction.h"

namespace coupling {
  namespace noisereduction {
    template<unsigned int dim>
    class NLM;

    template<unsigned int dim>
    class Cell;
  }
}

/** Noise reduction algorithm using non-local means (NLM) method
 *  See 'Fast Non Local Means Denoising for 3D MR Images' by Coupé et al. 2006.
 *
 *  @author Piet Jarmatz
 */
template<unsigned int dim>
class coupling::noisereduction::NLM:
public coupling::noisereduction::NoiseReduction<dim> {
  public:
    NLM(const coupling::IndexConversion<dim> &indexConversion, const tarch::utils::MultiMDService<dim>& multiMDService, int tws);
    virtual ~NLM();

    virtual void processInnerMacroscopicCell(
      coupling::datastructures::MacroscopicCell<dim> &cell, const unsigned int &index
    );
    virtual void beginProcessInnerMacroscopicCells();
    virtual void endProcessInnerMacroscopicCells();

    virtual void processOuterMacroscopicCell(
      coupling::datastructures::MacroscopicCell<dim> &cell, const unsigned int &index
    );
    virtual void beginProcessOuterMacroscopicCells();
    virtual void endProcessOuterMacroscopicCells();
  private:
  	const unsigned int _timeWindowSize; // number of snapshots / coupling cycles taken into consideration for noise reduction
    //unsigned int _timeModulo; // e.g. 2,4,8 ....
  	//const unsigned int _M; // search volume has size (2M+1)^4
  	const unsigned int _d; // local neighborhoods have size (2d+1)^4; this makes at least _d ghost layer necessary
    unsigned int _cycleCounter; // coupling cycle counter, indicates how many data snapshots are available already
    unsigned int _t; // active temporal index, iterates cyclic between zero and _timeWindowSize
    bool _firstTraversal;   // distinguishes between reading and writing cell traversals
    unsigned int _spatialSize;
    coupling::noisereduction::Cell<dim>* _data; 
};

// datastructure to store information about one filter cell:
// => mass / momentum from the microscopic for the macroscopic simulation
// => mean and standard deviation in a local neighborhood of size (2d+1)^4 centered at this cell
template<unsigned int dim>
class coupling::noisereduction::Cell{
  public:
    Cell(): _mass(0), _momentum(0), 
    _localMeanMass(0), _localMeanMomentum(0), 
    _localStandardDeviationMass(0), _localStandardDeviationMomentum(0) {}

    double _mass;
    tarch::la::Vector<dim,double> _momentum;
    double _localMeanMass;
    tarch::la::Vector<dim,double> _localMeanMomentum;
    double _localStandardDeviationMass;
    tarch::la::Vector<dim,double> _localStandardDeviationMomentum;
};

#include "coupling/noisereduction/NLM.cpph"

#endif // _MOLECULARDYNAMICS_COUPLING_NOISEREDUCTION_NLM_H_
