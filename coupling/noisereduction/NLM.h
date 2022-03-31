#ifndef _MOLECULARDYNAMICS_COUPLING_NOISEREDUCTION_NLM_H_
#define _MOLECULARDYNAMICS_COUPLING_NOISEREDUCTION_NLM_H_

#define NLM_DEBUG

#include "coupling/noisereduction/NoiseReduction.h"
#include "coupling/noisereduction/Datastructures.h"

namespace coupling {
namespace noisereduction { template <unsigned int dim> class NLM; }
}

/** Noise reduction algorithm using non-local means (NLM) method
 *  See 'Fast Non Local Means Denoising for 3D MR Images' by Coupé et al. 2006.
 *
 *  @author Piet Jarmatz
 */
template <unsigned int dim>
class coupling::noisereduction::NLM
    : public coupling::noisereduction::NoiseReduction<dim> {
public:
  NLM(const coupling::IndexConversion<dim> &indexConversion,
      const tarch::utils::MultiMDService<dim> &multiMDService, int tws);
  virtual ~NLM();

  virtual void processInnerMacroscopicCell(
      coupling::datastructures::MacroscopicCell<dim> &cell,
      const unsigned int &index);
  virtual void beginProcessInnerMacroscopicCells();
  virtual void endProcessInnerMacroscopicCells();

  virtual void processOuterMacroscopicCell(
      coupling::datastructures::MacroscopicCell<dim> &cell,
      const unsigned int &index);
  virtual void beginProcessOuterMacroscopicCells();
  virtual void endProcessOuterMacroscopicCells();

private:
  const unsigned int
      _timeWindowSize; // number of snapshots / coupling cycles taken into
                       // consideration for noise reduction
  //unsigned int _timeModulo; // e.g. 2,4,8 ....
  const unsigned int _M; // search volume has size (2M+1)^3
  const unsigned int _d; // patches have size (2d+1)^4; this makes at least _d
                         // ghost layer necessary
  unsigned int _cycleCounter; // coupling cycle counter, indicates how many data
                              // snapshots are available already
  unsigned int _t; // active temporal index, iterates cyclic between zero and
                   // _timeWindowSize
  bool _firstTraversal; // distinguishes between reading and writing cell
                        // traversals
  coupling::noisereduction::Flowfield<dim> _flowfield;
  coupling::noisereduction::Patchfield<dim> _patchfield;
};

#include "coupling/noisereduction/NLM.cpph"

#endif // _MOLECULARDYNAMICS_COUPLING_NOISEREDUCTION_NLM_H_
