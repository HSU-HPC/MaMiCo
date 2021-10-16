// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder

#pragma once

#define NLM_DEBUG

#include "coupling/filtering/interfaces/JunctorInterface.h"
#include "coupling/filtering/filters/Datastructures.h"
#include <array>

namespace coupling {
  template<unsigned int dim>
  class NLM;
}

/** 
 *  Noise reduction algorithm using non-local means (NLM) method
 *  See 'Fast Non Local Means Denoising for 3D MR Images' by Coupé et al. 2006
 *   and 'Non-Local Means Denoising' by Buades et al. 2011. 
 *
 *  @author Piet Jarmatz
 * 
 */
template<unsigned int dim>
class coupling::NLM : public coupling::JunctorInterface<dim,2,1> {
public:
  NLM(const std::vector<coupling::datastructures::MacroscopicCell<dim> *> inputCellVector_unfiltered,
    const std::vector<coupling::datastructures::MacroscopicCell<dim> *> inputCellVector_prefiltered,
    const std::vector<coupling::datastructures::MacroscopicCell<dim> *> outputCellVector,
    const std::vector<tarch::la::Vector<dim, unsigned int>> cellIndices,
    const std::array<bool, 7> filteredValues, 
    const coupling::IndexConversionMD2Macro<dim>* indexConversion,
    int tws,
    double sigsq,
    double hsq,
    int M = 2,
    int d = 1
    ):
  coupling::JunctorInterface<dim,2,1>( 
   { inputCellVector_unfiltered, inputCellVector_prefiltered },
   { outputCellVector }, 
   cellIndices, 
   filteredValues,
   "NLM"),
  _timeWindowSize(tws),
  _M(M),
  _d(d),
  _sigsq(sigsq),
  _hsq(hsq),
  _cycleCounter(0),
  _t(0),
  _ic(indexConversion),
  _flowfield(_ic->getLocalMD2MacroDomainSize(), tws),
  _flowfield_prefiltered(_ic->getLocalMD2MacroDomainSize(), tws),
  _patchfield(_ic->getLocalMD2MacroDomainSize() - tarch::la::Vector<dim,unsigned int>(2), tws),
  _innerCellIndices()
  {
    // Initialize flowfield
    for(int t = 0; t < tws; ++t)
      for(const tarch::la::Vector<dim, unsigned int>& idx : cellIndices) {
        coupling::filtering::Quantities<dim>& q = _flowfield(idx, t);
        q[0] = 1;
        for(unsigned int d = 1; d <= dim; ++d) q[d] = 0;
          coupling::filtering::Quantities<dim>& qp = _flowfield_prefiltered(idx, t);
        qp[0] = 1;
        for(unsigned int d = 1; d <= dim; ++d) qp[d] = 0;
      }

    // Initialize innerCellIndices
    auto domainSize = _ic->getLocalMD2MacroDomainSize();
    for(auto idx : cellIndices){ 
      for (unsigned int d = 0; d < dim; d++)  
        if(idx[d] == 0 || idx[d] == domainSize[d] - 1)
          goto continue_loop;
      _innerCellIndices.push_back(idx);
      continue_loop:;
    }

    #ifdef NLM_DEBUG
    std::cout << "    NLM: Created NLM instance." << std::endl;
    std::cout << "    WARNING: Regardless of configuration, NLM always filters macroscopic mass and momentum." << std::endl;
    #endif
  }

  virtual ~NLM(){
    #ifdef NLM_DEBUG
    std::cout << "    NLM: Destroyed NLM instance." << std::endl;
    #endif
  }

  void operator()();
private:
  /** 
   * Stores new values from inputCellVectors into _flowfield and _flowfield_prefiltered for this timestep _t.
   */
  void save_data();

  /**
   *  Reconstructs _patchfield(idx,_t) for all inner cell indices idx, 
   *  so that new information from updated flowfield will be included.
   */
  void update_patchfield();

  /**
   * Main loop of NLM algorithm, computes filtering output for one coupling cycle. 
   */
  void denoise();

  /**
   * Compute which filtering results we need here. 
   * @param idx_inner Index of this cell on the inner local MD 2 macro domain, i.e. the domain where the patchfield is defined
   * @return Relative offsets, of all neighboring boundary cells and this cell itself
   */
  std::vector<tarch::la::Vector<dim, int>> compute_boundary_neighbors(const tarch::la::Vector<dim, unsigned int>&);

  /**
   * Must be called to indicate that one coupling cycle is finished and timestep counters should be incremented.
   */
  void increment_time();

	const unsigned int _timeWindowSize; // number of snapshots / coupling cycles taken into consideration for noise reduction
  //unsigned int _timeModulo; // e.g. 2,4,8 ....
	const unsigned int _M; // search volume has size (2M+1)^3
	const unsigned int _d; // patches have size (2d+1)^4; this makes at least _d ghost layer necessary
  const double _sigsq;
  const double _hsq;

  unsigned int _cycleCounter; // coupling cycle counter, indicates how many data snapshots are available already
  unsigned int _t; // active temporal index, iterates cyclic between zero and _timeWindowSize
  const coupling::IndexConversionMD2Macro<dim>* _ic;
  coupling::filtering::Flowfield<dim> _flowfield;
  coupling::filtering::Flowfield<dim> _flowfield_prefiltered;
  coupling::filtering::Patchfield<dim> _patchfield;
  std::vector<tarch::la::Vector<dim, unsigned int>> _innerCellIndices;

  inline unsigned int posmod(int i, int n) const{
    return (i % n + n) % n;
  }
};

//include implementation of header
#include "NLM.cpph"
