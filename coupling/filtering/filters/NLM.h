// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder

#pragma once

#define NLM_DEBUG

#include "coupling/filtering/filters/Datastructures.h"
#include "coupling/filtering/interfaces/JunctorInterface.h"
#include <array>

namespace coupling {
namespace filtering {
template <unsigned int dim> class NLM;
}
} // namespace coupling

/**
 *  Noise reduction algorithm using non-local means (NLM) method
 *  See 'Fast Non Local Means Denoising for 3D MR Images' by Coupé et al. 2006
 *   and 'Non-Local Means Denoising' by Buades et al. 2011.
 *
 *  @author Piet Jarmatz
 *
 */
template <unsigned int dim>
class coupling::filtering::NLM
    : public coupling::filtering::JunctorInterface<dim, 2, 1> {
public:
  using CellIndex_T =
      coupling::indexing::CellIndex<dim, coupling::indexing::IndexTrait::vector,
                                    coupling::indexing::IndexTrait::local,
                                    coupling::indexing::IndexTrait::md2macro,
                                    coupling::indexing::IndexTrait::noGhost>;

  NLM(const std::vector<coupling::datastructures::MacroscopicCell<dim> *>
          inputCellVector_unfiltered,
      const std::vector<coupling::datastructures::MacroscopicCell<dim> *>
          inputCellVector_prefiltered,
      const std::vector<coupling::datastructures::MacroscopicCell<dim> *>
          outputCellVector,
      const std::array<bool, 7> filteredValues, int tws, double sigsq,
      double sigsq_rel, double hsq, double hsq_rel, int M = 2, int d = 1)
      : coupling::filtering::JunctorInterface<dim, 2, 1>(
            {inputCellVector_unfiltered, inputCellVector_prefiltered},
            {outputCellVector}, filteredValues, "NLM"),
        _timeWindowSize(tws), _M(M), _d(d), _sigsq(sigsq),
        _sigsq_rel(sigsq_rel), _hsq(hsq), _hsq_rel(hsq_rel), _cycleCounter(0),
        _t(0), _flowfield(CellIndex_T::numberCellsInDomain, tws),
        _flowfield_prefiltered(CellIndex_T::numberCellsInDomain, tws),
        _patchfield(CellIndex_T::numberCellsInDomain -
                        tarch::la::Vector<dim, unsigned int>(2),
                    tws),
        _innerCellIndices() {
    // Initialize flowfield
    for (int t = 0; t < tws; ++t)
      for (auto idx : CellIndex_T()) {
        tarch::la::Vector<dim, unsigned int> idxv(idx.get());
        coupling::filtering::Quantities<dim> &q = _flowfield(idxv, t);
        q[0] = 1;
        for (unsigned int d = 1; d <= dim; ++d)
          q[d] = 0;
        coupling::filtering::Quantities<dim> &qp =
            _flowfield_prefiltered(idxv, t);
        qp[0] = 1;
        for (unsigned int d = 1; d <= dim; ++d)
          qp[d] = 0;
      }

    // Initialize innerCellIndices
    auto domainSize = CellIndex_T::numberCellsInDomain;
    for (auto idx : CellIndex_T()) {
      tarch::la::Vector<dim, unsigned int> idxv(idx.get());
      for (unsigned int d = 0; d < dim; d++)
        if (idxv[d] == 0 || idxv[d] == domainSize[d] - 1)
          goto continue_loop;
      _innerCellIndices.push_back(idx);
    continue_loop:;
    }

#ifdef NLM_DEBUG
    std::cout << "    NLM: Created NLM instance." << std::endl;
    std::cout << "    WARNING: Regardless of configuration, NLM always filters "
                 "macroscopic mass and momentum."
              << std::endl;
#endif
  }

  virtual ~NLM() {
#ifdef NLM_DEBUG
    std::cout << "    NLM: Destroyed NLM instance." << std::endl;
#endif
  }

  void operator()();

private:
  /**
   * Stores new values from inputCellVectors into _flowfield and
   * _flowfield_prefiltered for this timestep _t.
   */
  void save_data();

  /**
   *  Reconstructs _patchfield(idx,_t) for all inner cell indices idx,
   *  so that new information from updated flowfield will be included.
   */
  void update_patchfield();

  /**
   * Main loop of NLM algorithm, computes filtering output for one coupling
   * cycle.
   */
  void denoise();

  /**
   * Compute which filtering results we need here.
   * @param idx_inner Index of this cell on the inner local MD 2 macro domain,
   * i.e. the domain where the patchfield is defined
   * @return Relative offsets, of all neighboring boundary cells and this cell
   * itself
   */
  std::vector<tarch::la::Vector<dim, int>>
  compute_boundary_neighbors(const tarch::la::Vector<dim, unsigned int> &);

  /**
   * Must be called to indicate that one coupling cycle is finished and timestep
   * counters should be incremented.
   */
  void increment_time();

  const unsigned int
      _timeWindowSize;   // number of snapshots / coupling cycles taken into
                         // consideration for noise reduction
                         // unsigned int _timeModulo; // e.g. 2,4,8 ....
  const unsigned int _M; // search volume has size (2M+1)^3
  const unsigned int _d; // patches have size (2d+1)^4; this makes at least _d
                         // ghost layer necessary
  double _sigsq, _sigsq_rel;
  double _hsq, _hsq_rel;

  unsigned int _cycleCounter; // coupling cycle counter, indicates how many data
                              // snapshots are available already
  unsigned int _t; // active temporal index, iterates cyclic between zero and
                   // _timeWindowSize
  coupling::filtering::Flowfield<dim> _flowfield;
  coupling::filtering::Flowfield<dim> _flowfield_prefiltered;
  coupling::filtering::Patchfield<dim> _patchfield;
  std::vector<CellIndex_T> _innerCellIndices;

  inline unsigned int posmod(int i, int n) const { return (i % n + n) % n; }
};

// include implementation of header
#include "NLM.cpph"
