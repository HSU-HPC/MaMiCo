#pragma once

#include "coupling/configurations/MaMiCoConfiguration.h"
#include "coupling/interface/MacroscopicSolverInterface.h"
#include "simplemd/configurations/MolecularDynamicsConfiguration.h"
// parallel topologies
#include "coupling/paralleltopology/ParallelTopology.h"
#include "coupling/paralleltopology/ParallelTopologyFactory.h"

// Include CellIndex template class definition
#include "CellIndex.h"

namespace coupling {
namespace indexing {

template <unsigned int dim> class IndexingService;

// template <unsigned int dim>
// std::vector<unsigned int> getRanksForGlobalIndex(const BaseIndex<dim>& globalCellIndex,
//                                                  const tarch::la::Vector<dim, unsigned int>& globalNumberMacroscopicCells);

} // namespace indexing
} // namespace coupling

// Include non-member functions operating on indexes
#include "Operations.h"

// enable/disable tests
//#define TEST_INDEXING

#ifdef TEST_INDEXING
// Inlcude index tests
#include "Testing.h"
#endif

/**
 * Singleton service class initialising lower and upper boundaries of all
 * possible CellIndex specialisations.
 *
 * @tparam dim number of dimensions of the coupled simulation
 * @param simpleMDConfig config object of SimpleMD instance used in coupling
 * @param mamicoConfig config object containg general information of coupling
 * process
 * @param msi pointer to interface of coupled macroscopic solver
 *
 * @author Felix Maurer
 */
template <unsigned int dim> class coupling::indexing::IndexingService {
public:
  static IndexingService& getInstance() {
    static IndexingService singleton{};
    return singleton;
  }

  void init(const simplemd::configurations::MolecularDynamicsConfiguration& simpleMDConfig,
            const coupling::configurations::MaMiCoConfiguration<dim>& mamicoConfig, coupling::interface::MacroscopicSolverInterface<dim>* msi,
            const unsigned int rank);

#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES) // parallel scenario
  /**
   * Determines all ranks that contain a certain global BaseIndex.
   * Ripped from deprecated IndexConversion.
   *
   * @param globalCellIndex index to be looked up
   * @param globalNumberMacroscopicCells global number of cells in BaseIndex
   * domain EXCLUDING global ghost layer cells.
   * @returns vector of all cells which contain the index
   */
  std::vector<unsigned int> getRanksForGlobalIndex(const BaseIndex<dim>& globalCellIndex) const;
#endif

  unsigned int getRank() const { return _rank; }

#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES) // parallel scenario
   /** If the cell is contained in the MD volume, the rank is uniquely chosen by
   * using the process which contains the macroscopic cell as non-ghost (i.e.
   * real inner) cell. If the cell is a global ghost cell, we choose the rank
   * according to the block-decomposition of the grid applied to the complete
   * grid incl. the ghost layer.
   *  @brief returns the unique rank for a macroscopic cell.
   *  @param globalCellIndex the global vector coordinates of a macroscopic cell
   *  @returns the rank for the given macroscopic cell */
  // TODO inline in getRanksForGlobalIndex()
  unsigned int getUniqueRankForGlobalIndex(tarch::la::Vector<dim, unsigned int> globalCellIndex) const;
#endif

private:
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES) // parallel scenario
  /*const*/ tarch::la::Vector<dim, unsigned int> _numberProcesses; // TODO: make const
  const coupling::paralleltopology::ParallelTopology<dim>* _parallelTopology;
#endif

  simplemd::configurations::MolecularDynamicsConfiguration _simpleMDConfig;
  coupling::configurations::MaMiCoConfiguration<dim> _mamicoConfig;
  coupling::interface::MacroscopicSolverInterface<dim>* _msi;
  unsigned int _rank;
};
