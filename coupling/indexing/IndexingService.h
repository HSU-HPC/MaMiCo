#pragma once

#include "coupling/configurations/MaMiCoConfiguration.h"
#include "coupling/interface/MacroscopicSolverInterface.h"
#include "simplemd/configurations/MolecularDynamicsConfiguration.h"
// parallel topologies
#include "coupling/paralleltopology/ParallelTopology.h"
#include "coupling/paralleltopology/ParallelTopologyFactory.h"

#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
#include <mpi.h>
#endif

namespace coupling {
namespace indexing {

template <unsigned int dim> class IndexingService;

}
} // namespace coupling

// Include CellIndex template class definition
#include "CellIndex.h"

namespace coupling {
namespace indexing {

template <unsigned int dim>
std::vector<unsigned int> getRanksForGlobalIndex(const BaseIndex<dim>& globalCellIndex,
                                                 const tarch::la::Vector<dim, unsigned int>& globalNumberMacroscopicCells);

} // namespace indexing
} // namespace coupling

// Include non-member functions operating on indexes
#include "Operations.h"

class IndexingServiceTest;

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

  void init(tarch::la::Vector<dim, unsigned int> globalNumberMacroscopicCells, tarch::la::Vector<dim, unsigned int> numberProcesses,
            coupling::paralleltopology::ParallelTopologyType type, unsigned int outerRegion, const unsigned int rank
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
            ,
            MPI_Comm comm = MPI_COMM_WORLD
#endif
  );

  // Config unpacking variant of init
  template <unsigned int mddim>
  typename std::enable_if<mddim == MD_DIM>::type init(const simplemd::configurations::MolecularDynamicsConfiguration& simpleMDConfig,
                                                      const coupling::configurations::MaMiCoConfiguration<mddim>& mamicoConfig,
                                                      coupling::interface::MacroscopicSolverInterface<mddim>* msi, const unsigned int rank
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
                                                      ,
                                                      MPI_Comm comm = MPI_COMM_WORLD
#endif
  ) {
    // read relevant data from configs
    const auto globalMDDomainSize{simpleMDConfig.getDomainConfiguration().getGlobalDomainSize()};
    const auto macroscopicCellSize{mamicoConfig.getMacroscopicCellConfiguration().getMacroscopicCellSize()};

    // calculate total number of macroscopic cells on all ranks in Base Domain
    tarch::la::Vector<dim, unsigned int> globalNumberMacroscopicCells(0);
    for (unsigned int d = 0; d < dim; d++) {
      globalNumberMacroscopicCells[d] = (unsigned int)floor(globalMDDomainSize[d] / macroscopicCellSize[d] + 0.5);

      if (fabs((globalNumberMacroscopicCells[d]) * macroscopicCellSize[d] - globalMDDomainSize[d]) > 1e-13)
        std::cout << "IndexingService: Deviation of domain size > 1e-13!" << std::endl;
    }

    init(globalNumberMacroscopicCells, simpleMDConfig.getMPIConfiguration().getNumberOfProcesses(),
         mamicoConfig.getParallelTopologyConfiguration().getParallelTopologyType(), msi, rank
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
         ,
         comm
#endif
    );
  }

  void init(tarch::la::Vector<dim, unsigned int> globalNumberMacroscopicCells, tarch::la::Vector<dim, unsigned int> numberProcesses,
            coupling::paralleltopology::ParallelTopologyType type, coupling::interface::MacroscopicSolverInterface<dim>* msi, const unsigned int rank
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
            ,
            MPI_Comm comm
#endif
  );

  void finalize() {
#if (COUPLING_MD_ERROR == COUPLING_MD_YES)
    _isInitialized = false;
#endif
  }

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

#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES) // parallel scenario
  MPI_Comm getComm() const {
#if (COUPLING_MD_ERROR == COUPLING_MD_YES)
    if (!_isInitialized) {
      throw std::runtime_error(std::string("IndexingService: Called index system getComm() before initalization! "));
    }
#endif

    return _comm;
  }
#endif

  unsigned int getRank() const {
#if (COUPLING_MD_ERROR == COUPLING_MD_YES)
    if (!_isInitialized) {
      throw std::runtime_error(std::string("IndexingService: Called index system getRank() before initalization! "));
    }
#endif

    return _rank;
  }

#if (COUPLING_MD_ERROR == COUPLING_MD_YES)
  bool isInitialized() const { return _isInitialized; }
#endif

private:
  /**
   * Helper function used by getRanksForGlobalIndex().
   */
  // TODO inline in getRanksForGlobalIndex()
  unsigned int getUniqueRankForMacroscopicCell(tarch::la::Vector<dim, unsigned int> globalCellIndex,
                                               const tarch::la::Vector<dim, unsigned int>& globalNumberMacroscopicCells) const;

  /*const*/ tarch::la::Vector<dim, unsigned int> _numberProcesses; // TODO: make const
  const coupling::paralleltopology::ParallelTopology<dim>* _parallelTopology;
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES) // parallel scenario
  MPI_Comm _comm;
#endif
  unsigned int _rank;
#if (COUPLING_MD_ERROR == COUPLING_MD_YES)
  bool _isInitialized = false;
#endif
  friend IndexingServiceTest;
};
