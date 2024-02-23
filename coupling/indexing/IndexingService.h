#pragma once

// parallel topologies
#include "coupling/CouplingMDDefinitions.h"
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
#include "coupling/indexing/IndexTypes.h"

namespace coupling {
namespace indexing {

template <unsigned int dim>
std::vector<unsigned int> getRanksForGlobalIndex(const BaseIndex<dim>& globalCellIndex, const tarch::la::Vector<dim, unsigned int>& globalNumberCouplingCells);

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
 * @param outer region
 *
 * @author Felix Maurer
 */
template <unsigned int dim> class coupling::indexing::IndexingService {
public:
  static IndexingService& getInstance() {
    static IndexingService singleton{};
    return singleton;
  }

  void initWithCells(tarch::la::Vector<dim, unsigned int> globalNumberCouplingCells, tarch::la::Vector<dim, unsigned int> numberProcesses,
                     const tarch::la::Vector<3, double>& couplingCellSize, coupling::paralleltopology::ParallelTopologyType type, unsigned int outerRegion,
                     const unsigned int rank
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
                     ,
                     MPI_Comm comm = MPI_COMM_WORLD
#endif
  );

  void initWithMDSize(const tarch::la::Vector<3, double>& globalMDDomainSize, const tarch::la::Vector<3, double>& globalMDDomainOffset,
                      const tarch::la::Vector<3, unsigned int>& mdNumberProcesses, const tarch::la::Vector<3, double>& couplingCellSize,
                      coupling::paralleltopology::ParallelTopologyType parallelTopologyType, unsigned int outerRegion, unsigned int rank
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
                      ,
                      MPI_Comm comm = MPI_COMM_WORLD
#endif
  ) {

    _globalMDDomainSize = globalMDDomainSize;
    _globalMDDomainOffset = globalMDDomainOffset;

    // calculate total number of coupling cells on all ranks in Base Domain
    tarch::la::Vector<dim, unsigned int> globalNumberCouplingCells(0);
    for (unsigned int d = 0; d < dim; d++) {
      globalNumberCouplingCells[d] = (unsigned int)floor(globalMDDomainSize[d] / couplingCellSize[d] + 0.5);

      if (fabs((globalNumberCouplingCells[d]) * couplingCellSize[d] - globalMDDomainSize[d]) > 1e-13)
        std::cout << "IndexingService: Deviation of domain size > 1e-13!" << std::endl;
    }

    initWithCells(globalNumberCouplingCells, mdNumberProcesses, couplingCellSize, parallelTopologyType, outerRegion, rank
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
                  ,
                  comm
#endif
    );
  }

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
   * @param globalNumberCouplingCells global number of cells in BaseIndex
   * domain EXCLUDING global ghost layer cells.
   * @returns vector of all cells which contain the index
   */
  std::vector<unsigned int> getRanksForGlobalIndex(const BaseIndex<dim>& globalCellIndex) const;

  unsigned int getUniqueRankForCouplingCell(const BaseIndex<dim>& globalCellIndex) const;

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

  /** returns the global domain size of the MD domain (excl. ghost layer which
   * naturally is not part of MD).
   *  @returns the total size of the md simulation domain (dimensional) */
  tarch::la::Vector<dim, double> getGlobalMDDomainSize() const {
#if (COUPLING_MD_ERROR == COUPLING_MD_YES)
    if (!_isInitialized) {
      throw std::runtime_error(std::string("IndexingService: Called getGlobalMDDomainSize() before initalization! "));
    }
#endif
    return _globalMDDomainSize;
  }

  /** @brief returns the offset, i.e. the lower,left... corner coordinate, of
   * the MD domain.
   *  @returns the offset of the MD domain */
  tarch::la::Vector<dim, double> getGlobalMDDomainOffset() const {
#if (COUPLING_MD_ERROR == COUPLING_MD_YES)
    if (!_isInitialized) {
      throw std::runtime_error(std::string("IndexingService: Called getGlobalMDDomainOffset() before initalization! "));
    }
#endif
    return _globalMDDomainOffset;
  }

  /** @brief returns the vector size of each coupling cell.
   *  @returns the size of the coupling cells (dimensional) */
  tarch::la::Vector<dim, double> getCouplingCellSize() const {
#if (COUPLING_MD_ERROR == COUPLING_MD_YES)
    if (!_isInitialized) {
      throw std::runtime_error(std::string("IndexingService: Called getCouplingCellSize() before initalization! "));
    }
#endif
    return _couplingCellSize;
  }

  unsigned int getScalarNumberProcesses() const {
#if (COUPLING_MD_ERROR == COUPLING_MD_YES)
    if (!_isInitialized) {
      throw std::runtime_error(std::string("IndexingService: Called getScalarNumberProcesses() before initalization! "));
    }
#endif
    return _scalarNumberProcesses;
  }

private:
  unsigned int getUniqueRankForCouplingCell(tarch::la::Vector<dim, unsigned int> globalCellIndex,
                                            const tarch::la::Vector<dim, unsigned int>& globalNumberCouplingCells) const;

  /*const*/ tarch::la::Vector<dim, unsigned int> _numberProcesses; // TODO: make const
  unsigned int _scalarNumberProcesses;
  const coupling::paralleltopology::ParallelTopology<dim>* _parallelTopology;
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES) // parallel scenario
  MPI_Comm _comm;
#endif
  unsigned int _rank;
#if (COUPLING_MD_ERROR == COUPLING_MD_YES)
  bool _isInitialized = false;
#endif
  tarch::la::Vector<dim, double> _globalMDDomainSize;
  tarch::la::Vector<dim, double> _globalMDDomainOffset;
  tarch::la::Vector<dim, double> _couplingCellSize;
  friend IndexingServiceTest;
};
