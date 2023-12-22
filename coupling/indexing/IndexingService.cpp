// Include header
#include "IndexingService.h"
#include "coupling/solvers/CouetteSolverInterface.h" // to create default msi

#include <algorithm>
#include <iterator>

/*
 * Define specialisations of CellIndex.
 * Differentiate between dim=2 and dim=3.
 *
 * Only declarations and default inits should be made here, lower/upper
 * boundaries must be determined at runtime using IndexingService
 */

// #DEFINE INDEXING_ENABLE_DIM2
#ifdef INDEXING_ENABLE_DIM2
// We must compile both, dim2 and dim3, so that both can be used (also in the same executable e.g. main_lammps.cpp)
// Dim2    ///////////
namespace coupling {
namespace indexing {

/*
 * Declare specialisations of CellIndex.
 * Define their static members.
 */

/*
 * NON-MD-TO-MACRO aka MAMICO INDEXING, INCL GHOST LAYER
 */

// scalar, global, !md2macro, !noGL
template <> BaseIndex<2> CellIndex<2>::lowerBoundary{};
template <> BaseIndex<2> CellIndex<2>::upperBoundary{};
template <> tarch::la::Vector<2, unsigned int> CellIndex<2>::numberCellsInDomain{};
template <> unsigned int CellIndex<2>::linearNumberCellsInDomain{};
template <> tarch::la::Vector<2, unsigned int> CellIndex<2>::divisionFactor{};

// BaseIndex
template <> BaseIndex<2> BaseIndex<2>::lowerBoundary{};
template <> BaseIndex<2> BaseIndex<2>::upperBoundary{};
template <> tarch::la::Vector<2, unsigned int> BaseIndex<2>::numberCellsInDomain{};
template <> unsigned int BaseIndex<2>::linearNumberCellsInDomain{};
template <> tarch::la::Vector<2, unsigned int> BaseIndex<2>::divisionFactor{};

// scalar, local, !md2macro, !noGL
template <> BaseIndex<2> CellIndex<2, IndexTrait::local>::lowerBoundary{};
template <> BaseIndex<2> CellIndex<2, IndexTrait::local>::upperBoundary{};
template <> tarch::la::Vector<2, unsigned int> CellIndex<2, IndexTrait::local>::numberCellsInDomain{};
template <> unsigned int CellIndex<2, IndexTrait::local>::linearNumberCellsInDomain{};
template <> tarch::la::Vector<2, unsigned int> CellIndex<2, IndexTrait::local>::divisionFactor{};

// vector, local, !md2macro, !noGL
template <> BaseIndex<2> CellIndex<2, IndexTrait::vector, IndexTrait::local>::lowerBoundary{};
template <> BaseIndex<2> CellIndex<2, IndexTrait::vector, IndexTrait::local>::upperBoundary{};
template <> tarch::la::Vector<2, unsigned int> CellIndex<2, IndexTrait::vector, IndexTrait::local>::numberCellsInDomain{};
template <> unsigned int CellIndex<2, IndexTrait::vector, IndexTrait::local>::linearNumberCellsInDomain{};
template <> tarch::la::Vector<2, unsigned int> CellIndex<2, IndexTrait::vector, IndexTrait::local>::divisionFactor{};

/*
 * MD TO MACRO, INCL GHOST LAYER
 */

// scalar, global, md2macro, !noGL
template <> BaseIndex<2> CellIndex<2, IndexTrait::md2macro>::lowerBoundary{};
template <> BaseIndex<2> CellIndex<2, IndexTrait::md2macro>::upperBoundary{};
template <> tarch::la::Vector<2, unsigned int> CellIndex<2, IndexTrait::md2macro>::numberCellsInDomain{};
template <> unsigned int CellIndex<2, IndexTrait::md2macro>::linearNumberCellsInDomain{};
template <> tarch::la::Vector<2, unsigned int> CellIndex<2, IndexTrait::md2macro>::divisionFactor{};

// vector, global, md2macro, !noGL
template <> BaseIndex<2> CellIndex<2, IndexTrait::vector, IndexTrait::md2macro>::lowerBoundary{};
template <> BaseIndex<2> CellIndex<2, IndexTrait::vector, IndexTrait::md2macro>::upperBoundary{};
template <> tarch::la::Vector<2, unsigned int> CellIndex<2, IndexTrait::vector, IndexTrait::md2macro>::numberCellsInDomain{};
template <> unsigned int CellIndex<2, IndexTrait::vector, IndexTrait::md2macro>::linearNumberCellsInDomain{};
template <> tarch::la::Vector<2, unsigned int> CellIndex<2, IndexTrait::vector, IndexTrait::md2macro>::divisionFactor{};

// scalar, local, md2macro, !noGL
template <> BaseIndex<2> CellIndex<2, IndexTrait::local, IndexTrait::md2macro>::lowerBoundary{};
template <> BaseIndex<2> CellIndex<2, IndexTrait::local, IndexTrait::md2macro>::upperBoundary{};
template <> tarch::la::Vector<2, unsigned int> CellIndex<2, IndexTrait::local, IndexTrait::md2macro>::numberCellsInDomain{};
template <> unsigned int CellIndex<2, IndexTrait::local, IndexTrait::md2macro>::linearNumberCellsInDomain{};
template <> tarch::la::Vector<2, unsigned int> CellIndex<2, IndexTrait::local, IndexTrait::md2macro>::divisionFactor{};

// vector, local, md2macro, !noGL
template <> BaseIndex<2> CellIndex<2, IndexTrait::vector, IndexTrait::local, IndexTrait::md2macro>::lowerBoundary{};
template <> BaseIndex<2> CellIndex<2, IndexTrait::vector, IndexTrait::local, IndexTrait::md2macro>::upperBoundary{};
template <> tarch::la::Vector<2, unsigned int> CellIndex<2, IndexTrait::vector, IndexTrait::local, IndexTrait::md2macro>::numberCellsInDomain{};
template <> unsigned int CellIndex<2, IndexTrait::vector, IndexTrait::local, IndexTrait::md2macro>::linearNumberCellsInDomain{};
template <> tarch::la::Vector<2, unsigned int> CellIndex<2, IndexTrait::vector, IndexTrait::local, IndexTrait::md2macro>::divisionFactor{};

/*
 * !MD TO MACRO aka MAMICO INDEXING, EXCL GHOST LAYER
 */

// scalar, global, !md2macro, noGL
template <> BaseIndex<2> CellIndex<2, IndexTrait::noGhost>::lowerBoundary{};
template <> BaseIndex<2> CellIndex<2, IndexTrait::noGhost>::upperBoundary{};
template <> tarch::la::Vector<2, unsigned int> CellIndex<2, IndexTrait::noGhost>::numberCellsInDomain{};
template <> unsigned int CellIndex<2, IndexTrait::noGhost>::linearNumberCellsInDomain{};
template <> tarch::la::Vector<2, unsigned int> CellIndex<2, IndexTrait::noGhost>::divisionFactor{};

// vector, global, !md2macro, noGL
template <> BaseIndex<2> CellIndex<2, IndexTrait::vector, IndexTrait::noGhost>::lowerBoundary{};
template <> BaseIndex<2> CellIndex<2, IndexTrait::vector, IndexTrait::noGhost>::upperBoundary{};
template <> tarch::la::Vector<2, unsigned int> CellIndex<2, IndexTrait::vector, IndexTrait::noGhost>::numberCellsInDomain{};
template <> unsigned int CellIndex<2, IndexTrait::vector, IndexTrait::noGhost>::linearNumberCellsInDomain{};
template <> tarch::la::Vector<2, unsigned int> CellIndex<2, IndexTrait::vector, IndexTrait::noGhost>::divisionFactor{};

// scalar, local, !md2macro, noGL
template <> BaseIndex<2> CellIndex<2, IndexTrait::local, IndexTrait::noGhost>::lowerBoundary{};
template <> BaseIndex<2> CellIndex<2, IndexTrait::local, IndexTrait::noGhost>::upperBoundary{};
template <> tarch::la::Vector<2, unsigned int> CellIndex<2, IndexTrait::local, IndexTrait::noGhost>::numberCellsInDomain{};
template <> unsigned int CellIndex<2, IndexTrait::local, IndexTrait::noGhost>::linearNumberCellsInDomain{};
template <> tarch::la::Vector<2, unsigned int> CellIndex<2, IndexTrait::local, IndexTrait::noGhost>::divisionFactor{};

// vector, local, !md2macro, noGL
template <> BaseIndex<2> CellIndex<2, IndexTrait::vector, IndexTrait::local, IndexTrait::noGhost>::lowerBoundary{};
template <> BaseIndex<2> CellIndex<2, IndexTrait::vector, IndexTrait::local, IndexTrait::noGhost>::upperBoundary{};
template <> tarch::la::Vector<2, unsigned int> CellIndex<2, IndexTrait::vector, IndexTrait::local, IndexTrait::noGhost>::numberCellsInDomain{};
template <> unsigned int CellIndex<2, IndexTrait::vector, IndexTrait::local, IndexTrait::noGhost>::linearNumberCellsInDomain{};
template <> tarch::la::Vector<2, unsigned int> CellIndex<2, IndexTrait::vector, IndexTrait::local, IndexTrait::noGhost>::divisionFactor{};

/*
 * MD TO MACRO, EXCL GHOST LAYER
 */

// scalar, global, md2macro, noGL
template <> BaseIndex<2> CellIndex<2, IndexTrait::md2macro, IndexTrait::noGhost>::lowerBoundary{};
template <> BaseIndex<2> CellIndex<2, IndexTrait::md2macro, IndexTrait::noGhost>::upperBoundary{};
template <> tarch::la::Vector<2, unsigned int> CellIndex<2, IndexTrait::md2macro, IndexTrait::noGhost>::numberCellsInDomain{};
template <> unsigned int CellIndex<2, IndexTrait::md2macro, IndexTrait::noGhost>::linearNumberCellsInDomain{};
template <> tarch::la::Vector<2, unsigned int> CellIndex<2, IndexTrait::md2macro, IndexTrait::noGhost>::divisionFactor{};

// vector, global, md2macro, noGL
template <> BaseIndex<2> CellIndex<2, IndexTrait::vector, IndexTrait::md2macro, IndexTrait::noGhost>::lowerBoundary{};
template <> BaseIndex<2> CellIndex<2, IndexTrait::vector, IndexTrait::md2macro, IndexTrait::noGhost>::upperBoundary{};
template <> tarch::la::Vector<2, unsigned int> CellIndex<2, IndexTrait::vector, IndexTrait::md2macro, IndexTrait::noGhost>::numberCellsInDomain{};
template <> unsigned int CellIndex<2, IndexTrait::vector, IndexTrait::md2macro, IndexTrait::noGhost>::linearNumberCellsInDomain{};
template <> tarch::la::Vector<2, unsigned int> CellIndex<2, IndexTrait::vector, IndexTrait::md2macro, IndexTrait::noGhost>::divisionFactor{};

// scalar, local, md2macro, noGL
template <> BaseIndex<2> CellIndex<2, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::lowerBoundary{};
template <> BaseIndex<2> CellIndex<2, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::upperBoundary{};
template <> tarch::la::Vector<2, unsigned int> CellIndex<2, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::numberCellsInDomain{};
template <> unsigned int CellIndex<2, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::linearNumberCellsInDomain{};
template <> tarch::la::Vector<2, unsigned int> CellIndex<2, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::divisionFactor{};

// vector, local, md2macro, noGL
template <> BaseIndex<2> CellIndex<2, IndexTrait::vector, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::lowerBoundary{};
template <> BaseIndex<2> CellIndex<2, IndexTrait::vector, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::upperBoundary{};
template <>
tarch::la::Vector<2, unsigned int> CellIndex<2, IndexTrait::vector, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::numberCellsInDomain{};
template <> unsigned int CellIndex<2, IndexTrait::vector, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::linearNumberCellsInDomain{};
template <> tarch::la::Vector<2, unsigned int> CellIndex<2, IndexTrait::vector, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::divisionFactor{};
} // namespace indexing
} // namespace coupling
#endif

// Dim3   //////////////////////////
namespace coupling {
namespace indexing {

/*
 * Declare specialisations of CellIndex.
 * Define their static members.
 */

/*
 * NON-MD-TO-MACRO aka MAMICO INDEXING, INCL GHOST LAYER
 */

// scalar, global, !md2macro, !noGL
template <> BaseIndex<3> CellIndex<3>::lowerBoundary{};
template <> BaseIndex<3> CellIndex<3>::upperBoundary{};
template <> tarch::la::Vector<3, unsigned int> CellIndex<3>::numberCellsInDomain{};
template <> unsigned int CellIndex<3>::linearNumberCellsInDomain{};
template <> tarch::la::Vector<3, unsigned int> CellIndex<3>::divisionFactor{};
template <> const char CellIndex<3>::TNAME[] = "CellIndex<3>";

// BaseIndex
template <> BaseIndex<3> BaseIndex<3>::lowerBoundary{};
template <> BaseIndex<3> BaseIndex<3>::upperBoundary{};
template <> tarch::la::Vector<3, unsigned int> BaseIndex<3>::numberCellsInDomain{};
template <> unsigned int BaseIndex<3>::linearNumberCellsInDomain{};
template <> tarch::la::Vector<3, unsigned int> BaseIndex<3>::divisionFactor{};
template <> const char CellIndex<3, IndexTrait::vector>::TNAME[] = "CellIndex<3, vector>";

// scalar, local, !md2macro, !noGL
template <> BaseIndex<3> CellIndex<3, IndexTrait::local>::lowerBoundary{};
template <> BaseIndex<3> CellIndex<3, IndexTrait::local>::upperBoundary{};
template <> tarch::la::Vector<3, unsigned int> CellIndex<3, IndexTrait::local>::numberCellsInDomain{};
template <> unsigned int CellIndex<3, IndexTrait::local>::linearNumberCellsInDomain{};
template <> tarch::la::Vector<3, unsigned int> CellIndex<3, IndexTrait::local>::divisionFactor{};
template <> const char CellIndex<3, IndexTrait::local>::TNAME[] = "CellIndex<3, local>";

// vector, local, !md2macro, !noGL
template <> BaseIndex<3> CellIndex<3, IndexTrait::vector, IndexTrait::local>::lowerBoundary{};
template <> BaseIndex<3> CellIndex<3, IndexTrait::vector, IndexTrait::local>::upperBoundary{};
template <> tarch::la::Vector<3, unsigned int> CellIndex<3, IndexTrait::vector, IndexTrait::local>::numberCellsInDomain{};
template <> unsigned int CellIndex<3, IndexTrait::vector, IndexTrait::local>::linearNumberCellsInDomain{};
template <> tarch::la::Vector<3, unsigned int> CellIndex<3, IndexTrait::vector, IndexTrait::local>::divisionFactor{};
template <> const char CellIndex<3, IndexTrait::vector, IndexTrait::local>::TNAME[] = "CellIndex<3, vector, local>";

/*
 * MD TO MACRO, INCL GHOST LAYER
 */

// scalar, global, md2macro, !noGL
template <> BaseIndex<3> CellIndex<3, IndexTrait::md2macro>::lowerBoundary{};
template <> BaseIndex<3> CellIndex<3, IndexTrait::md2macro>::upperBoundary{};
template <> tarch::la::Vector<3, unsigned int> CellIndex<3, IndexTrait::md2macro>::numberCellsInDomain{};
template <> unsigned int CellIndex<3, IndexTrait::md2macro>::linearNumberCellsInDomain{};
template <> tarch::la::Vector<3, unsigned int> CellIndex<3, IndexTrait::md2macro>::divisionFactor{};
template <> const char CellIndex<3, IndexTrait::md2macro>::TNAME[] = "CellIndex<3, md2macro>";

// vector, global, md2macro, !noGL
template <> BaseIndex<3> CellIndex<3, IndexTrait::vector, IndexTrait::md2macro>::lowerBoundary{};
template <> BaseIndex<3> CellIndex<3, IndexTrait::vector, IndexTrait::md2macro>::upperBoundary{};
template <> tarch::la::Vector<3, unsigned int> CellIndex<3, IndexTrait::vector, IndexTrait::md2macro>::numberCellsInDomain{};
template <> unsigned int CellIndex<3, IndexTrait::vector, IndexTrait::md2macro>::linearNumberCellsInDomain{};
template <> tarch::la::Vector<3, unsigned int> CellIndex<3, IndexTrait::vector, IndexTrait::md2macro>::divisionFactor{};
template <> const char CellIndex<3, IndexTrait::vector, IndexTrait::md2macro>::TNAME[] = "CellIndex<3, vector, md2macro>";

// scalar, local, md2macro, !noGL
template <> BaseIndex<3> CellIndex<3, IndexTrait::local, IndexTrait::md2macro>::lowerBoundary{};
template <> BaseIndex<3> CellIndex<3, IndexTrait::local, IndexTrait::md2macro>::upperBoundary{};
template <> tarch::la::Vector<3, unsigned int> CellIndex<3, IndexTrait::local, IndexTrait::md2macro>::numberCellsInDomain{};
template <> unsigned int CellIndex<3, IndexTrait::local, IndexTrait::md2macro>::linearNumberCellsInDomain{};
template <> tarch::la::Vector<3, unsigned int> CellIndex<3, IndexTrait::local, IndexTrait::md2macro>::divisionFactor{};
template <> const char CellIndex<3, IndexTrait::local, IndexTrait::md2macro>::TNAME[] = "CellIndex<3, local, md2macro>";

// vector, local, md2macro, !noGL
template <> BaseIndex<3> CellIndex<3, IndexTrait::vector, IndexTrait::local, IndexTrait::md2macro>::lowerBoundary{};
template <> BaseIndex<3> CellIndex<3, IndexTrait::vector, IndexTrait::local, IndexTrait::md2macro>::upperBoundary{};
template <> tarch::la::Vector<3, unsigned int> CellIndex<3, IndexTrait::vector, IndexTrait::local, IndexTrait::md2macro>::numberCellsInDomain{};
template <> unsigned int CellIndex<3, IndexTrait::vector, IndexTrait::local, IndexTrait::md2macro>::linearNumberCellsInDomain{};
template <> tarch::la::Vector<3, unsigned int> CellIndex<3, IndexTrait::vector, IndexTrait::local, IndexTrait::md2macro>::divisionFactor{};
template <> const char CellIndex<3, IndexTrait::vector, IndexTrait::local, IndexTrait::md2macro>::TNAME[] = "CellIndex<3, vector, local, md2macro>";

/*
 * !MD TO MACRO aka MAMICO INDEXING, EXCL GHOST LAYER
 */

// scalar, global, !md2macro, noGL
template <> BaseIndex<3> CellIndex<3, IndexTrait::noGhost>::lowerBoundary{};
template <> BaseIndex<3> CellIndex<3, IndexTrait::noGhost>::upperBoundary{};
template <> tarch::la::Vector<3, unsigned int> CellIndex<3, IndexTrait::noGhost>::numberCellsInDomain{};
template <> unsigned int CellIndex<3, IndexTrait::noGhost>::linearNumberCellsInDomain{};
template <> tarch::la::Vector<3, unsigned int> CellIndex<3, IndexTrait::noGhost>::divisionFactor{};
template <> const char CellIndex<3, IndexTrait::noGhost>::TNAME[] = "CellIndex<3, noGhost>";

// vector, global, !md2macro, noGL
template <> BaseIndex<3> CellIndex<3, IndexTrait::vector, IndexTrait::noGhost>::lowerBoundary{};
template <> BaseIndex<3> CellIndex<3, IndexTrait::vector, IndexTrait::noGhost>::upperBoundary{};
template <> tarch::la::Vector<3, unsigned int> CellIndex<3, IndexTrait::vector, IndexTrait::noGhost>::numberCellsInDomain{};
template <> unsigned int CellIndex<3, IndexTrait::vector, IndexTrait::noGhost>::linearNumberCellsInDomain{};
template <> tarch::la::Vector<3, unsigned int> CellIndex<3, IndexTrait::vector, IndexTrait::noGhost>::divisionFactor{};
template <> const char CellIndex<3, IndexTrait::vector, IndexTrait::noGhost>::TNAME[] = "CellIndex<3, vector, noGhost>";

// scalar, local, !md2macro, noGL
template <> BaseIndex<3> CellIndex<3, IndexTrait::local, IndexTrait::noGhost>::lowerBoundary{};
template <> BaseIndex<3> CellIndex<3, IndexTrait::local, IndexTrait::noGhost>::upperBoundary{};
template <> tarch::la::Vector<3, unsigned int> CellIndex<3, IndexTrait::local, IndexTrait::noGhost>::numberCellsInDomain{};
template <> unsigned int CellIndex<3, IndexTrait::local, IndexTrait::noGhost>::linearNumberCellsInDomain{};
template <> tarch::la::Vector<3, unsigned int> CellIndex<3, IndexTrait::local, IndexTrait::noGhost>::divisionFactor{};
template <> const char CellIndex<3, IndexTrait::local, IndexTrait::noGhost>::TNAME[] = "CellIndex<3, local, noGhost>";

// vector, local, !md2macro, noGL
template <> BaseIndex<3> CellIndex<3, IndexTrait::vector, IndexTrait::local, IndexTrait::noGhost>::lowerBoundary{};
template <> BaseIndex<3> CellIndex<3, IndexTrait::vector, IndexTrait::local, IndexTrait::noGhost>::upperBoundary{};
template <> tarch::la::Vector<3, unsigned int> CellIndex<3, IndexTrait::vector, IndexTrait::local, IndexTrait::noGhost>::numberCellsInDomain{};
template <> unsigned int CellIndex<3, IndexTrait::vector, IndexTrait::local, IndexTrait::noGhost>::linearNumberCellsInDomain{};
template <> tarch::la::Vector<3, unsigned int> CellIndex<3, IndexTrait::vector, IndexTrait::local, IndexTrait::noGhost>::divisionFactor{};
template <> const char CellIndex<3, IndexTrait::vector, IndexTrait::local, IndexTrait::noGhost>::TNAME[] = "CellIndex<3, vector, local, noGhost>";

/*
 * MD TO MACRO, EXCL GHOST LAYER
 */

// scalar, global, md2macro, noGL
template <> BaseIndex<3> CellIndex<3, IndexTrait::md2macro, IndexTrait::noGhost>::lowerBoundary{};
template <> BaseIndex<3> CellIndex<3, IndexTrait::md2macro, IndexTrait::noGhost>::upperBoundary{};
template <> tarch::la::Vector<3, unsigned int> CellIndex<3, IndexTrait::md2macro, IndexTrait::noGhost>::numberCellsInDomain{};
template <> unsigned int CellIndex<3, IndexTrait::md2macro, IndexTrait::noGhost>::linearNumberCellsInDomain{};
template <> tarch::la::Vector<3, unsigned int> CellIndex<3, IndexTrait::md2macro, IndexTrait::noGhost>::divisionFactor{};
template <> const char CellIndex<3, IndexTrait::md2macro, IndexTrait::noGhost>::TNAME[] = "CellIndex<3, md2macro, noGhost>";

// vector, global, md2macro, noGL
template <> BaseIndex<3> CellIndex<3, IndexTrait::vector, IndexTrait::md2macro, IndexTrait::noGhost>::lowerBoundary{};
template <> BaseIndex<3> CellIndex<3, IndexTrait::vector, IndexTrait::md2macro, IndexTrait::noGhost>::upperBoundary{};
template <> tarch::la::Vector<3, unsigned int> CellIndex<3, IndexTrait::vector, IndexTrait::md2macro, IndexTrait::noGhost>::numberCellsInDomain{};
template <> unsigned int CellIndex<3, IndexTrait::vector, IndexTrait::md2macro, IndexTrait::noGhost>::linearNumberCellsInDomain{};
template <> tarch::la::Vector<3, unsigned int> CellIndex<3, IndexTrait::vector, IndexTrait::md2macro, IndexTrait::noGhost>::divisionFactor{};
template <> const char CellIndex<3, IndexTrait::vector, IndexTrait::md2macro, IndexTrait::noGhost>::TNAME[] = "CellIndex<3, vector, md2macro, noGhost>";

// scalar, local, md2macro, noGL
template <> BaseIndex<3> CellIndex<3, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::lowerBoundary{};
template <> BaseIndex<3> CellIndex<3, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::upperBoundary{};
template <> tarch::la::Vector<3, unsigned int> CellIndex<3, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::numberCellsInDomain{};
template <> unsigned int CellIndex<3, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::linearNumberCellsInDomain{};
template <> tarch::la::Vector<3, unsigned int> CellIndex<3, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::divisionFactor{};
template <> const char CellIndex<3, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::TNAME[] = "CellIndex<3, local, md2macro, noGhost>";

// vector, local, md2macro, noGL
template <> BaseIndex<3> CellIndex<3, IndexTrait::vector, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::lowerBoundary{};
template <> BaseIndex<3> CellIndex<3, IndexTrait::vector, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::upperBoundary{};
template <>
tarch::la::Vector<3, unsigned int> CellIndex<3, IndexTrait::vector, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::numberCellsInDomain{};
template <> unsigned int CellIndex<3, IndexTrait::vector, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::linearNumberCellsInDomain{};
template <> tarch::la::Vector<3, unsigned int> CellIndex<3, IndexTrait::vector, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::divisionFactor{};
template <>
const char CellIndex<3, IndexTrait::vector, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::TNAME[] =
    "CellIndex<3, vector, local, md2macro, noGhost>";
} // namespace indexing
} // namespace coupling

// raw date based variant of init
template <unsigned int dim>
void coupling::indexing::IndexingService<dim>::init(tarch::la::Vector<dim, unsigned int> globalNumberMacroscopicCells,
                                                    tarch::la::Vector<dim, unsigned int> numberProcesses, coupling::paralleltopology::ParallelTopologyType type,
                                                    unsigned int outerRegion, const unsigned int rank
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
                                                    ,
                                                    MPI_Comm comm
#endif
) {
  coupling::interface::MacroscopicSolverInterface<dim>* msi = new coupling::solvers::CouetteSolverInterface<dim>(globalNumberMacroscopicCells, outerRegion);

  init(globalNumberMacroscopicCells, numberProcesses, type, msi, rank
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
       ,
       comm
#endif
  );
}

// delegated init, this does the main work
template <unsigned int dim>
void coupling::indexing::IndexingService<dim>::init(tarch::la::Vector<dim, unsigned int> globalNumberMacroscopicCells,
                                                    tarch::la::Vector<dim, unsigned int> numberProcesses,
                                                    coupling::paralleltopology::ParallelTopologyType parallelTopologyType,
                                                    coupling::interface::MacroscopicSolverInterface<dim>* msi, const unsigned int rank
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
                                                    ,
                                                    MPI_Comm comm
#endif
) {

  _rank = rank;
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
  _comm = comm;
#endif

#if (COUPLING_MD_ERROR == COUPLING_MD_YES)
  if (_isInitialized) {
    std::cout << "IndexingService: WARNING: Initializing twice! " << std::endl;
  }
#endif

  // TODO: make this globalNumberMacroscopicCells and remove all usages of the
  // old meaning (seen above)
  const auto globalNumberMacroscopicCellsInclGL{globalNumberMacroscopicCells + tarch::la::Vector<dim, unsigned int>{2}};

  // init boundaries of all global, non-m2m, GL including indexing types
  CellIndex<dim>::lowerBoundary = {0};
  CellIndex<dim>::upperBoundary = tarch::la::Vector<dim, int>{globalNumberMacroscopicCellsInclGL - tarch::la::Vector<dim, unsigned int>{1}};
  CellIndex<dim>::setDomainParameters();

  CellIndex<dim, IndexTrait::vector>::lowerBoundary = CellIndex<dim>::lowerBoundary;
  CellIndex<dim, IndexTrait::vector>::upperBoundary = CellIndex<dim>::upperBoundary;
  CellIndex<dim, IndexTrait::vector>::setDomainParameters();

  // init boundaries of all global, non-m2m, GL excluding indexing types
  CellIndex<dim, IndexTrait::noGhost>::lowerBoundary = {1};
  CellIndex<dim, IndexTrait::noGhost>::upperBoundary =
      tarch::la::Vector<dim, int>{globalNumberMacroscopicCellsInclGL - tarch::la::Vector<dim, unsigned int>{2}};
  CellIndex<dim, IndexTrait::noGhost>::setDomainParameters();

  CellIndex<dim, IndexTrait::vector, IndexTrait::noGhost>::lowerBoundary = CellIndex<dim, IndexTrait::noGhost>::lowerBoundary;
  CellIndex<dim, IndexTrait::vector, IndexTrait::noGhost>::upperBoundary = CellIndex<dim, IndexTrait::noGhost>::upperBoundary;
  CellIndex<dim, IndexTrait::vector, IndexTrait::noGhost>::setDomainParameters();

// CellIndex<dim> and CellIndex<dim, IndexTrait::vector> have been set up, so from here on it is ok to do some basic conversions
#if (COUPLING_MD_ERROR == COUPLING_MD_YES)
  _isInitialized = true;
#endif

  // init boundaries of all global, m2m, GL excluding indexing types
  {
    CellIndex<dim> lowerBoundary{BaseIndex<dim>::lowerBoundary};
    while (not msi->receiveMacroscopicQuantityFromMDSolver(tarch::la::Vector<dim, unsigned int>{CellIndex<dim, IndexTrait::vector>{lowerBoundary}.get()})) {
      // sanity check: empty m2m domain
      if (lowerBoundary == BaseIndex<dim>::upperBoundary) {
        std::cout << "IndexingService: WARNING: Empty MD-To-Macro domain!" << std::endl;
        break;
      }

      // increment by one if above is too low to be in md-to-macro domain
      ++lowerBoundary;
    }
    CellIndex<dim> upperBoundary{BaseIndex<dim>::upperBoundary};
    while (not msi->receiveMacroscopicQuantityFromMDSolver(tarch::la::Vector<dim, unsigned int>{CellIndex<dim, IndexTrait::vector>{upperBoundary}.get()})) {
      // sanity check: empty m2m domain
      if (upperBoundary < lowerBoundary) {
        std::cout << "IndexingService: WARNING: Empty MD-To-Macro domain!" << std::endl;
        break;
      }

      // decrement by one if above is too high to be in md-to-macro domain
      --upperBoundary;
    }

    CellIndex<dim, IndexTrait::md2macro, IndexTrait::noGhost>::lowerBoundary = lowerBoundary;
    CellIndex<dim, IndexTrait::md2macro, IndexTrait::noGhost>::upperBoundary = upperBoundary;
    CellIndex<dim, IndexTrait::md2macro, IndexTrait::noGhost>::setDomainParameters();
  }

  CellIndex<dim, IndexTrait::vector, IndexTrait::md2macro, IndexTrait::noGhost>::lowerBoundary =
      CellIndex<dim, IndexTrait::md2macro, IndexTrait::noGhost>::lowerBoundary;
  CellIndex<dim, IndexTrait::vector, IndexTrait::md2macro, IndexTrait::noGhost>::upperBoundary =
      CellIndex<dim, IndexTrait::md2macro, IndexTrait::noGhost>::upperBoundary;
  CellIndex<dim, IndexTrait::vector, IndexTrait::md2macro, IndexTrait::noGhost>::setDomainParameters();

  // init boundaries of all global, m2m, GL including indexing types
  CellIndex<dim, IndexTrait::md2macro>::lowerBoundary =
      CellIndex<dim, IndexTrait::md2macro, IndexTrait::noGhost>::lowerBoundary.get() - tarch::la::Vector<dim, int>{1};
  CellIndex<dim, IndexTrait::md2macro>::upperBoundary =
      CellIndex<dim, IndexTrait::md2macro, IndexTrait::noGhost>::upperBoundary.get() + tarch::la::Vector<dim, int>{1};
  CellIndex<dim, IndexTrait::md2macro>::setDomainParameters();

  CellIndex<dim, IndexTrait::vector, IndexTrait::md2macro>::lowerBoundary = CellIndex<dim, IndexTrait::md2macro>::lowerBoundary;
  CellIndex<dim, IndexTrait::vector, IndexTrait::md2macro>::upperBoundary = CellIndex<dim, IndexTrait::md2macro>::upperBoundary;
  CellIndex<dim, IndexTrait::vector, IndexTrait::md2macro>::setDomainParameters();

// handle all local indexing types
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES) // parallel scenario

  _numberProcesses = numberProcesses;

  unsigned int scalarNumberProcesses = _numberProcesses[0];
  for (unsigned int d = 1; d < dim; d++)
    scalarNumberProcesses *= _numberProcesses[d];

  const unsigned int parallelTopologyOffset = (_rank / scalarNumberProcesses) * scalarNumberProcesses; // copied from IndexConversion
  _parallelTopology =
      coupling::paralleltopology::ParallelTopologyFactory::getParallelTopology<dim>(parallelTopologyType, _numberProcesses, parallelTopologyOffset);

  std::vector<unsigned int> ranks; // used to store ranks in which certain indices occur

  // init boundaries of all local, non-m2m, GL including indexing types
  {
    BaseIndex<dim> lowerBoundary{CellIndex<dim /*global*/>::lowerBoundary};
    BaseIndex<dim> upperBoundary{CellIndex<dim /*global*/>::upperBoundary};

    // TODO: determine these two analyticaly (i.e. calculate domain bounds): you
    // dont need to iterate over all global indices.
    while (true) {
      ranks = getRanksForGlobalIndex(lowerBoundary);
      if (std::find(ranks.begin(), ranks.end(), _rank) != ranks.end()) /*if _rank is found in ranks in which the tested index
                                                                          occurs...*/
        break;

      //...increment by one if above is too high to be in local domain
      ++lowerBoundary;
    }
    while (true) {
      ranks = getRanksForGlobalIndex(upperBoundary);
      if (std::find(ranks.begin(), ranks.end(), _rank) != ranks.end()) /*if _rank is found in ranks in which the tested index
                                                                          occurs...*/
        break;

      //...decrement by one if above is too high to be in local domain
      --upperBoundary;
    }

    CellIndex<dim, IndexTrait::local>::lowerBoundary = lowerBoundary;
    CellIndex<dim, IndexTrait::local>::upperBoundary = upperBoundary;
    CellIndex<dim, IndexTrait::local>::setDomainParameters();
  }

  CellIndex<dim, IndexTrait::vector, IndexTrait::local>::lowerBoundary = CellIndex<dim, IndexTrait::local>::lowerBoundary;
  CellIndex<dim, IndexTrait::vector, IndexTrait::local>::upperBoundary = CellIndex<dim, IndexTrait::local>::upperBoundary;
  CellIndex<dim, IndexTrait::vector, IndexTrait::local>::setDomainParameters();

  // init boundaries of all local, non-m2m, GL excluding indexing types
  CellIndex<dim, IndexTrait::local, IndexTrait::noGhost>::lowerBoundary =
      CellIndex<dim, IndexTrait::local>::lowerBoundary.get() + tarch::la::Vector<dim, int>{1};
  CellIndex<dim, IndexTrait::local, IndexTrait::noGhost>::upperBoundary =
      CellIndex<dim, IndexTrait::local>::upperBoundary.get() - tarch::la::Vector<dim, int>{1};
  CellIndex<dim, IndexTrait::local, IndexTrait::noGhost>::setDomainParameters();

  CellIndex<dim, IndexTrait::vector, IndexTrait::local, IndexTrait::noGhost>::lowerBoundary =
      CellIndex<dim, IndexTrait::local, IndexTrait::noGhost>::lowerBoundary;
  CellIndex<dim, IndexTrait::vector, IndexTrait::local, IndexTrait::noGhost>::upperBoundary =
      CellIndex<dim, IndexTrait::local, IndexTrait::noGhost>::upperBoundary;
  CellIndex<dim, IndexTrait::vector, IndexTrait::local, IndexTrait::noGhost>::setDomainParameters();

  // init boundaries of all local, m2m, GL excluding indexing types
  {
    using LocalNoGlIndex = CellIndex<dim, IndexTrait::local, IndexTrait::noGhost>;
    using GlobalMD2MacroNoGlIndex = CellIndex<dim, IndexTrait::md2macro, IndexTrait::noGhost>;

    tarch::la::Vector<dim, int> lowerBoundary;
    tarch::la::Vector<dim, int> upperBoundary;
    for (unsigned int d = 0; d < dim; d++) {
      lowerBoundary[d] = std::max(LocalNoGlIndex::lowerBoundary.get()[d], GlobalMD2MacroNoGlIndex::lowerBoundary.get()[d]);
      upperBoundary[d] = std::min(LocalNoGlIndex::upperBoundary.get()[d], GlobalMD2MacroNoGlIndex::upperBoundary.get()[d]);
    }

    CellIndex<dim, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::lowerBoundary = lowerBoundary;
    CellIndex<dim, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::upperBoundary = upperBoundary;
    CellIndex<dim, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::setDomainParameters();
  }

  CellIndex<dim, IndexTrait::vector, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::lowerBoundary =
      CellIndex<dim, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::lowerBoundary;
  CellIndex<dim, IndexTrait::vector, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::upperBoundary =
      CellIndex<dim, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::upperBoundary;
  CellIndex<dim, IndexTrait::vector, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::setDomainParameters();

  // init boundaries of all local, m2m, GL including indexing types
  CellIndex<dim, IndexTrait::local, IndexTrait::md2macro>::lowerBoundary =
      CellIndex<dim, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::lowerBoundary.get() - tarch::la::Vector<dim, int>{1};
  CellIndex<dim, IndexTrait::local, IndexTrait::md2macro>::upperBoundary =
      CellIndex<dim, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::upperBoundary.get() + tarch::la::Vector<dim, int>{1};
  CellIndex<dim, IndexTrait::local, IndexTrait::md2macro>::setDomainParameters();

  CellIndex<dim, IndexTrait::vector, IndexTrait::local, IndexTrait::md2macro>::lowerBoundary =
      CellIndex<dim, IndexTrait::local, IndexTrait::md2macro>::lowerBoundary;
  CellIndex<dim, IndexTrait::vector, IndexTrait::local, IndexTrait::md2macro>::upperBoundary =
      CellIndex<dim, IndexTrait::local, IndexTrait::md2macro>::upperBoundary;
  CellIndex<dim, IndexTrait::vector, IndexTrait::local, IndexTrait::md2macro>::setDomainParameters();

#else // sequential scenario
      // Copy all local indexing from global
  CellIndex<dim, IndexTrait::local>::lowerBoundary = CellIndex<dim>::lowerBoundary;
  CellIndex<dim, IndexTrait::local>::upperBoundary = CellIndex<dim>::upperBoundary;
  CellIndex<dim, IndexTrait::local>::setDomainParameters();

  CellIndex<dim, IndexTrait::vector, IndexTrait::local>::lowerBoundary = CellIndex<dim, IndexTrait::vector>::lowerBoundary;
  CellIndex<dim, IndexTrait::vector, IndexTrait::local>::upperBoundary = CellIndex<dim, IndexTrait::vector>::upperBoundary;
  CellIndex<dim, IndexTrait::vector, IndexTrait::local>::setDomainParameters();

  CellIndex<dim, IndexTrait::local, IndexTrait::noGhost>::lowerBoundary = CellIndex<dim, IndexTrait::noGhost>::lowerBoundary;
  CellIndex<dim, IndexTrait::local, IndexTrait::noGhost>::upperBoundary = CellIndex<dim, IndexTrait::noGhost>::upperBoundary;
  CellIndex<dim, IndexTrait::local, IndexTrait::noGhost>::setDomainParameters();

  CellIndex<dim, IndexTrait::vector, IndexTrait::local, IndexTrait::noGhost>::lowerBoundary =
      CellIndex<dim, IndexTrait::vector, IndexTrait::noGhost>::lowerBoundary;
  CellIndex<dim, IndexTrait::vector, IndexTrait::local, IndexTrait::noGhost>::upperBoundary =
      CellIndex<dim, IndexTrait::vector, IndexTrait::noGhost>::upperBoundary;
  CellIndex<dim, IndexTrait::vector, IndexTrait::local, IndexTrait::noGhost>::setDomainParameters();

  CellIndex<dim, IndexTrait::local, IndexTrait::md2macro>::lowerBoundary = CellIndex<dim, IndexTrait::md2macro>::lowerBoundary;
  CellIndex<dim, IndexTrait::local, IndexTrait::md2macro>::upperBoundary = CellIndex<dim, IndexTrait::md2macro>::upperBoundary;
  CellIndex<dim, IndexTrait::local, IndexTrait::md2macro>::setDomainParameters();

  CellIndex<dim, IndexTrait::vector, IndexTrait::local, IndexTrait::md2macro>::lowerBoundary =
      CellIndex<dim, IndexTrait::vector, IndexTrait::md2macro>::lowerBoundary;
  CellIndex<dim, IndexTrait::vector, IndexTrait::local, IndexTrait::md2macro>::upperBoundary =
      CellIndex<dim, IndexTrait::vector, IndexTrait::md2macro>::upperBoundary;
  CellIndex<dim, IndexTrait::vector, IndexTrait::local, IndexTrait::md2macro>::setDomainParameters();

  CellIndex<dim, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::lowerBoundary =
      CellIndex<dim, IndexTrait::md2macro, IndexTrait::noGhost>::lowerBoundary;
  CellIndex<dim, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::upperBoundary =
      CellIndex<dim, IndexTrait::md2macro, IndexTrait::noGhost>::upperBoundary;
  CellIndex<dim, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::setDomainParameters();

  CellIndex<dim, IndexTrait::vector, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::lowerBoundary =
      CellIndex<dim, IndexTrait::vector, IndexTrait::md2macro, IndexTrait::noGhost>::lowerBoundary;
  CellIndex<dim, IndexTrait::vector, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::upperBoundary =
      CellIndex<dim, IndexTrait::vector, IndexTrait::md2macro, IndexTrait::noGhost>::upperBoundary;
  CellIndex<dim, IndexTrait::vector, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::setDomainParameters();
#endif
}

#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES) // unused in sequential scenario
/*
 * This was in large parts stolen from IndexConversion.
 */
template <unsigned int dim>
std::vector<unsigned int> coupling::indexing::IndexingService<dim>::getRanksForGlobalIndex(const BaseIndex<dim>& globalCellIndex) const {

#if (COUPLING_MD_ERROR == COUPLING_MD_YES)
  if (!_isInitialized) {
    throw std::runtime_error(std::string("coupling::indexing::convertToVector: IndexingService not initialized! "));
  }
#endif

  std::vector<unsigned int> ranks;
  // using the old meaning of 'globalNumberMacroscopicCells' from
  // IndexConversion
  const auto globalNumberMacroscopicCells = BaseIndex<dim>::numberCellsInDomain - tarch::la::Vector<dim, unsigned int>{2};

  // start and end coordinates of neighboured cells.
  tarch::la::Vector<3, unsigned int> start(0);
  tarch::la::Vector<3, unsigned int> end(0);
  tarch::la::Vector<3, unsigned int> loopIndex(0);

  // determine up to 3^dim neighboured cells in the surrounding of
  // globalCellIndex; reduce this number if globalCellIndex lies on the global
  // boundary
  for (unsigned int d = 0; d < dim; d++) {
    if ((unsigned int)globalCellIndex.get()[d] > 0) {
      start[d] = (unsigned int)globalCellIndex.get()[d] - 1;
    }
    end[d] = globalNumberMacroscopicCells[d] + 1;
    if ((unsigned int)globalCellIndex.get()[d] < end[d]) {
      end[d] = (unsigned int)globalCellIndex.get()[d] + 1;
    }
  }

  // loop over neighbouring regions
  for (loopIndex[2] = start[2]; loopIndex[2] <= end[2]; loopIndex[2]++) {
    for (loopIndex[1] = start[1]; loopIndex[1] <= end[1]; loopIndex[1]++) {
      for (loopIndex[0] = start[0]; loopIndex[0] <= end[0]; loopIndex[0]++) {

        // determine the global cell index of this particular grid cell
        tarch::la::Vector<dim, unsigned int> thisGlobalCellIndex(0);
        for (unsigned int d = 0; d < dim; d++) {
          thisGlobalCellIndex[d] = loopIndex[d];
        }

        // determine the unique rank for this cell
        const unsigned int rank = getUniqueRankForMacroscopicCell(thisGlobalCellIndex, globalNumberMacroscopicCells);

        // add this rank to the vector with all ranks if we did not add this one
        // before
        bool isContained = false;
        const unsigned int thisSize = (unsigned int)ranks.size();
        for (unsigned int i = 0; i < thisSize; i++) {
          if (ranks[i] == rank) {
            isContained = true;
            break;
          }
        }
        if (!isContained) {
          ranks.push_back(rank);
        }
      }
    }
  }

  return ranks;
}
#endif

#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES) // unused in sequential scenario
/*
 * This was in large parts stolen from IndexConversion.
 * Note that this uses the globalNumberMacroscopicCells definition excl. the
 * ghost layer.
 */
template <unsigned int dim>
unsigned int
coupling::indexing::IndexingService<dim>::getUniqueRankForMacroscopicCell(tarch::la::Vector<dim, unsigned int> globalCellIndex,
                                                                          const tarch::la::Vector<dim, unsigned int>& globalNumberMacroscopicCells) const {
  // vector containing avg number of macro cells, not counting global GL.
  tarch::la::Vector<dim, unsigned int> averageLocalNumberMacroscopicCells{0};
  for (unsigned int d = 0; d < dim; d++) {
    if (globalCellIndex[d] >= globalNumberMacroscopicCells[d] + 2) { // greater or equal to the total global number incl GL (+2)
      using namespace std::string_literals;
      throw std::runtime_error("IndexingService: getUniqueRankForMacroscopicCell(): Global cell index greater than global size in dim "s + std::to_string(d));
    }
    if (globalNumberMacroscopicCells[d] % _numberProcesses[d] != 0) {
      std::stringstream ss;
      ss << "IndexingService: getUniqueRankForMacroscopicCell(): ERROR: Number "
            "of macroscopic cells must be divisible by number of processes! ";
      ss << "globalNumberMacroscopicCells = " << globalNumberMacroscopicCells;
      ss << ", numberProcesses = " << _numberProcesses;
      throw std::runtime_error(ss.str());
    }
    averageLocalNumberMacroscopicCells[d] = globalNumberMacroscopicCells[d] / _numberProcesses[d];
  }

  tarch::la::Vector<dim, unsigned int> processCoords(0);
  for (unsigned int d = 0; d < dim; d++) {
    // special case: cell in first section
    if (globalCellIndex[d] < averageLocalNumberMacroscopicCells[d] + 1) {
      processCoords[d] = 0;
      // special case: cell in last section
    } else if (globalCellIndex[d] > averageLocalNumberMacroscopicCells[d] * (_numberProcesses[d] - 1)) {
      processCoords[d] = _numberProcesses[d] - 1;
      // all other cases
    } else {
      // remove ghost layer contribution from vector index (...-1)
      processCoords[d] = (globalCellIndex[d] - 1) / averageLocalNumberMacroscopicCells[d];
    }
  }

  return _parallelTopology->getRank(processCoords);
}
#endif

// declare specialisation of IndexingService
#ifdef INDEXING_ENABLE_DIM2
template class coupling::indexing::IndexingService<2>;
#endif
template class coupling::indexing::IndexingService<3>;
