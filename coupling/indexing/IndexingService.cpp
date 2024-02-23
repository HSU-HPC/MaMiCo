// Include header
#include "IndexingService.h"

#include <algorithm>
#include <iterator>
#include <numeric>

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
template <> BaseIndex<3> I00::lowerBoundary{};
template <> BaseIndex<3> I00::upperBoundary{};
template <> tarch::la::Vector<3, unsigned int> I00::numberCellsInDomain{};
template <> unsigned int I00::linearNumberCellsInDomain{};
template <> tarch::la::Vector<3, unsigned int> I00::divisionFactor{};
template <> const char I00::TNAME[] = "CellIndex<3>";

// BaseIndex
template <> BaseIndex<3> BaseIndex<3>::lowerBoundary{};
template <> BaseIndex<3> BaseIndex<3>::upperBoundary{};
template <> tarch::la::Vector<3, unsigned int> BaseIndex<3>::numberCellsInDomain{};
template <> unsigned int BaseIndex<3>::linearNumberCellsInDomain{};
template <> tarch::la::Vector<3, unsigned int> BaseIndex<3>::divisionFactor{};
template <> const char CellIndex<3, IndexTrait::vector>::TNAME[] = "CellIndex<3, vector>";

// scalar, local, !md2macro, !noGL
template <> BaseIndex<3> I02::lowerBoundary{};
template <> BaseIndex<3> I02::upperBoundary{};
template <> tarch::la::Vector<3, unsigned int> I02::numberCellsInDomain{};
template <> unsigned int I02::linearNumberCellsInDomain{};
template <> tarch::la::Vector<3, unsigned int> I02::divisionFactor{};
template <> const char I02::TNAME[] = "CellIndex<3, local>";

// vector, local, !md2macro, !noGL
template <> BaseIndex<3> I03::lowerBoundary{};
template <> BaseIndex<3> I03::upperBoundary{};
template <> tarch::la::Vector<3, unsigned int> I03::numberCellsInDomain{};
template <> unsigned int I03::linearNumberCellsInDomain{};
template <> tarch::la::Vector<3, unsigned int> I03::divisionFactor{};
template <> const char I03::TNAME[] = "CellIndex<3, vector, local>";

/*
 * MD TO MACRO, INCL GHOST LAYER
 */

// scalar, global, md2macro, !noGL
template <> BaseIndex<3> I04::lowerBoundary{};
template <> BaseIndex<3> I04::upperBoundary{};
template <> tarch::la::Vector<3, unsigned int> I04::numberCellsInDomain{};
template <> unsigned int I04::linearNumberCellsInDomain{};
template <> tarch::la::Vector<3, unsigned int> I04::divisionFactor{};
template <> const char I04::TNAME[] = "CellIndex<3, md2macro>";

// vector, global, md2macro, !noGL
template <> BaseIndex<3> I05::lowerBoundary{};
template <> BaseIndex<3> I05::upperBoundary{};
template <> tarch::la::Vector<3, unsigned int> I05::numberCellsInDomain{};
template <> unsigned int I05::linearNumberCellsInDomain{};
template <> tarch::la::Vector<3, unsigned int> I05::divisionFactor{};
template <> const char I05::TNAME[] = "CellIndex<3, vector, md2macro>";

// scalar, local, md2macro, !noGL
template <> BaseIndex<3> I06::lowerBoundary{};
template <> BaseIndex<3> I06::upperBoundary{};
template <> tarch::la::Vector<3, unsigned int> I06::numberCellsInDomain{};
template <> unsigned int I06 ::linearNumberCellsInDomain{};
template <> tarch::la::Vector<3, unsigned int> I06::divisionFactor{};
template <> const char I06::TNAME[] = "CellIndex<3, local, md2macro>";

// vector, local, md2macro, !noGL
template <> BaseIndex<3> I07::lowerBoundary{};
template <> BaseIndex<3> I07::upperBoundary{};
template <> tarch::la::Vector<3, unsigned int> I07::numberCellsInDomain{};
template <> unsigned int I07::linearNumberCellsInDomain{};
template <> tarch::la::Vector<3, unsigned int> I07::divisionFactor{};
template <> const char I07::TNAME[] = "CellIndex<3, vector, local, md2macro>";

/*
 * !MD TO MACRO aka MAMICO INDEXING, EXCL GHOST LAYER
 */

// scalar, global, !md2macro, noGL
template <> BaseIndex<3> I08::lowerBoundary{};
template <> BaseIndex<3> I08::upperBoundary{};
template <> tarch::la::Vector<3, unsigned int> I08::numberCellsInDomain{};
template <> unsigned int I08::linearNumberCellsInDomain{};
template <> tarch::la::Vector<3, unsigned int> I08::divisionFactor{};
template <> const char I08::TNAME[] = "CellIndex<3, noGhost>";

// vector, global, !md2macro, noGL
template <> BaseIndex<3> I09::lowerBoundary{};
template <> BaseIndex<3> I09::upperBoundary{};
template <> tarch::la::Vector<3, unsigned int> I09::numberCellsInDomain{};
template <> unsigned int I09::linearNumberCellsInDomain{};
template <> tarch::la::Vector<3, unsigned int> I09::divisionFactor{};
template <> const char I09::TNAME[] = "CellIndex<3, vector, noGhost>";

// scalar, local, !md2macro, noGL
template <> BaseIndex<3> I10::lowerBoundary{};
template <> BaseIndex<3> I10::upperBoundary{};
template <> tarch::la::Vector<3, unsigned int> I10::numberCellsInDomain{};
template <> unsigned int I10::linearNumberCellsInDomain{};
template <> tarch::la::Vector<3, unsigned int> I10::divisionFactor{};
template <> const char I10::TNAME[] = "CellIndex<3, local, noGhost>";

// vector, local, !md2macro, noGL
template <> BaseIndex<3> I11::lowerBoundary{};
template <> BaseIndex<3> I11::upperBoundary{};
template <> tarch::la::Vector<3, unsigned int> I11::numberCellsInDomain{};
template <> unsigned int I11::linearNumberCellsInDomain{};
template <> tarch::la::Vector<3, unsigned int> I11::divisionFactor{};
template <> const char I11::TNAME[] = "CellIndex<3, vector, local, noGhost>";

/*
 * MD TO MACRO, EXCL GHOST LAYER
 */

// scalar, global, md2macro, noGL
template <> BaseIndex<3> I12::lowerBoundary{};
template <> BaseIndex<3> I12::upperBoundary{};
template <> tarch::la::Vector<3, unsigned int> I12::numberCellsInDomain{};
template <> unsigned int I12::linearNumberCellsInDomain{};
template <> tarch::la::Vector<3, unsigned int> I12::divisionFactor{};
template <> const char I12::TNAME[] = "CellIndex<3, md2macro, noGhost>";

// vector, global, md2macro, noGL
template <> BaseIndex<3> I13::lowerBoundary{};
template <> BaseIndex<3> I13::upperBoundary{};
template <> tarch::la::Vector<3, unsigned int> I13::numberCellsInDomain{};
template <> unsigned int I13::linearNumberCellsInDomain{};
template <> tarch::la::Vector<3, unsigned int> I13::divisionFactor{};
template <> const char I13::TNAME[] = "CellIndex<3, vector, md2macro, noGhost>";

// scalar, local, md2macro, noGL
template <> BaseIndex<3> I14::lowerBoundary{};
template <> BaseIndex<3> I14::upperBoundary{};
template <> tarch::la::Vector<3, unsigned int> I14::numberCellsInDomain{};
template <> unsigned int I14::linearNumberCellsInDomain{};
template <> tarch::la::Vector<3, unsigned int> I14::divisionFactor{};
template <> const char I14::TNAME[] = "CellIndex<3, local, md2macro, noGhost>";

// vector, local, md2macro, noGL
template <> BaseIndex<3> I15::lowerBoundary{};
template <> BaseIndex<3> I15::upperBoundary{};
template <> tarch::la::Vector<3, unsigned int> I15::numberCellsInDomain{};
template <> unsigned int I15::linearNumberCellsInDomain{};
template <> tarch::la::Vector<3, unsigned int> I15::divisionFactor{};
template <> const char I15::TNAME[] = "CellIndex<3, vector, local, md2macro, noGhost>";
} // namespace indexing
} // namespace coupling

// delegated init, this does the main work
template <unsigned int dim>
void coupling::indexing::IndexingService<dim>::initWithCells(tarch::la::Vector<dim, std::vector<unsigned int>>& subdomainWeights,
                                                             tarch::la::Vector<dim, unsigned int> globalNumberCouplingCells,
                                                             tarch::la::Vector<dim, unsigned int> numberProcesses,
                                                             const tarch::la::Vector<3, double>& couplingCellSize,
                                                             coupling::paralleltopology::ParallelTopologyType parallelTopologyType, unsigned int outerRegion,
                                                             const unsigned int rank
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

  _couplingCellSize = couplingCellSize;

  // TODO: make this globalNumberCouplingCells and remove all usages of the
  // old meaning (seen above)
  const auto globalNumberCouplingCellsInclGL{globalNumberCouplingCells + tarch::la::Vector<dim, unsigned int>{2}};

  // init boundaries of all global, non-m2m, GL including indexing types
  CellIndex<dim>::lowerBoundary = {0};
  CellIndex<dim>::upperBoundary = tarch::la::Vector<dim, int>{globalNumberCouplingCellsInclGL - tarch::la::Vector<dim, unsigned int>{1}};
  CellIndex<dim>::setDomainParameters();

  CellIndex<dim, IndexTrait::vector>::lowerBoundary = CellIndex<dim>::lowerBoundary;
  CellIndex<dim, IndexTrait::vector>::upperBoundary = CellIndex<dim>::upperBoundary;
  CellIndex<dim, IndexTrait::vector>::setDomainParameters();

  // init boundaries of all global, non-m2m, GL excluding indexing types
  CellIndex<dim, IndexTrait::noGhost>::lowerBoundary = {1};
  CellIndex<dim, IndexTrait::noGhost>::upperBoundary = tarch::la::Vector<dim, int>{globalNumberCouplingCellsInclGL - tarch::la::Vector<dim, unsigned int>{2}};
  CellIndex<dim, IndexTrait::noGhost>::setDomainParameters();

  CellIndex<dim, IndexTrait::vector, IndexTrait::noGhost>::lowerBoundary = CellIndex<dim, IndexTrait::noGhost>::lowerBoundary;
  CellIndex<dim, IndexTrait::vector, IndexTrait::noGhost>::upperBoundary = CellIndex<dim, IndexTrait::noGhost>::upperBoundary;
  CellIndex<dim, IndexTrait::vector, IndexTrait::noGhost>::setDomainParameters();

// CellIndex<dim> and CellIndex<dim, IndexTrait::vector> have been set up, so from here on it is ok to do some basic conversions
#if (COUPLING_MD_ERROR == COUPLING_MD_YES)
  _isInitialized = true;
#endif

  // init boundaries of all global, m2m, GL excluding indexing types

  CellIndex<dim, IndexTrait::md2macro, IndexTrait::noGhost>::lowerBoundary = BaseIndex<dim>{tarch::la::Vector<dim, int>{(int)(outerRegion + 1)}};
  CellIndex<dim, IndexTrait::md2macro, IndexTrait::noGhost>::upperBoundary =
      BaseIndex<dim>::upperBoundary - BaseIndex<dim>{tarch::la::Vector<dim, int>{(int)(outerRegion + 1)}};
  CellIndex<dim, IndexTrait::md2macro, IndexTrait::noGhost>::setDomainParameters();

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

  _numberProcesses = numberProcesses;

  _scalarNumberProcesses = _numberProcesses[0];
  for (unsigned int d = 1; d < dim; d++)
    _scalarNumberProcesses *= _numberProcesses[d];

  const unsigned int parallelTopologyOffset = (_rank / _scalarNumberProcesses) * _scalarNumberProcesses; // copied from IndexConversion
  _parallelTopology =
      coupling::paralleltopology::ParallelTopologyFactory::getParallelTopology<dim>(parallelTopologyType, _numberProcesses, parallelTopologyOffset);

  auto coords = _parallelTopology->getProcessCoordinates(_rank);
  tarch::la::Vector<3, int> boxMin, boxMax;
  // init boundaries of all local, non-m2m, GL including indexing types
  {
    for (unsigned int i = 0; i < dim; i++) {
      const auto backWeight = std::reduce(subdomainWeights[i].begin(), subdomainWeights[i].begin() + coords[i], 0u);
      const auto totalWeight = std::reduce(subdomainWeights[i].begin() + coords[i], subdomainWeights[i].end(), backWeight);
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
      std::cout << "Dim: " << i << " totalWeight: " << totalWeight << " backWeight: " << backWeight << " coords: " << coords[0] << ", " << coords[1] << ", "
                << coords[2] << std::endl;
#endif
      // calculate box bounds from cumulative weights of previous ranks, and the
      // weight of the current rank
      boxMin[i] = backWeight * globalNumberCouplingCells[i] / totalWeight;
      boxMax[i] = boxMin[i] + (subdomainWeights[i][coords[i]] * globalNumberCouplingCells[i] / totalWeight);
    }
    CellIndex<dim, IndexTrait::local>::lowerBoundary = BaseIndex<dim>{boxMin + tarch::la::Vector<dim, int>{1}};
    CellIndex<dim, IndexTrait::local>::upperBoundary = BaseIndex<dim>{boxMax + tarch::la::Vector<dim, int>{1}};
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
}

/*
 * This was in large parts stolen from IndexConversion.
 */
template <unsigned int dim>
std::vector<unsigned int> coupling::indexing::IndexingService<dim>::getRanksForGlobalIndex(const BaseIndex<dim>& globalCellIndex) const {

#if (COUPLING_MD_ERROR == COUPLING_MD_YES)
  if (!_isInitialized) {
    throw std::runtime_error(std::string("coupling::indexing::getRanksForGlobalIndex: IndexingService not initialized! "));
  }
#endif

  std::vector<unsigned int> ranks;
  // using the old meaning of 'globalNumberCouplingCells' from
  // IndexConversion
  const auto globalNumberCouplingCells = BaseIndex<dim>::numberCellsInDomain - tarch::la::Vector<dim, unsigned int>{2};

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
    end[d] = globalNumberCouplingCells[d] + 1;
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
        const unsigned int rank = getUniqueRankForCouplingCell(thisGlobalCellIndex, globalNumberCouplingCells);

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

/*
 * This was in large parts stolen from IndexConversion.
 * Note that this uses the globalNumberCouplingCells definition excl. the
 * ghost layer.
 */
template <unsigned int dim>
unsigned int
coupling::indexing::IndexingService<dim>::getUniqueRankForCouplingCell(tarch::la::Vector<dim, unsigned int> globalCellIndex,
                                                                       const tarch::la::Vector<dim, unsigned int>& globalNumberCouplingCells) const {
  // vector containing avg number of macro cells, not counting global GL.
  tarch::la::Vector<dim, unsigned int> averageLocalNumberCouplingCells{0};
  for (unsigned int d = 0; d < dim; d++) {
    if (globalCellIndex[d] >= globalNumberCouplingCells[d] + 2) { // greater or equal to the total global number incl GL (+2)
      using namespace std::string_literals;
      throw std::runtime_error("IndexingService: getUniqueRankForCouplingCell(): Global cell index greater than global size in dim "s + std::to_string(d));
    }
    if (globalNumberCouplingCells[d] % _numberProcesses[d] != 0) {
      std::stringstream ss;
      ss << "IndexingService: getUniqueRankForCouplingCell(): ERROR: Number "
            "of coupling cells must be divisible by number of processes! ";
      ss << "globalNumberCouplingCells = " << globalNumberCouplingCells;
      ss << ", numberProcesses = " << _numberProcesses;
      throw std::runtime_error(ss.str());
    }
    averageLocalNumberCouplingCells[d] = globalNumberCouplingCells[d] / _numberProcesses[d];
  }

  tarch::la::Vector<dim, unsigned int> processCoords(0);
  for (unsigned int d = 0; d < dim; d++) {
    // special case: cell in first section
    if (globalCellIndex[d] < averageLocalNumberCouplingCells[d] + 1) {
      processCoords[d] = 0;
      // special case: cell in last section
    } else if (globalCellIndex[d] > averageLocalNumberCouplingCells[d] * (_numberProcesses[d] - 1)) {
      processCoords[d] = _numberProcesses[d] - 1;
      // all other cases
    } else {
      // remove ghost layer contribution from vector index (...-1)
      processCoords[d] = (globalCellIndex[d] - 1) / averageLocalNumberCouplingCells[d];
    }
  }

  return _parallelTopology->getRank(processCoords);
}

template <unsigned int dim> unsigned int coupling::indexing::IndexingService<dim>::getUniqueRankForCouplingCell(const BaseIndex<dim>& globalCellIndex) const {

#if (COUPLING_MD_ERROR == COUPLING_MD_YES)
  if (!_isInitialized) {
    throw std::runtime_error(std::string("coupling::indexing::getUniqueRankForCouplingCell: IndexingService not initialized! "));
  }
#endif

  return getUniqueRankForCouplingCell((tarch::la::Vector<dim, unsigned int>)(globalCellIndex.get()), I09::numberCellsInDomain);
}

// declare specialisation of IndexingService
#ifdef INDEXING_ENABLE_DIM2
template class coupling::indexing::IndexingService<2>;
#endif
template class coupling::indexing::IndexingService<3>;
