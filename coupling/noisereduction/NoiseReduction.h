#ifndef _MOLECULARDYNAMICS_COUPLING_NOISEREDUCTION_NOISEREDUCTION_H_
#define _MOLECULARDYNAMICS_COUPLING_NOISEREDUCTION_NOISEREDUCTION_H_

#include "coupling/datastructures/MacroscopicCell.h"
#include "coupling/IndexConversion.h"
#include "tarch/utils/MultiMDService.h"

namespace coupling {
namespace noisereduction { template <unsigned int dim> class NoiseReduction; }
}

/**
 *  Interface for noise reduction algorithm, that is for smoothing operations to
 * be applied on data from MD
 *  before it is sent to the macroscopic solver.
 *  @author Piet Jarmatz
 */
template <unsigned int dim> class coupling::noisereduction::NoiseReduction {
public:
  NoiseReduction(const coupling::IndexConversion<dim> &indexConversion,
                 const tarch::utils::MultiMDService<dim> &multiMDService,
                 const bool doubleTraversal = false)
      : _doubleTraversal(doubleTraversal), _indexConversion(indexConversion),
        _multiMDService(multiMDService) {}
  virtual ~NoiseReduction() {}

  /** Is called for every macroscopic cell right before sending the
   * macroscopicMass and -Momentum
   *  data to the macroscopic solver, but after the call to
   * TransferStrategy::processInnerMacroscopicCellBeforeSendingMDSolverData
   *  This method is only applied to macroscopic cells that cover parts of the
   * MD domain; it is not applied
   *  in the outer macroscopic cells.
   */
  virtual void processInnerMacroscopicCell(
      coupling::datastructures::MacroscopicCell<dim> &cell,
      const unsigned int &index) {}
  virtual void beginProcessInnerMacroscopicCells() {}
  virtual void endProcessInnerMacroscopicCells() {}

  /** Is called for every macroscopic cell right before sending the
   * macroscopicMass and -Momentum
   *  data to the macroscopic solver, but after the call to
   * TransferStrategy::processOuterMacroscopicCellBeforeSendingMDSolverData
   *  This method is only applied to outer macroscopic cells, that is cells that
   * are located outside the MD domain.
   */
  virtual void processOuterMacroscopicCell(
      coupling::datastructures::MacroscopicCell<dim> &cell,
      const unsigned int &index) {}
  virtual void beginProcessOuterMacroscopicCells() {}
  virtual void endProcessOuterMacroscopicCells() {}

  bool _doubleTraversal;

protected:
  const coupling::IndexConversion<dim> &_indexConversion;
  const tarch::utils::MultiMDService<dim> &_multiMDService;
};
#endif // _MOLECULARDYNAMICS_COUPLING_NOISEREDUCTION_NOISEREDUCTION_H_
