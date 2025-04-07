// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_SENDRECV_DATAEXCHANGEFROMALLMACRO2MD_H_
#define _MOLECULARDYNAMICS_COUPLING_SENDRECV_DATAEXCHANGEFROMALLMACRO2MD_H_

#include "coupling/CouplingMDDefinitions.h"
#include "coupling/datastructures/CouplingCell.h"
#include "coupling/indexing/IndexingService.h"
#include "coupling/interface/MacroscopicSolverInterface.h"
#include "coupling/sendrecv/DataExchange.h"

namespace coupling {
namespace sendrecv {
template <unsigned int dim> class DataExchangeFromAllMacro2MD;
}
} // namespace coupling

/** Same as _MOLECULARDYNAMICS_COUPLING_SENDRECV_DATAEXCHANGEFROMMACRO2MD_H_ but also sends on inner cells
 */
template <unsigned int dim>
class coupling::sendrecv::DataExchangeFromAllMacro2MD : public coupling::sendrecv::DataExchange<coupling::datastructures::CouplingCell<dim>, dim> {

public:
  /** Constructor
   * 	@param interface macroscopic solver interface
   * 	@param tagoffset 0 per default
   */
  DataExchangeFromAllMacro2MD(coupling::interface::MacroscopicSolverInterface<dim>* interface, unsigned int topologyOffset, unsigned int tagoffset = 0)
      : coupling::sendrecv::DataExchange<coupling::datastructures::CouplingCell<dim>, dim>(TAG_FROM_MACRO2MD + tagoffset), _msi(interface),
        _topologyOffset(topologyOffset) {
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
    std::cout << "DataExchangeFromAllMacro2MD initialised..." << std::endl;
#endif
  }
  /** Destructor */
  virtual ~DataExchangeFromAllMacro2MD() {}

  /** returns the ranks to which a particular cell (at index idx)
   *should be sent.
   * 	@param  idx
   *	@return the corresponding ranks, if we need
   *information on MD side, otherwise empty vector
   */
  std::vector<unsigned int> getTargetRanks(I01 idx) override {
    // if we need information on MD side, return the respective ranks via
    if (I08::contains(idx)) {
      return IDXS.getRanksForGlobalIndex(idx, _topologyOffset);
      // otherwise return empty vector
    }
    return std::vector<unsigned int>();
  }

  /** returns all ranks from which a particular cell (at index idx)
   *is sent.
   * 	@param  idx
   *	@return the corresponding ranks via MacroscopicSolverInterface, if we
   *need information on MD side, otherwise empty vector
   */
  std::vector<unsigned int> getSourceRanks(I01 idx) override {
    if (I08::contains(idx)) {
      return _msi->getSourceRanks(idx);
    }
    return std::vector<unsigned int>();
  }

  /** local rule to read from coupling cell and write data to (e.g. send)
   * buffer. We only send the macroscopic mass and macroscopic momentum from MD
   * to the macroscopic solver.
   * 	@param buffer
   * 	@param cell
   */
  void readFromCell(double* const buffer, const coupling::datastructures::CouplingCell<dim>& cell) override {
    buffer[0] = cell.getMicroscopicMass();
    for (unsigned int d = 0; d < dim; d++) {
      buffer[d + 1] = cell.getMicroscopicMomentum()[d];
    }
  }

  /** local rule to read from receive buffer and write data to coupling cell
   * 	@param buffer
   * 	@param cell
   */
  void writeToCell(const double* const buffer, coupling::datastructures::CouplingCell<dim>& cell) override {
    tarch::la::Vector<dim, double> microscopicMomentum(0.0);
    for (unsigned int d = 0; d < dim; d++) {
      microscopicMomentum[d] = buffer[1 + d];
    }
    cell.setMicroscopicMomentum(microscopicMomentum);
    cell.setMicroscopicMass(buffer[0]);
  }

  /** returns the number of doubles that are sent per coupling cell. @return
   * 1+dim  */
  unsigned int getDoublesPerCell() const override {
    // 1 double: microscopic mass; dim doubles: microscopic momentum
    return 1 + dim;
  }

private:
  coupling::interface::MacroscopicSolverInterface<dim>* _msi;
  unsigned int _topologyOffset;
};
#endif // _MOLECULARDYNAMICS_COUPLING_SENDRECV_DATAEXCHANGEFROMALLMACRO2MD_H_
