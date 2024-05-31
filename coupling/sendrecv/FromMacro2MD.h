// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_SENDRECV_FROMMACRO2MD_H_
#define _MOLECULARDYNAMICS_COUPLING_SENDRECV_FROMMACRO2MD_H_

#include "coupling/sendrecv/DataExchangeFromMacro2MD.h"
#include "coupling/sendrecv/SendReceiveBuffer.h"
#include <vector>

namespace coupling {
namespace sendrecv {
template <class CouplingCell, unsigned int dim> class FromMacro2MD;
}
} // namespace coupling

/** extends the generic SendReceiveBuffer for transfer of quantities from a
 *macroscopic solver to the coupling cells of the coupling tool (incl. ghost
 *cells).
 *	@brief SendReceiveBuffer for transfer of quantities from a macroscopic
 *solver to the coupling cells. Derived from the class
 *coupling::sendrecv::SendReceiveBuffer
 *	@tparam CouplingCell cell type
 *	@tparam dim Number of dimensions; it can be 1, 2 or 3
 *  @author Philipp Neumann
 */
template <class CellContainerT, class CouplingCell, unsigned int dim> class coupling::sendrecv::FromMacro2MD : public coupling::sendrecv::SendReceiveBuffer<CouplingCell, dim> {

public:
  /** Constructor, just calling the constructor of the
   * coupling::sendrecv::SendReceiveBuffer  */
  FromMacro2MD() : coupling::sendrecv::SendReceiveBuffer<CellContainerT, CouplingCell, dim>() {}
  /** Destructor */
  virtual ~FromMacro2MD() {}

  /** sends information from the local coupling cells of a macroscopic solver
   * to the coupling cells of MaMico (ghost cells are also allowed). Since
   * the macroscopic solver can have an arbitrary distribution of cells on the
   * processes, the buffer for sending the cell data is provided to this
   * function in terms of an array of coupling cells including the respective
   * global cell indices. This function calls sendFromMacro2MDNonBlocking(..)
   * and immediately subsequently wait4SendFromMacro2MD(...). Those two methods
   * may alternatively be used, e.g., for "non-blocking" communication.
   * 	@param dataExchange
   * 	@param couplingCellsFromMamico
   * 	@param couplingCellsFromMacroscopicSolver
   * 	@param globalCellIndicesFromMacroscopicSolver
   */
  void sendFromMacro2MD(coupling::sendrecv::DataExchange<CouplingCell, dim>& dataExchange, const coupling::datastructures::CellContainer<dim> couplingCellContainer,
                        const coupling::datastructures::FlexibleCellContainer<dim>& couplingCellContainerFromMacroscopicSolver);

  void bcastFromMacro2MD(std::vector<coupling::sendrecv::DataExchangeFromMacro2MD<dim>*>& dataExchangeFromCouplingCellServices,
                         const std::vector<CouplingCell*>& couplingCellsFromMacroscopicSolver, const I00* const globalCellIndicesFromMacroscopicSolver,
                         std::vector<std::vector<CouplingCell*>> allCouplingCellsFromCouplingCellServices);

  /** sends data from macro to MD. After returning, the data transfer may not be
   * completely finished, similar to a IRecv/ISend-call by MPI. Please use
   * wait4SendFromMacro2MD(...) to guarantee that the data transfer has been
   * finished.
   * 	@param dataExchange
   * 	@param couplingCellsFromMacroscopicSolver
   * 	@param globalCellIndicesFromMacroscopicSolver
   */
  void sendFromMacro2MDNonBlocking(coupling::sendrecv::DataExchange<CouplingCell, dim>& dataExchange,
                                   const std::vector<CouplingCell*>& couplingCellsFromMacroscopicSolver,
                                   const I00* const globalCellIndicesFromMacroscopicSolver);

  /** waits for the data transfer--instantiated by
   * sendFromMacro2MDNonBlocking(..)--to be finished and fills the information
   * into the coupling cells from Mamico.
   * 	@param dataExchange
   * 	@param couplingCellsFromMamico
   */
  void wait4SendFromMacro2MD(coupling::sendrecv::DataExchange<CouplingCell, dim>& dataExchange, const std::vector<CouplingCell*>& couplingCellsFromMamico);

private:
  /** given a list of coupling cells (from the macroscopic solver), the data
   * from these cells are written to the send buffer.
   * 	@param dataExchange
   * 	@param couplingCells
   * 	@param globalCellIndices
   */
  void writeToSendBuffer(coupling::sendrecv::DataExchange<CouplingCell, dim>& dataExchange, const std::vector<CouplingCell*>& couplingCells,
                         const I00* const globalCellIndices);

  void writeToSendBufferCollective(coupling::sendrecv::DataExchange<CouplingCell, dim>& dataExchange, const std::vector<CouplingCell*>& couplingCells,
                                   const I00* const globalCellIndices);

  /** allocates the receive buffers for the macroscopic solver. Since we want to
   * obtain data on the side of MaMiCo, we can just loop over all local
   * coupling cells of the coupling tool and call allocateReceiveBuffers(...)
   * of the SendReceiveBuffer for each respective cell.
   * 	@param dataExchange
   */
  void allocateReceiveBuffers(coupling::sendrecv::DataExchange<CouplingCell, dim>& dataExchange);

  void allocateReceiveBuffersCollective(coupling::sendrecv::DataExchange<CouplingCell, dim>& dataExchange);

  /** reads information from the receive buffer and stores the result in the
   * coupling cells. Since this is a receive for the coupling cells on the
   * side of the coupling tool, we can consider a Cartesian grid topology for
   * the coupling cells. For each cell, readFromReceiveBuffer(...) of
   * SendReceiveBuffer is called.
   * 	@param dataExchange
   * 	@param couplingCells
   */
  void readFromReceiveBuffer(coupling::sendrecv::DataExchange<CouplingCell, dim>& dataExchange, const std::vector<CouplingCell*>& couplingCells);

  void readFromCollectiveBuffer(coupling::sendrecv::DataExchange<CouplingCell, dim>& dataExchange, const std::vector<CouplingCell*>& couplingCells);
};

#include "FromMacro2MD.cpph"

#endif // _MOLECULARDYNAMICS_COUPLING_SENDRECV_FROMMACRO2MD_H_
