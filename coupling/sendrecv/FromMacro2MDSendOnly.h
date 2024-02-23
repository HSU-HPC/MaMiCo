// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_SENDRECV_FROMMACRO2MDSENDONLY_H_
#define _MOLECULARDYNAMICS_COUPLING_SENDRECV_FROMMACRO2MDSENDONLY_H_

#include "coupling/sendrecv/SendReceiveBuffer.h"
#include <vector>

namespace coupling {
namespace sendrecv {
template <class CouplingCell, unsigned int dim> class FromMacro2MDSendOnly;
}
} // namespace coupling

/** extends the generic SendReceiveBuffer for transfer of quantities from a
 *macroscopic solver to the coupling cells of the coupling tool (incl. ghost
 *cells). This version only sends data but expects not to receive any data, e.g.
 *due to MD running on different ranks. Derived from the class
 *coupling::sendrecv::SendReceiveBuffer
 *	@brief SendReceiveBuffer for transfer of quantities from a macroscopic
 *solver to the coupling cells. Only sends data but expects not to receive
 *any data. Derived from the class coupling::sendrecv::SendReceiveBuffer
 *	@tparam CouplingCell cell type
 *	@tparam dim Number of dimensions; it can be 1, 2 or 3
 *  @author Philipp Neumann
 */
template <class CouplingCell, unsigned int dim>
class coupling::sendrecv::FromMacro2MDSendOnly : public coupling::sendrecv::SendReceiveBuffer<CouplingCell, dim> {

public:
  /** Constructor, just calling the constructor of the
   * coupling::sendrecv::SendReceiveBuffer  */
  FromMacro2MDSendOnly() : coupling::sendrecv::SendReceiveBuffer<CouplingCell, dim>() {}
  /** Destructor */
  virtual ~FromMacro2MDSendOnly() {}

  /** sends information from the local coupling cells of a macroscopic solver
   * to the coupling cells of MaMico (ghost cells are also allowed). Since
   * the macroscopic solver can have an arbitrary distribution of cells on the
   * processes, the buffer for sending the cell data is provided to this
   * function in terms of an array of coupling cells including the respective
   * global cell indices. Basically calls sendFromMacro2MDNonBlocking(...) and
   * wait4SendFromMacro2MD(...) in a row.
   * 	@param dataExchange
   * 	@param couplingCellsFromMacroscopicSolver
   * 	@param globalCellIndicesFromMacroscopicSolver
   */
  void sendFromMacro2MD(coupling::sendrecv::DataExchange<CouplingCell, dim>& dataExchange, const std::vector<CouplingCell*>& couplingCellsFromMacroscopicSolver,
                        const I00* const globalCellIndicesFromMacroscopicSolver);

  /** sends data from macro to MD. After returning, the data transfer may not be
   *completely finished, similar to a IRecv/ISend-call by MPI. Please use
   *wait4SendFromMacro2MD(...) to guarantee that the data transfer has been
   *finished.
   * 	@param dataExchange
   * 	@param couplingCellsFromMacroscopicSolver
   * 	@param globalCellIndicesFromMacroscopicSolver
   *	@sa  descriptions of same function of class FromMacro2MD.
   */
  void sendFromMacro2MDNonBlocking(coupling::sendrecv::DataExchange<CouplingCell, dim>& dataExchange,
                                   const std::vector<CouplingCell*>& couplingCellsFromMacroscopicSolver,
                                   const I00* const globalCellIndicesFromMacroscopicSolver);

  /** waits for the data transfer--instantiated by
   *sendFromMacro2MDNonBlocking(..)--to be finished and fills the information
   *into the coupling cells from Mamico.
   *	@sa class FromMacro2MD.
   */
  void wait4SendFromMacro2MD();

protected:
  /** given a list of coupling cells (from the macroscopic solver), the data
   * from these cells are written to the send buffer.
   * 	@param dataExchange
   * 	@param couplingCells
   * 	@param globalCellIndices
   */
  void writeToSendBuffer(coupling::sendrecv::DataExchange<CouplingCell, dim>& dataExchange, const std::vector<CouplingCell*>& couplingCells,
                         const I00* const globalCellIndices);
};

#include "FromMacro2MDSendOnly.cpph"

#endif // _MOLECULARDYNAMICS_COUPLING_SENDRECV_FROMMACRO2MDSENDONLY_H_
