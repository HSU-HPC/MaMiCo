// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_SENDRECV_MD2MACRO_H_
#define _MOLECULARDYNAMICS_COUPLING_SENDRECV_MD2MACRO_H_

#include "coupling/CouplingMDDefinitions.h"
#include "coupling/datastructures/CellContainer.h"
#include "coupling/datastructures/FlexibleCellContainer.h"
#include "coupling/sendrecv/DataExchangeFromMD2Macro.h"
#include "coupling/sendrecv/SendReceiveBuffer.h"
#include <vector>

namespace coupling {
namespace sendrecv {
template <class Cell_T, unsigned int dim> class FromMD2Macro;
}
} // namespace coupling

/** extends the generic SendReceiveBuffer for send coupling cell information
 *from MaMiCo to the macroscopic solver.
 *	@brief sends coupling cell information from MaMiCo to the macroscopic
 *solver. Derived from the class coupling::sendrecv::SendReceiveBuffer
 *	@tparam Cell_T cell type
 *	@tparam dim Number of dimensions; it can be 1, 2 or 3
 *  @author Philipp Neumann
 */
template <class Cell_T, unsigned int dim> class coupling::sendrecv::FromMD2Macro : public coupling::sendrecv::SendReceiveBuffer<Cell_T, dim> {

public:
  /** Constructor, just calling the constructor of the
   * coupling::sendrecv::SendReceiveBuffer  */
  FromMD2Macro() : coupling::sendrecv::SendReceiveBuffer<Cell_T, dim>() {}
  /** Destructor */
  virtual ~FromMD2Macro() {}

  using Macro_Container_T = coupling::datastructures::FlexibleCellContainer<dim>;
  using Local_Container_T = coupling::datastructures::CellContainer<I02, dim>;

  /** sends information from the local coupling cells of MaMiCo (only inner
   * non-ghost cells of this process) to the macroscopic solver. Since the
   * macroscopic solver can have an arbitrary distribution of cells on the
   * processes, the buffer for receiving the cell data is provided to this
   * function in terms of an array of coupling cells including the respective
   * global cell indices. This function basically calls for these purposes
   * sendFromMD2MacroNonBlocking(...) and wait4SendFromMD2Macro(...) in a row.
   * 	@param dataExchange
   * 	@param src
   *  @param dst
   */
  void sendFromMD2Macro(coupling::sendrecv::DataExchange<Cell_T, dim>& dataExchange, const Local_Container_T& src, const Macro_Container_T& dst);

  /** reduces information from the local coupling cells of MaMiCo (only inner
   * non-ghost cells of this process) and sends it to the macroscopic solver.
   * 	@param dataExchange
   * 	@param src
   *  @param dst
   */
  template <class Container_T1, class Container_T2>
  void reduceFromMD2Macro(std::vector<coupling::sendrecv::DataExchangeFromMD2Macro<dim>*>& dataExchange, const Container_T1& src, const Container_T2& dst);

  /** triggers the send/recv operations for data transfer. After returning,
   * these data transfers do not necessarily need to be finished, according to
   * ISend/IRecv in MPI.
   * 	@param dataExchange
   * 	@param src
   *  @param dst
   */
  void sendFromMD2MacroNonBlocking(coupling::sendrecv::DataExchange<Cell_T, dim>& dataExchange, const Local_Container_T& src, const Macro_Container_T& dst);

  /** waits for the send operation--instantiated by
   * sendFromMD2MacroNonBlocking(...)--to be finished and writes the data to
   * couplingCellsFromMacroscopicSolver.
   * 	@param dataExchange
   * 	@param cells
   */
  void wait4SendFromMD2Macro(coupling::sendrecv::DataExchange<Cell_T, dim>& dataExchange, const Macro_Container_T& cells);

private:
  /** loops over the whole local Cartesian grid (only non-ghost cells!) and
   * writes respective cells to send buffer. For each cell,
   * writeToSendBuffer(...) of SendReceiveBuffer is used.
   * 	@param dataExchange
   * 	@param cells
   */
  void writeToSendBuffer(coupling::sendrecv::DataExchange<Cell_T, dim>& dataExchange, const Local_Container_T& cells);

  /** loops over the whole local Cartesian grid (only non-ghost cells!) and
   * writes respective cells to reduce buffer.
   * 	@param dataExchange
   * 	@param cells
   */
  void writeToReduceBuffer(coupling::sendrecv::DataExchangeFromMD2Macro<dim>& dataExchange, const Local_Container_T& cells);

  /** allocates the receive buffers for the macroscopic solver. Since we do not
   * know anything about the macroscopic solver, we only have a list of global
   * vector cell indices available for possible cells to be received on this
   * rank. For each global cell which is in the list,
   * allocateReceiveBuffers(...) of SendReceiveBuffer is called.
   * 	@param dataExchange
   * 	@param cells
   */
  void allocateReceiveBuffers(coupling::sendrecv::DataExchange<Cell_T, dim>& dataExchange, const Macro_Container_T& cells);

  /** allocates the reduce receive buffers for the macroscopic solver. Since we do not
   * know anything about the macroscopic solver, we only have a list of global
   * vector cell indices available for possible cells to be received on this
   * rank. For each global cell which is in the list,
   * allocateReduceBufferForReceiving(...) of SendReceiveBuffer is called.
   * 	@param dataExchange
   * 	@param cells
   */
  template <class Container_T> void allocateReduceBufferForReceiving(coupling::sendrecv::DataExchangeFromMD2Macro<dim>& dataExchange, const Container_T& cells);

  /** reads information from the receive buffer and stores the result in the
   * list of coupling cells. Since this is a receive for the coupling
   * cells on the side of the macroscopic solver, we just have a list of global
   * cell indices and corresponding coupling cell buffers. For each cell in
   * this list, readFromReceiveBuffer(...) of SendReceiveBuffer is called.
   * 	@param dataExchange
   * 	@param cells
   */
  void readFromReceiveBuffer(coupling::sendrecv::DataExchange<Cell_T, dim>& dataExchange, const Macro_Container_T& cells);

  /** reads information from the reduce buffer and stores the result in the
   * list of coupling cells. Since this is a receive for the coupling
   * cells on the side of the macroscopic solver, we just have a list of global
   * cell indices and corresponding coupling cell buffers. For each cell in
   * this list, readFromReduceBuffer(...) of SendReceiveBuffer is called.
   * 	@param dataExchange
   * 	@param cells
   */
  template <class Container_T> void readFromReduceBuffer(coupling::sendrecv::DataExchangeFromMD2Macro<dim>& dataExchange, const Container_T& cells);
};

#include "FromMD2Macro.cpph"

#endif // _MOLECULARDYNAMICS_COUPLING_SENDRECV_MD2MACRO_H_
