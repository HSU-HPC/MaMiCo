// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_SENDRECV_MD2MACRORECVONLY_H_
#define _MOLECULARDYNAMICS_COUPLING_SENDRECV_MD2MACRORECVONLY_H_

#include "coupling/CouplingMDDefinitions.h"
#include "coupling/sendrecv/SendReceiveBuffer.h"
#include <vector>

namespace coupling {
namespace sendrecv {
template <class Cell_T, unsigned int dim> class FromMD2MacroRecvOnly;
}
} // namespace coupling

/** sends coupling cell information from MaMiCo to the macroscopic solver.
 *This version expects to only receive data from, but will not send any data to
 *a macroscopic solver process. This is used for multi-MD simulations where the
 *macroscopic solver process may run on a completely different process than a
 *single MD simulation.
 *	@brief sends coupling cell information from MaMiCo to the macroscopic
 *solver. Expects to only receive data from, but will not send any data. Derived
 *from the class coupling::sendrecv::SendReceiveBuffer
 *	@tparam Cell_T cell type
 *	@tparam dim Number of dimensions; it can be 1, 2 or 3
 *  @author Philipp Neumann
 */
template <class Cell_T, unsigned int dim> class coupling::sendrecv::FromMD2MacroRecvOnly : public coupling::sendrecv::SendReceiveBuffer<Cell_T, dim> {

public:
  /** Constructor, just calling the constructor of the
   * coupling::sendrecv::SendReceiveBuffer  */
  FromMD2MacroRecvOnly() : coupling::sendrecv::SendReceiveBuffer<Cell_T, dim>() {}
  /** Destructor */
  virtual ~FromMD2MacroRecvOnly() {}

  /** sends information from the local coupling cells of MaMiCo (only inner
   * non-ghost cells of this process) to the macroscopic solver. Since the
   * macroscopic solver can have an arbitrary distribution of cells on the
   * processes, the buffer for receiving the cell data is provided to this
   * function in terms of an array of coupling cells including the respective
   * global cell indices. Basically calls sendFromMD2MacroNonBlocking(...) and
   * wait4SendFromMD2Macro(...) in a row.
   * 	@param dataExchange
   * 	@param cells
   */
  template <class Container_T> void sendFromMD2Macro(coupling::sendrecv::DataExchange<Cell_T, dim>& dataExchange, const Container_T& cells);

  /** triggers the send/recv operations for data transfer. After returning,
   *these data transfers do not necessarily need to be finished, according to
   *ISend/IRecv in MPI.
   * 	@param dataExchange
   * 	@param cells
   *	@sa  class FromMD2Macro. */
  template <class Container_T> void sendFromMD2MacroNonBlocking(coupling::sendrecv::DataExchange<Cell_T, dim>& dataExchange, const Container_T& cells);

  /** waits for the send operation--instantiated by
   *sendFromMD2MacroNonBlocking(...)--to be finished and writes the data to
   *couplingCellsFromMacroscopicSolver.
   * 	@param dataExchange
   * 	@param cells
   *	@sa class FromMD2Macro */
  template <class Container_T> void wait4SendFromMD2Macro(coupling::sendrecv::DataExchange<Cell_T, dim>& dataExchange, const Container_T& cells);

private:
  /** allocates the receive buffers for the macroscopic solver. Since we do not
   * know anything about the macroscopic solver, we only have a list of global
   * vector cell indices available for possible cells to be received on this
   * rank. For each global cell which is in the list,
   * allocateReceiveBuffers(...) of SendReceiveBuffer is called.
   * 	@param dataExchange
   *  @param cells
   */
  template <class Container_T> void allocateReceiveBuffers(coupling::sendrecv::DataExchange<Cell_T, dim>& dataExchange, const Container_T& cells);

  /** reads information from the receive buffer and stores the result in the
   * list of coupling cells. Since this is a receive for the coupling
   * cells on the side of the macroscopic solver, we just have a list of global
   * cell indices and corresponding macrocsopic cell buffers. For each cell in
   * this list, readFromReceiveBuffer(...) of SendReceiveBuffer is called.
   * 	@param dataExchange
   * 	@param cells
   */
  template <class Container_T> void readFromReceiveBuffer(coupling::sendrecv::DataExchange<Cell_T, dim>& dataExchange, const Container_T& cells);
};

#include "FromMD2MacroRecvOnly.cpph"

#endif // _MOLECULARDYNAMICS_COUPLING_SENDRECV_MD2MACRORECVONLY_H_
