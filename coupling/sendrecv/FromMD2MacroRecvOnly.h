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
template <class MacroscopicCell, unsigned int dim> class FromMD2MacroRecvOnly;
}
} // namespace coupling

/** sends macroscopic cell information from MaMiCo to the macroscopic solver.
 *This version expects to only receive data from, but will not send any data to
 *a macroscopic solver process. This is used for multi-MD simulations where the
 *macroscopic solver process may run on a completely different process than a
 *single MD simulation.
 *	@brief sends macroscopic cell information from MaMiCo to the macroscopic
 *solver. Expects to only receive data from, but will not send any data. Derived
 *from the class coupling::sendrecv::SendReceiveBuffer
 *	@tparam MacroscopicCell cell type
 *	@tparam dim Number of dimensions; it can be 1, 2 or 3
 *  @author Philipp Neumann
 */
template <class MacroscopicCell, unsigned int dim>
class coupling::sendrecv::FromMD2MacroRecvOnly : public coupling::sendrecv::SendReceiveBuffer<MacroscopicCell, dim> {

public:
  /** Constructor, just calling the constructor of the
   * coupling::sendrecv::SendReceiveBuffer  */
  FromMD2MacroRecvOnly() : coupling::sendrecv::SendReceiveBuffer<MacroscopicCell, dim>() {}
  /** Destructor */
  virtual ~FromMD2MacroRecvOnly() {}

  /** sends information from the local macroscopic cells of MaMiCo (only inner
   * non-ghost cells of this process) to the macroscopic solver. Since the
   * macroscopic solver can have an arbitrary distribution of cells on the
   * processes, the buffer for receiving the cell data is provided to this
   * function in terms of an array of macroscopic cells including the respective
   * global cell indices. Basically calls sendFromMD2MacroNonBlocking(...) and
   * wait4SendFromMD2Macro(...) in a row.
   * 	@param indexConversion
   * 	@param dataExchange
   * 	@param macroscopicCellsFromMacroscopicSolver
   * 	@param globalCellIndicesFromMacroscopicSolver
   */
  void sendFromMD2Macro(const coupling::IndexConversion<dim> &indexConversion, coupling::sendrecv::DataExchange<MacroscopicCell, dim> &dataExchange,
                        const std::vector<MacroscopicCell *> &macroscopicCellsFromMacroscopicSolver,
                        const unsigned int *const globalCellIndicesFromMacroscopicSolver);

  /** triggers the send/recv operations for data transfer. After returning,
   *these data transfers do not necessarily need to be finished, according to
   *ISend/IRecv in MPI.
   * 	@param indexConversion
   * 	@param dataExchange
   * 	@param macroscopicCellsFromMacroscopicSolver
   * 	@param globalCellIndicesFromMacroscopicSolver
   *	@sa  class FromMD2Macro. */
  void sendFromMD2MacroNonBlocking(const coupling::IndexConversion<dim> &indexConversion, coupling::sendrecv::DataExchange<MacroscopicCell, dim> &dataExchange,
                                   const std::vector<MacroscopicCell *> &macroscopicCellsFromMacroscopicSolver,
                                   const unsigned int *const globalCellIndicesFromMacroscopicSolver);

  /** waits for the send operation--instantiated by
   *sendFromMD2MacroNonBlocking(...)--to be finished and writes the data to
   *macroscopicCellsFromMacroscopicSolver.
   * 	@param indexConversion
   * 	@param dataExchange
   * 	@param macroscopicCellsFromMacroscopicSolver
   * 	@param globalCellIndicesFromMacroscopicSolver
   *	@sa class FromMD2Macro */
  void wait4SendFromMD2Macro(const coupling::IndexConversion<dim> &indexConversion, coupling::sendrecv::DataExchange<MacroscopicCell, dim> &dataExchange,
                             const std::vector<MacroscopicCell *> &macroscopicCellsFromMacroscopicSolver,
                             const unsigned int *const globalCellIndicesFromMacroscopicSolver);

private:
  /** allocates the receive buffers for the macroscopic solver. Since we do not
   * know anything about the macroscopic solver, we only have a list of global
   * vector cell indices available for possible cells to be received on this
   * rank. For each global cell which is in the list,
   * allocateReceiveBuffers(...) of SendReceiveBuffer is called.
   * 	@param indexConversion
   * 	@param dataExchange
   * 	@param globalCellIndices
   * 	@param numberCells
   */
  void allocateReceiveBuffers(const coupling::IndexConversion<dim> &indexConversion, coupling::sendrecv::DataExchange<MacroscopicCell, dim> &dataExchange,
                              const unsigned int *const globalCellIndices, unsigned int numberCells);

  /** reads information from the receive buffer and stores the result in the
   * list of macroscopic cells. Since this is a receive for the macroscopic
   * cells on the side of the macroscopic solver, we just have a list of global
   * cell indices and corresponding macrocsopic cell buffers. For each cell in
   * this list, readFromReceiveBuffer(...) of SendReceiveBuffer is called.
   * 	@param indexConversion
   * 	@param dataExchange
   * 	@param macroscopicCells
   * 	@param globalCellIndices
   */
  void readFromReceiveBuffer(const coupling::IndexConversion<dim> &indexConversion, coupling::sendrecv::DataExchange<MacroscopicCell, dim> &dataExchange,
                             const std::vector<MacroscopicCell *> &macroscopicCells, const unsigned int *const globalCellIndices);
};

#include "FromMD2MacroRecvOnly.cpph"

#endif // _MOLECULARDYNAMICS_COUPLING_SENDRECV_MD2MACRORECVONLY_H_
