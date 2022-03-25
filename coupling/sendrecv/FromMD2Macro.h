// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_SENDRECV_MD2MACRO_H_
#define _MOLECULARDYNAMICS_COUPLING_SENDRECV_MD2MACRO_H_

#include "coupling/CouplingMDDefinitions.h"
#include "coupling/sendrecv/SendReceiveBuffer.h"
#include <vector>

namespace coupling {
namespace sendrecv {
template <class MacroscopicCell, unsigned int dim> class FromMD2Macro;
}
} // namespace coupling

/** extends the generic SendReceiveBuffer for send macroscopic cell information
 *from MaMiCo to the macroscopic solver.
 *	@brief sends macroscopic cell information from MaMiCo to the macroscopic
 *solver. Derived from the class coupling::sendrecv::SendReceiveBuffer
 *	@tparam MacroscopicCell cell type
 *	@tparam dim Number of dimensions; it can be 1, 2 or 3
 *  @author Philipp Neumann
 */
template <class MacroscopicCell, unsigned int dim> class coupling::sendrecv::FromMD2Macro : public coupling::sendrecv::SendReceiveBuffer<MacroscopicCell, dim> {

public:
  /** Constructor, just calling the constructor of the
   * coupling::sendrecv::SendReceiveBuffer  */
  FromMD2Macro() : coupling::sendrecv::SendReceiveBuffer<MacroscopicCell, dim>() {}
  /** Destructor */
  virtual ~FromMD2Macro() {}

  /** sends information from the local macroscopic cells of MaMiCo (only inner
   * non-ghost cells of this process) to the macroscopic solver. Since the
   * macroscopic solver can have an arbitrary distribution of cells on the
   * processes, the buffer for receiving the cell data is provided to this
   * function in terms of an array of macroscopic cells including the respective
   * global cell indices. This function basically calls for these purposes
   * sendFromMD2MacroNonBlocking(...) and wait4SendFromMD2Macro(...) in a row.
   * 	@param indexConversion
   * 	@param dataExchange
   * 	@param macroscopicCellsFromMamico
   * 	@param macroscopicCellsFromMacroscopicSolver
   * 	@param globalCellIndicesFromMacroscopicSolver
   */
  void sendFromMD2Macro(const coupling::IndexConversion<dim> &indexConversion, coupling::sendrecv::DataExchange<MacroscopicCell, dim> &dataExchange,
                        const std::vector<MacroscopicCell *> &macroscopicCellsFromMamico,
                        const std::vector<MacroscopicCell *> &macroscopicCellsFromMacroscopicSolver,
                        const unsigned int *const globalCellIndicesFromMacroscopicSolver);

  /** triggers the send/recv operations for data transfer. After returning,
   * these data transfers do not necessarily need to be finished, according to
   * ISend/IRecv in MPI.
   * 	@param indexConversion
   * 	@param dataExchange
   * 	@param macroscopicCellsFromMamico
   * 	@param macroscopicCellsFromMacroscopicSolver
   * 	@param globalCellIndicesFromMacroscopicSolver
   */
  void sendFromMD2MacroNonBlocking(const coupling::IndexConversion<dim> &indexConversion, coupling::sendrecv::DataExchange<MacroscopicCell, dim> &dataExchange,
                                   const std::vector<MacroscopicCell *> &macroscopicCellsFromMamico,
                                   const std::vector<MacroscopicCell *> &macroscopicCellsFromMacroscopicSolver,
                                   const unsigned int *const globalCellIndicesFromMacroscopicSolver);

  /** waits for the send operation--instantiated by
   * sendFromMD2MacroNonBlocking(...)--to be finished and writes the data to
   * macroscopicCellsFromMacroscopicSolver.
   * 	@param indexConversion
   * 	@param dataExchange
   * 	@param macroscopicCellsFromMacroscopicSolver
   * 	@param globalCellIndicesFromMacroscopicSolver
   */
  void wait4SendFromMD2Macro(const coupling::IndexConversion<dim> &indexConversion, coupling::sendrecv::DataExchange<MacroscopicCell, dim> &dataExchange,
                             const std::vector<MacroscopicCell *> &macroscopicCellsFromMacroscopicSolver,
                             const unsigned int *const globalCellIndicesFromMacroscopicSolver);

private:
  /** loops over the whole local Cartesian grid (only non-ghost cells!) and
   * writes respective cells to send buffer. For each cell,
   * writeToSendBuffer(...) of SendReceiveBuffer is used.
   * 	@param indexConversion
   * 	@param dataExchange
   * 	@param macroscopicCells
   */
  void writeToSendBuffer(const coupling::IndexConversion<dim> &indexConversion, coupling::sendrecv::DataExchange<MacroscopicCell, dim> &dataExchange,
                         const std::vector<MacroscopicCell *> &macroscopicCells);

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

#include "FromMD2Macro.cpph"

#endif // _MOLECULARDYNAMICS_COUPLING_SENDRECV_MD2MACRO_H_
