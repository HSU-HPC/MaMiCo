// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_SENDRECV_FROMMACRO2MD_H_
#define _MOLECULARDYNAMICS_COUPLING_SENDRECV_FROMMACRO2MD_H_

#include "coupling/sendrecv/SendReceiveBuffer.h"
#include <vector>

namespace coupling {
  namespace sendrecv {
    template<class MacroscopicCell, unsigned int dim>
    class FromMacro2MD;
  }
}


/** extends the generic SendReceiveBuffer for transfer of quantities from a macroscopic solver to
 *  the macroscopic cells of the coupling tool (incl. ghost cells).
 *  @author Philipp Neumann
 */
template<class MacroscopicCell,unsigned int dim>
class coupling::sendrecv::FromMacro2MD: public coupling::sendrecv::SendReceiveBuffer<MacroscopicCell,dim> {

  public:
    FromMacro2MD(): coupling::sendrecv::SendReceiveBuffer<MacroscopicCell,dim>(){}
    virtual ~FromMacro2MD(){}

    /** sends information from the local macroscopic cells of a macroscopic solver to the macroscopic
     *  cells of MaMico (ghost cells are also allowed). Since the macroscopic solver can have an arbitrary distribution of cells on the processes,
     *  the buffer for sending the cell data is provided to this function in terms of an array of macroscopic cells including the
     *  respective global cell indices.
     *  This function calls sendFromMacro2MDNonBlocking(..) and immediately subsequently wait4SendFromMacro2MD(...). Those two methods may alternatively be used, e.g., for "non-blocking" communication.
     */
    void sendFromMacro2MD(
      const coupling::IndexConversion<dim> &indexConversion,
      coupling::sendrecv::DataExchange<MacroscopicCell,dim> &dataExchange,
      const std::vector<MacroscopicCell*> &macroscopicCellsFromMamico,
      const std::vector<MacroscopicCell*> &macroscopicCellsFromMacroscopicSolver,
      const unsigned int * const globalCellIndicesFromMacroscopicSolver
    );

    void sendFromMacro2MDCollective(
      const coupling::IndexConversion<dim> &indexConversion,
      coupling::sendrecv::DataExchange<MacroscopicCell,dim> &dataExchange,
      const std::vector<MacroscopicCell*> &macroscopicCellsFromMacroscopicSolver,
      const unsigned int * const globalCellIndicesFromMacroscopicSolver
    );

    /** sends data from macro to MD. After returning, the data transfer may not be completely finished, similar to a IRecv/ISend-call by MPI. Please use wait4SendFromMacro2MD(...) to guarantee
     *  that the data transfer has been finished.
     */
    void sendFromMacro2MDNonBlocking(
      const coupling::IndexConversion<dim> &indexConversion,
      coupling::sendrecv::DataExchange<MacroscopicCell,dim> &dataExchange,
      const std::vector<MacroscopicCell*> &macroscopicCellsFromMacroscopicSolver,
      const unsigned int * const globalCellIndicesFromMacroscopicSolver
    );

    /** waits for the data transfer--instantiated by sendFromMacro2MDNonBlocking(..)--to be finished and fills the information into the macroscopic cells from Mamico. */
    void wait4SendFromMacro2MD(
      const coupling::IndexConversion<dim> &indexConversion,
      coupling::sendrecv::DataExchange<MacroscopicCell,dim> &dataExchange,
      const std::vector<MacroscopicCell*> &macroscopicCellsFromMamico
    );

  private:
    /** given a list of macroscopic cells (from the macroscopic solver), the data from these cells are written
     *  to the send buffer.
     */
    void writeToSendBuffer(
      const coupling::IndexConversion<dim>& indexConversion,
      coupling::sendrecv::DataExchange<MacroscopicCell,dim> &dataExchange,
      const std::vector<MacroscopicCell*> &macroscopicCells,
      const unsigned int * const globalCellIndices
    );

    void writeToSendBufferCollective(
      const coupling::IndexConversion<dim>& indexConversion,
      coupling::sendrecv::DataExchange<MacroscopicCell,dim> &dataExchange,
      const std::vector<MacroscopicCell*> &macroscopicCells,
      const unsigned int * const globalCellIndices
    );


    /** allocates the receive buffers for the macroscopic solver. Since we want to obtain data on the side of MaMiCo,
     *  we can just loop over all local macroscopic cells of the coupling tool and call allocateReceiveBuffers(...) of the
     *  SendReceiveBuffer for each respective cell.
     */
    void allocateReceiveBuffers(
      const coupling::IndexConversion<dim> &indexConversion,
      coupling::sendrecv::DataExchange<MacroscopicCell,dim> &dataExchange
    );

    void allocateReceiveBuffersCollective(
      const coupling::IndexConversion<dim> &indexConversion,
      coupling::sendrecv::DataExchange<MacroscopicCell,dim> &dataExchange
    );

    /** reads information from the receive buffer and stores the result in the macroscopic cells. Since this is a receive
     *  for the macroscopic cells on the side of the coupling tool, we can consider a Cartesian grid topology for the
     *  macroscopic cells. For each cell, readFromReceiveBuffer(...) of SendReceiveBuffer is called.
     */
    void readFromReceiveBuffer(
      const coupling::IndexConversion<dim> &indexConversion,
      coupling::sendrecv::DataExchange<MacroscopicCell,dim> &dataExchange,
      const std::vector<MacroscopicCell *> &macroscopicCells
    );
};

#include "FromMacro2MD.cpph"

#endif // _MOLECULARDYNAMICS_COUPLING_SENDRECV_FROMMACRO2MD_H_
