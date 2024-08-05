// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_SENDRECV_SENDRECEIVEBUFFER_H_
#define _MOLECULARDYNAMICS_COUPLING_SENDRECV_SENDRECEIVEBUFFER_H_

#include "DataExchangeFromMD2Macro.h"
#include "coupling/CouplingMDDefinitions.h"
#include "coupling/indexing/IndexingService.h"
#include "coupling/sendrecv/DataExchange.h"
#include "tarch/la/Vector.h"
#include <algorithm>
#include <map>
#include <set>
#include <string.h>
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
#include <mpi.h>
#endif

namespace coupling {
namespace sendrecv {
template <class Cell_T, unsigned int dim> class SendReceiveBuffer;
}
} // namespace coupling

/** generic class for send-/ receive methodology including buffer
 *implementations. The access to the buffers is prescribed by the DataExchange
 *object.
 *	@brief generic class for send-/ receive methodology.
 *	@tparam CouplingCell cell type
 *	@tparam dim Number of dimensions; it can be 1, 2 or 3
 *  @author Philipp Neumann
 */
template <class Cell_T, unsigned int dim> class coupling::sendrecv::SendReceiveBuffer {
public:
  /** Constructor */
  SendReceiveBuffer();
  /** Destructor */
  virtual ~SendReceiveBuffer();

protected:
  /** @brief deletes the buffers */
  void deleteBuffers();

  /** @brief fills all information that needs to be sent from a coupling cell
   * into the send-buffer.
   * 	@param dataExchange
   * 	@param cells
   */
  template<class Container_T> void writeToSendBuffer(coupling::sendrecv::DataExchange<Cell_T, dim>& dataExchange, const Container_T& cells);

  /** @brief fills all information that needs to be broadcast from a coupling cell
   * into the broadcast-buffer.
   * 	@param dataExchange
   * 	@param cell
   * 	@param idx
   */
  void writeToBcastBuffer(coupling::sendrecv::DataExchange<Cell_T, dim>& dataExchange, const Cell_T& cell, I01 idx);

  /** @brief fills all information that needs to be reduced to a coupling cell
   * into the reduce-buffer.
   * 	@param dataExchange
   * 	@param cell
   * 	@param idx
   */
  void writeToReduceBuffer(coupling::sendrecv::DataExchange<Cell_T, dim>& dataExchange, const Cell_T& cell, I01 idx);

  /** reads the information from the receive-buffer and fills it into a
   * coupling cell.
   * 	@param dataExchange
   * 	@param couplingCell
   * 	@param idx
   */
  template<class Container_T> void readFromReceiveBuffer(coupling::sendrecv::DataExchange<Cell_T, dim>& dataExchange, const Container_T& cells);

  void readFromCollectiveBuffer(coupling::sendrecv::DataExchange<Cell_T, dim>& dataExchange, Cell_T& couplingCell, I01 idx);

  /** reads the information from the reduce-buffer and fills it into a
   * coupling cell.
   * 	@param dataExchange
   * 	@param couplingCell
   * 	@param idx
   */
  void readFromReduceBuffer(coupling::sendrecv::DataExchangeFromMD2Macro<dim>& dataExchange, Cell_T& couplingCell, I01 idx);

  /** according to rule by dataExchange, the receive buffers are allocated. This
   * function adds a contribution for the cell at idx.
   * 	@param dataExchange
   * 	@param idx
   */
  void allocateReceiveBuffers(coupling::sendrecv::DataExchange<Cell_T, dim>& dataExchange, I01 idx);

  /** Allocates buffer for receiving in the context of the broadcast operation
   * 	@param dataExchange
   * 	@param idx
   */
  void allocateBcastBufferForReceiving(coupling::sendrecv::DataExchange<Cell_T, dim>& dataExchange, I01 idx);

  /** Allocates buffer for receiving in the context of the reduce operation
   * 	@param dataExchange
   * 	@param idx
   */
  void allocateReduceBufferForReceiving(coupling::sendrecv::DataExchange<Cell_T, dim>& dataExchange, I01 idx);

  /** triggers the MPI-sending on the respective buffers. No sending for
   * information transfer from/ to this rank.
   * 	@param dataExchange
   */
  void triggerSending(coupling::sendrecv::DataExchange<Cell_T, dim>& dataExchange);

  /** triggers the MPI-broadcast on the respective buffers.
   * 	@param rank
   */
  void triggerBcasts(unsigned int rank);

  /** triggers the MPI-receiving on the respective buffers. No receiving of
   * information from/to this rank.
   * 	@param dataExchange
   */
  void triggerReceiving(coupling::sendrecv::DataExchange<Cell_T, dim>& dataExchange);

  /** triggers the MPI-reduce on the respective buffers.
   * 	@param rank
   */
  void triggerReduce(unsigned int rank);

  /** wait for all send and receive operations to complete.
   */
  void waitAllOperations();

  /** wait for all broadcast or reduce operations to complete */
  void waitAllCollectiveOperations();

  /** allocates send and receive requests
   */
  void allocateRequests();

  /** allocates broadcast request
   * 	@param thisRank
   */
  void allocateBcastRequests(unsigned int thisRank);

  /** allocates reduce request
   * 	@param thisRank
   */
  void allocateReduceRequests(unsigned int thisRank);

private:
  /** data structure for send- and receive-buffer. */
  struct BufferWithID {
    double* buffer;
    unsigned int bufferSize;

    BufferWithID() : buffer(NULL), bufferSize(0) {}
  };

  struct BufferCollective {
    std::vector<double> buffer;
    std::set<unsigned int> nonRootRanks;
    std::set<unsigned int> cellIndices;
    unsigned int rootRank;

    BufferCollective() : buffer(), nonRootRanks(), cellIndices(), rootRank(-1) {}
  };

  /** deletes everything inside a given buffer
   * 	@param buffer
   */
  void deleteBuffer(std::map<unsigned int, BufferWithID>& buffer);

  /** buffer for storing all received messages from MD. Each map entry is
   * identified by a respective rank. */
  std::map<unsigned int, BufferWithID> _receiveBuffer;
  std::map<unsigned int, BufferWithID> _sendBuffer;

  /** members for collective communication */
  std::map<unsigned int, BufferCollective> _bcastBuffer;
  std::map<unsigned int, BufferCollective> _reduceBuffer;

#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
  bool _requestsAllocated; /** flag that will always be reset after every send
                              operation. Triggers instantiation of requests */
  MPI_Request* _requests;
  int _receiveSize; /** number of receive requests */
  int _sendSize;    /** number of send requests */

  std::vector<MPI_Comm> _subComms;
  std::vector<MPI_Group> _subGroups;
  int _bcastOrReduceSize;

  static void elementWiseSum(void* in, void* inout, int* len, MPI_Datatype* datatype) {
    auto* output = (double*)inout;
    auto* input = (double*)in;
    for (int i = 0; i < *len; ++i) {
      output[i] += input[i];
    }
  }

  MPI_Op elementWiseSumOperation;
#endif
};

#include "SendReceiveBuffer.cpph"

#endif // _MOLECULARDYNAMICS_COUPLING_SENDRECV_SENDRECEIVEBUFFER_H_
