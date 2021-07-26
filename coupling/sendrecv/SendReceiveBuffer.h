// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_SENDRECV_SENDRECEIVEBUFFER_H_
#define _MOLECULARDYNAMICS_COUPLING_SENDRECV_SENDRECEIVEBUFFER_H_

#include <string.h>
#include <map>
#include <set>
#include "tarch/la/Vector.h"
#include "coupling/IndexConversion.h"
#include "coupling/sendrecv/DataExchange.h"
#include "coupling/CouplingMDDefinitions.h"
#include "DataExchangeFromMD2Macro.h"

#if (COUPLING_MD_PARALLEL==COUPLING_MD_YES)
#include <mpi.h>
#endif


namespace coupling {
  namespace sendrecv {
    template<class MacroscopicCell,unsigned int dim>
    class SendReceiveBuffer;
  }
}



/** send-/ receive methodology including buffer implementations. The access to the buffers is prescribed by the DataExchange object.
 *  @author Philipp Neumann
 */
template<class MacroscopicCell,unsigned int dim>
class coupling::sendrecv::SendReceiveBuffer {
  public:
    SendReceiveBuffer();
    virtual ~SendReceiveBuffer();

  protected:
    /** deletes the buffers */
    void deleteBuffers();


    /** fills all information that needs to be sent from a macroscopic cell into the send-buffer. */
    void writeToSendBuffer(
      const coupling::IndexConversion<dim>& indexConversion,
      coupling::sendrecv::DataExchange<MacroscopicCell,dim> &dataExchange,
      const MacroscopicCell& cell,
      tarch::la::Vector<dim,unsigned int> globalVectorIndex
    );


    void writeToBcastBuffer(
      const coupling::IndexConversion<dim>& indexConversion,
      coupling::sendrecv::DataExchange<MacroscopicCell,dim> &dataExchange,
      const MacroscopicCell& cell,
      tarch::la::Vector<dim, unsigned  int> globalVectorIndex
    );


    void writeToReduceBuffer(
      const coupling::IndexConversion<dim>& indexConversion,
      coupling::sendrecv::DataExchange<MacroscopicCell,dim> &dataExchange,
      const MacroscopicCell& cell,
      tarch::la::Vector<dim, unsigned  int> globalVectorIndex
    );


    /** reads the information from the receive-buffer and fills it into a macroscopic cell. */
    void readFromReceiveBuffer(
      const coupling::IndexConversion<dim> &indexConversion,
      coupling::sendrecv::DataExchange<MacroscopicCell,dim> &dataExchange,
      MacroscopicCell &macroscopicCell,
      tarch::la::Vector<dim,unsigned int> globalVectorIndex
    );


    void readFromCollectiveBuffer(
      const coupling::IndexConversion<dim> &indexConversion,
      coupling::sendrecv::DataExchange<MacroscopicCell,dim> &dataExchange,
      MacroscopicCell &macroscopicCell,
      tarch::la::Vector<dim,unsigned int> globalVectorIndex
    );

    void readFromReduceBuffer(
      const coupling::IndexConversion<dim> & indexConversion,
      coupling::sendrecv::DataExchangeFromMD2Macro<dim> &dataExchange,
      MacroscopicCell & macroscopicCell,
      tarch::la::Vector<dim, unsigned int> globalVectorIndex
    );


    /** according to rule by dataExchange, the receive buffers are allocated. This function adds a contribution for the cell at globalVectorIndex.
     */
    void allocateReceiveBuffers(
      const coupling::IndexConversion<dim> &indexConversion,
      coupling::sendrecv::DataExchange<MacroscopicCell,dim> &dataExchange,
      tarch::la::Vector<dim,unsigned int> globalVectorIndex
    );

    void allocateBcastBufferForReceiving(
      const coupling::IndexConversion<dim> &indexConversion,
      coupling::sendrecv::DataExchange<MacroscopicCell,dim> &dataExchange,
      tarch::la::Vector<dim,unsigned int> globalVectorIndex
    );


    void allocateReduceBufferForReceiving(
      const coupling::IndexConversion<dim> &indexConversion,
      coupling::sendrecv::DataExchange<MacroscopicCell,dim> &dataExchange,
      tarch::la::Vector<dim,unsigned int> globalVectorIndex
    );


    /** triggers the MPI-sending on the respective buffers. No sending for information transfer from/ to this rank. */
    void triggerSending(
      const coupling::IndexConversion<dim>& indexConversion,
      coupling::sendrecv::DataExchange<MacroscopicCell,dim>& dataExchange
    );

    void triggerBcasts(
      unsigned int rank
    );

    /** triggers the MPI-receiving on the respective buffers. No receiving of information from/to this rank. */
    void triggerReceiving(
      const coupling::IndexConversion<dim>& indexConversion,
      coupling::sendrecv::DataExchange<MacroscopicCell,dim>& dataExchange
    );

    void triggerReduce(
      unsigned int rank
    );

    /** wait for all send and receive operations to complete. */
    void waitAllOperations(const coupling::IndexConversion<dim>& indexConversion);

    void waitAllBcasts(const coupling::IndexConversion<dim>& indexConversion);

    /** allocates send and receive requests */
    void allocateRequests(const coupling::IndexConversion<dim>& indexConversion);

    void allocateBcastRequests(const unsigned int thisRank);

  private:
    /** data structure for send- and receive-buffer. */
    struct BufferWithID {
      double *buffer;
      unsigned int bufferSize;

      BufferWithID(): buffer(NULL),bufferSize(0){}
    };


    struct BufferCollective {
        std::vector<double> buffer;
        std::set<unsigned int> nonRootRanks;
        std::set<unsigned int> cellIndices;
        unsigned int rootRank;

        BufferCollective(): buffer(), nonRootRanks(), cellIndices(), rootRank(-1){}
    };


    /** deletes everything inside a given buffer */
    void deleteBuffer(std::map<unsigned int,BufferWithID>& buffer);

    /** buffer for storing all received messages from MD. Each map entry is identified by a respective rank. */
    std::map<unsigned int,BufferWithID > _receiveBuffer;
    std::map<unsigned int,BufferWithID > _sendBuffer;


    /** members for collective communication */
    std::map<unsigned int, BufferCollective> _bcastBuffer;
    std::map<unsigned int, BufferCollective> _reduceBuffer;


    #if (COUPLING_MD_PARALLEL==COUPLING_MD_YES)
    bool _requestsAllocated;  /** flag that will always be reset after every send operation. Triggers instantiation of requests */
    MPI_Request *_requests;
    int _receiveSize;  /** number of receive requests */
    int _sendSize;     /** number of send requests */

    std::vector<MPI_Comm> _subComms;
    std::vector<MPI_Group> _subGroups;
    int _bcastSize;

    static void elementWiseSum(void * in, void * inout, int * len, MPI_Datatype *datatype) {
      auto * output = (double *) inout;
      auto * input = (double *) in;
      for(int i = 0; i < *len; ++i) {
        output[i] += input[i];
      }
    }

    MPI_Op elementWiseSumOperation;
    #endif
};

#include "SendReceiveBuffer.cpph"


#endif // _MOLECULARDYNAMICS_COUPLING_SENDRECV_SENDRECEIVEBUFFER_H_
