// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_SERVICES_PARALLELANDLOCALBUFFERSERVICE_H_
#define _MOLECULARDYNAMICS_SERVICES_PARALLELANDLOCALBUFFERSERVICE_H_

#include "simplemd/MolecularDynamicsDefinitions.h"
#include "simplemd/Molecule.h"
#include <cmath>
#include <cstdlib>
#include <iostream>

/** Class for managing buffers for sending messages between processors and
 * "local" buffer for managing periodic boundary conditions, which need
 *  to be handled locally
 *
 * @author Nikola Tchipev
 */
namespace simplemd {
namespace services {
class ParallelAndLocalBufferService {
public:
  // We want a class for managing buffers, but we don't want
  // that class to exist out of the ParallelAndLocalBufferService
  // so we use a nested class

  // ========================= //
  // nested class SimpleBuffer //
  // ========================= //

  /** Class that mimics an stl vector, but is guaranteed to be
   * contiguous in memory, a requirement for MPI calls.
   *
   * @author Nikola Tchipev
   */
  class SimpleBuffer {
  public:
    SimpleBuffer();
    ~SimpleBuffer();

    /** allocate storage */
    bool initialise(const unsigned int doublesPerMolecule, const unsigned int upperBoundOnNumberOfMolecules);

    void shutdown();

    /** push data (molecule) into local or send buffer
     * if _capacity was not exceeded by pushing
     * or if reallocation was permitted (for local buffers),
     * function returns true
     * otherwise - false
     */
    bool pushData(const tarch::la::Vector<MD_DIM, double> position, const tarch::la::Vector<MD_DIM, double> velocity,
                  const tarch::la::Vector<MD_DIM, double> force, const double isFixed, const bool permitReallocation);

    /** clear buffer values for use at next iteration
     * capacity is preserved
     */
    void clearBuffer() { _length = 0; };

    void setLength(const unsigned int len) { _length = len; }

    unsigned int getLength() const { return _length; }
    unsigned int getCapacity() const { return _capacity; }
    double* getValues() const { return _values; }

  private:
    /* METHODS */

    /** double the capacity */
    bool reallocate();

    /* FIELDS */

    /** the values */
    double* _values;

    /** the number of values currently stored */
    unsigned int _length;

    /** allocated capacity */
    unsigned int _capacity;
  };

  // ========================= //
  // end of class SimpleBuffer //
  // ========================= //

public:
  ParallelAndLocalBufferService() {}
  ~ParallelAndLocalBufferService() {}

  /** allocate all necessary buffers
   * return false if allocation of a buffer fails
   */
  bool initialise(const unsigned int numUniqueNeighbours, const unsigned int numCellsPerBuffer[], const double avMoleculesPerCell);

  void shutdown();

  bool pushMoleculeToLocalBuffer(const tarch::la::Vector<MD_DIM, double>& position, const Molecule* mol);

  unsigned int getLocalBufferLength(unsigned int i_buf) const { return _localBuffer.getLength(); }

  SimpleBuffer* getLocalBuffer() { return &_localBuffer; }

#if (MD_PARALLEL == MD_YES)
  bool pushMoleculeToSendBuffer(const tarch::la::Vector<MD_DIM, double>& position, const Molecule* mol, const unsigned int i_buffer);

  unsigned int getBufferCapacity(unsigned int i_buf) const { return _sendBuffers[i_buf].getCapacity(); }
  unsigned int getSendBufferLength(unsigned int i_buf) const { return _sendBuffers[i_buf].getLength(); }
  SimpleBuffer* getSendBuffer(const unsigned int& i_buf) { return &(_sendBuffers[i_buf]); }
  SimpleBuffer* getReceiveBuffer(const unsigned int& i_buf) { return &(_receiveBuffers[i_buf]); }

  void setReceiveBufferLength(const unsigned int i_buf, const unsigned int count) { _receiveBuffers[i_buf].setLength(count); };
#endif

private:
  /* Methods */

  /** method that calculates upper bound on number of molecules to be stored in the respective buffer.
   * Determined from average number of molecules per linked cell, number of cells to be transferred
   * via respective buffer and a special function A:
   * upperBound[bufferIndex] = ceil(numCells[bufferIndex] * avMolPerCell * A(numCells[bufferIndex]))
   *
   * The purpose of this function A is to deal with the following:
   * when we send only one linked cell to a process, the probability that the average number of molecules is exceeded is
   * high, whereas, when we send thousands of cells, the probability that the average number of molecules is exceeded in all
   * of them is much lower. We also need to keep buffer size low,  so we choose the function A to satisfy the following:
   * 1. A: N -> R
   * 2. A(1) is "high", in our choice, A(1) = 5.5
   * 3. lim_{k to Inf} A(k) is "low", in our choice lim_{k to Inf} A(k) = 1.5
   * 4. A - monotonically decreasing
   *
   */
  unsigned int computeBufferUpperBound(const unsigned int numCells, const double avMoleculesPerCell) const;

  /* Fields */

  /**
   * buffer needed for handling local periodic boundaries in
   * broadcast-broadcast methods. Reallocation is permitted.
   */
  SimpleBuffer _localBuffer;

#if (MD_PARALLEL == MD_YES)
  unsigned int _numberActiveParallelBuffers;

  /** buffers for sending (only _numberActiveParallelBuffers will be used)
   * Reallocation is not permitted. It means the upper bound on number
   * of molecules was exceeded, causing termination.
   */
  SimpleBuffer _sendBuffers[MD_LINKED_CELL_NEIGHBOURS];

  /** buffers for receiving (only _numberActiveParallelBuffers will be used)
   * Reallocation is not permitted. It means the upper bound on number
   * of molecules was exceeded, causing termination.
   */
  SimpleBuffer _receiveBuffers[MD_LINKED_CELL_NEIGHBOURS];
#endif
};
} // namespace services
} // namespace simplemd

#endif
