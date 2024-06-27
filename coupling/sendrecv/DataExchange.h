// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_SENDRECV_DATAEXCHANGE_H_
#define _MOLECULARDYNAMICS_COUPLING_SENDRECV_DATAEXCHANGE_H_

#include "tarch/la/Vector.h"
#include "coupling/datastructures/CouplingCell.h"

namespace coupling {
namespace sendrecv {
template <unsigned int dim> class DataExchange;
}
} // namespace coupling

/** This class holds information on the data exchange between processes. In
 *contrast to the SendReceiveBuffer which solely handles the technical part of
 *sending and receiving, this class needs to guarantee the validity in
 *communication, e.g. make sure that each source rank matches to a target rank
 *etc.
 *	@brief generic class for the the data exchange purposes.
 *	@tparam CouplingCell cell type
 *	@tparam dim Number of dimensions; it can be 1, 2 or 3
 *  @author Philipp Neumann
 */
template <unsigned int dim> class coupling::sendrecv::DataExchange {
public:
  /** Constructor: assign an tag (_tag) to the DataExchange.
   * @param tag
   */
  DataExchange(unsigned int tag) : _tag(tag) {}
  /** Destructor */
  virtual ~DataExchange() {}

  /** returns the tag associated to messages for this DataExchange. @return
   * _tag*/
  unsigned int getTag() const { return _tag; }

  /** returns the ranks to which a particular cell (at index idx)
   * should be sent. This method should globally define the target ranks for
   * each global cell index of a coupling cell.
   * 	@param idx unique global cell index
   */
  virtual std::vector<unsigned int> getTargetRanks(I01 idx) = 0;

  /** returns all ranks from which a particular cell (at index idx)
   * is sent. This method should globally define the source ranks for each
   * global cell index of a coupling cell.
   * @param idx unique global cell index
   */
  virtual std::vector<unsigned int> getSourceRanks(I01 idx) = 0;

  /** local rule to read from a coupling cell and write data to (e.g. send)
   * buffer
   * 	@param buffer
   * 	@param cell
   */
  virtual void readFromCell(double* const buffer, const coupling::datastructures::CouplingCell<dim>& cell) = 0;

  /** local rule to read from receive buffer and write data to coupling cell
   * 	@param buffer
   * 	@param cell
   */
  virtual void writeToCell(const double* const buffer, const coupling::datastructures::CouplingCell<dim>& cell) = 0;

  /** returns the number of doubles that are sent per coupling cell. */
  virtual unsigned int getDoublesPerCell() const = 0;

private:
  /** tag to uniquely determine these kinds of messages that are sent between
   * processes. */
  const unsigned int _tag;
};

#endif // _MOLECULARDYNAMICS_COUPLING_SENDRECV_DATAEXCHANGE_H_
