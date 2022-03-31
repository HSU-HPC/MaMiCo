// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_SENDRECV_DATAEXCHANGEFROMMACRO2MD_H_
#define _MOLECULARDYNAMICS_COUPLING_SENDRECV_DATAEXCHANGEFROMMACRO2MD_H_

#include "coupling/CouplingMDDefinitions.h"
#include "coupling/datastructures/MacroscopicCell.h"
#include "coupling/interface/MacroscopicSolverInterface.h"
#include "coupling/sendrecv/DataExchange.h"

namespace coupling {
namespace sendrecv {
template <unsigned int dim> class DataExchangeFromMacro2MD;
}
} // namespace coupling

/** data exchange from the macroscopic solver to the MD solver, that is to the
 * coupling tool. We transfer the buffers
 *  microscopicMass and microscopicMomentum of the MacroscopicCell.
 *  The target ranks are determined by the getRanksForMacroscopicCell() method
 * of IndexConversion: since macroscopic cells
 *  may exist on different ranks (due to ghost layers), this method determines
 * all ranks with a particular macroscopic cell
 *  and returns a vector with all required ranks (see also documentation of
 * IndexConversion). The source ranks
 *  are determined via the macroscopic solver interface's method getRanks()
 *	@brief data exchange from the macroscopic solver to the MD solver.
 *Derived from the class coupling::sendrecv::DataExchange
 *	@tparam dim Number of dimensions; it can be 1, 2 or 3
 *	@sa  DataExchangeFromMD2Macro
 *  @author Philipp Neumann
 */
template <unsigned int dim>
class coupling::sendrecv::DataExchangeFromMacro2MD : public coupling::sendrecv::DataExchange<coupling::datastructures::MacroscopicCell<dim>, dim> {

public:
  /** Constructor
   * 	@param interface macroscopic solver interface
   * 	@param indexConversion index conversion
   * 	@param tagoffset 0 per default
   */
  DataExchangeFromMacro2MD(coupling::interface::MacroscopicSolverInterface<dim> *interface, coupling::IndexConversion<dim> *indexConversion,
                           unsigned int tagoffset = 0)
      : coupling::sendrecv::DataExchange<coupling::datastructures::MacroscopicCell<dim>, dim>(TAG_FROM_MACRO2MD + tagoffset), _interface(interface),
        _indexConversion(indexConversion) {
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
    std::cout << "DataExchangeFromMacro2MD initialised..." << std::endl;
#endif
  }
  /** Destructor */
  virtual ~DataExchangeFromMacro2MD() {}

  /** returns the ranks to which a particular cell (at index globalCellIndex)
should be sent.
         * 	@param  globalCellIndex
         *	@return the corresponding ranks via IndexConversion, if we need
information on MD side, otherwise empty vector
         */
  virtual std::vector<unsigned int> getTargetRanks(tarch::la::Vector<dim, unsigned int> globalCellIndex) {
    // if we need information on MD side, return the respective ranks via
    // IndexConversion
    if (_interface->sendMacroscopicQuantityToMDSolver(globalCellIndex)) {
      return _indexConversion->getRanksForMacroscopicCell(globalCellIndex);
      // otherwise return empty vector
    }
    return std::vector<unsigned int>();
  }

  /** returns all ranks from which a particular cell (at index globalCellIndex)
is sent.
         * 	@param  globalCellIndex
         *	@return the corresponding ranks via MacroscopicSolverInterface,
if we need information on MD side, otherwise empty vector
         */
  virtual std::vector<unsigned int> getSourceRanks(tarch::la::Vector<dim, unsigned int> globalCellIndex) {
    if (_interface->sendMacroscopicQuantityToMDSolver(globalCellIndex)) {
      return _interface->getSourceRanks(globalCellIndex);
    }
    return std::vector<unsigned int>();
  }

  /** local rule to read from macroscopic cell and write data to (e.g. send)
buffer.
     *  We only send the macroscopic mass and macroscopic momentum from MD to
the macroscopic solver.
         * 	@param buffer
         * 	@param cell
     */
  virtual void readFromCell(double *const buffer, const coupling::datastructures::MacroscopicCell<dim> &cell) {
    buffer[0] = cell.getMicroscopicMass();
    for (unsigned int d = 0; d < dim; d++) {
      buffer[d + 1] = cell.getMicroscopicMomentum()[d];
    }
  }

  /** local rule to read from receive buffer and write data to macroscopic cell
   * 	@param buffer
   * 	@param cell
   */
  virtual void writeToCell(const double *const buffer, coupling::datastructures::MacroscopicCell<dim> &cell) {
    tarch::la::Vector<dim, double> microscopicMomentum(0.0);
    for (unsigned int d = 0; d < dim; d++) {
      microscopicMomentum[d] = buffer[1 + d];
    }
    cell.setMicroscopicMomentum(microscopicMomentum);
    cell.setMicroscopicMass(buffer[0]);
  }

  /** returns the number of doubles that are sent per macroscopic cell. @return
   * 1+dim  */
  virtual unsigned int getDoublesPerCell() const {
    // 1 double: microscopic mass; dim doubles: microscopic momentum
    return 1 + dim;
  }

  void setIndexConversion(coupling::IndexConversion<dim> *indexConversion) { _indexConversion = indexConversion; }

private:
  coupling::interface::MacroscopicSolverInterface<dim> *_interface;
  coupling::IndexConversion<dim> *_indexConversion;
};
#endif // _MOLECULARDYNAMICS_COUPLING_SENDRECV_DATAEXCHANGEFROMMACRO2MD_H_
