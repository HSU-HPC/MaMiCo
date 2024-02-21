// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_INTERFACE_MACROSCOPICSOLVERINTERFACE_H_
#define _MOLECULARDYNAMICS_COUPLING_INTERFACE_MACROSCOPICSOLVERINTERFACE_H_

#include "coupling/CouplingMDDefinitions.h"
#include "tarch/la/Vector.h"
#include <vector>

namespace coupling {
namespace interface {
template <unsigned int dim> class MacroscopicSolverInterface;
}
} // namespace coupling

/** This class provides
 *	@brief interface for the macroscopic, i.e. continuum solver
 *	@tparam dim Number of dimensions; it can be 1, 2 or 3
 *  @author Philipp Neumann
 */
template <unsigned int dim> class coupling::interface::MacroscopicSolverInterface {

public:
  /** Costructor */
  MacroscopicSolverInterface() {}
  /** Destructor */
  virtual ~MacroscopicSolverInterface() {}

  /** This function specifies if the cell at position globalCellIndex shall be
   *received from the MD solver It does not send the cell, but only steers the
   *process.
   *	@param globalCellIndex
   *	@returns true if the cell at position globalCellIndex shall be received
   *from the MD solver, false otherwise.
   */
  virtual bool receiveMacroscopicQuantityFromMDSolver(tarch::la::Vector<dim, unsigned int> globalCellIndex) = 0;

  /** This function specifies if the cell at position globalCellIndex shall be
   *sent to the MD solver. It does not send the cell, but only steers the
   *process.
   *	@param globalCellIndex
   *	@returns true if the cell at position globalCellIndex shall be sent to
   *the MD solver, false otherwise.
   */
  virtual bool sendMacroscopicQuantityToMDSolver(tarch::la::Vector<dim, unsigned int> globalCellIndex) = 0;

  /** This function determines all the ranks on which the macroscopic solver
   *holds data of the coupling cell at index globalCellIndex. By default,
   *this method is used for send/recv operations. However, there are situations
   *where the target and source rank definitions (target rank=rank of
   *macroscopic solver that shall receive data from MD; source rank=rank of
   *macroscopic solver that shall send data to MD) may differ. For example,
   *assume a domain decomposition using ghost layers. You may want to send data
   *from cells inside a process-local domain and want to receive data in all
   *copies of a given coupling cell (incl. ghost cells), but you may not be
   *able to provide valid data on every process for each cell instance (e.g., if
   *this cell is part of a ghost layer). For this case, the method
   *getSourceRanks() and getTargetRanks() can be implemented accordingly.
   *	@param globalCellIndex
   *	@return all the ranks on which the macroscopic solver holds data of the
   *coupling cell at index globalCellIndex.
   */
  virtual std::vector<unsigned int> getRanks(tarch::la::Vector<dim, unsigned int> globalCellIndex) = 0;

  /** This function determines the source ranks for the global coupling cell
   *at index globalCellIndex. A "source rank" is defined as a rank of the
   *macroscopic solver that provides a valid copy of the cell at globalCellIndex
   *to the coupling tool. You may have multiple source ranks; however, each
   *source rank should hold a valid copy. Default: return all ranks that the
   *cell at globalCellIndex is associated to.
   *	@param globalCellIndex
   *	@return the source ranks for the global coupling cell at index
   *globalCellIndex.
   */
  virtual std::vector<unsigned int> getSourceRanks(tarch::la::Vector<dim, unsigned int> globalCellIndex) { return getRanks(globalCellIndex); }

  /** This function determines the target ranks for the global coupling cell
   *at index globalCellIndex. A "target rank" is defined as a rank of the
   *macroscopic solver that obtains a valid copy of the cell at globalCellIndex
   *from the coupling tool. You may have multiple target ranks. Default: return
   *all ranks that the cell at globalCellIndex is associated to.
   *	@param globalCellIndex
   *	@return the target ranks for the global coupling cell at index
   *globalCellIndex
   */
  virtual std::vector<unsigned int> getTargetRanks(tarch::la::Vector<dim, unsigned int> globalCellIndex) { return getRanks(globalCellIndex); }
};
#endif // _MOLECULARDYNAMICS_COUPLING_INTERFACE_MACROSCOPICSOLVERINTERFACE_H_
