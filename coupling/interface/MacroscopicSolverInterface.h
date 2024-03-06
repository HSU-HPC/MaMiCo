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

  /* This function defines an offset of cells which is considered to be the outer region.
   * It replaces the legacy functions 'receiveMacroscopicQuantityFromMDSolver' and 'sendMacroscopicQuantityToMDSolver'
   * from older versions of MaMiCo. Data are sent to MD solver for cells that are not in the ghost layer
   * and not part of the inner region. Data needs so be send from micro to macro solver for all inner cells.
   *  @return number of cells in outer region per direction at each boundary, e.g. 3
   */
  virtual unsigned int getOuterRegion() = 0;

  /** This function determines all the ranks on which the macroscopic solver
   *holds data of the coupling cell at index idx. By default,
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
   *	@param idx
   *	@return all the ranks on which the macroscopic solver holds data of the
   *coupling cell at index idx.
   */
  virtual std::vector<unsigned int> getRanks(I01 idx) = 0;

  /** This function determines the source ranks for the global coupling cell
   *at index idx. A "source rank" is defined as a rank of the
   *macroscopic solver that provides a valid copy of the cell at idx
   *to the coupling tool. You may have multiple source ranks; however, each
   *source rank should hold a valid copy. Default: return all ranks that the
   *cell at idx is associated to.
   *	@param idx
   *	@return the source ranks for the global coupling cell at index
   *idx.
   */
  virtual std::vector<unsigned int> getSourceRanks(I01 idx) { return getRanks(idx); }

  /** This function determines the target ranks for the global coupling cell
   *at index idx. A "target rank" is defined as a rank of the
   *macroscopic solver that obtains a valid copy of the cell at idx
   *from the coupling tool. You may have multiple target ranks. Default: return
   *all ranks that the cell at idx is associated to.
   *	@param idx
   *	@return the target ranks for the global coupling cell at index
   *idx
   */
  virtual std::vector<unsigned int> getTargetRanks(I01 idx) { return getRanks(idx); }
};
#endif // _MOLECULARDYNAMICS_COUPLING_INTERFACE_MACROSCOPICSOLVERINTERFACE_H_
