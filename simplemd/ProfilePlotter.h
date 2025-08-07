// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_PROFILEPLOTTER_H_
#define _MOLECULARDYNAMICS_PROFILEPLOTTER_H_

#include "simplemd/MolecularDynamicsUserInput.h"
#include "simplemd/cell-mappings/ProfilePlotterMapping.h"
#include "simplemd/configurations/ProfilePlotterConfiguration.h"
#include "simplemd/services/ParallelTopologyService.h"
#include "tarch/la/Vector.h"
#include <vector>

namespace simplemd {
class ProfilePlotter;
}

/** computes and evaluates flow profiles, such as density,momentum.
 *  @author Philipp Neumann
 */
class simplemd::ProfilePlotter {
public:
  ProfilePlotter(const std::vector<simplemd::configurations::ProfilePlotterConfiguration>& configurations,
                 const simplemd::services::ParallelTopologyService& parallelTopologyService, simplemd::MoleculeContainer& moleculeContainer,
                 const double& linkedCellVolume, const unsigned int& localMDSimulation);
  ~ProfilePlotter();

  /** accumulates information from the respective linked cells and - in case that this is the respective timestep -
   *  plots the information to a file.
   */
  void accumulateAndPlotInformation(const unsigned int& t);

private:
  simplemd::MoleculeContainer& _moleculeContainer;
#if (MD_PARALLEL == MD_YES)
  /** in parallel, the _startCell and _range need to be adjusted to the intersection of the global profile position
   *  with the local domain. If a profile does not intersect the local domain, it is marked as inactive
   *  via the _isValid variable.
   */
  void adjustProfilesInParallel(const simplemd::services::ParallelTopologyService& parallelTopologyService);
#endif

  std::vector<simplemd::cellmappings::ProfilePlotterMapping*> _plotters;
  std::vector<tarch::la::Vector<MD_DIM, unsigned int>> _startCells;
  std::vector<tarch::la::Vector<MD_DIM, unsigned int>> _ranges;

#if (MD_PARALLEL == MD_YES)
  /** in parallel, some profiles may not intersect the local domain of a process.
   *  In that case, the ProfilePlotterMapping will be allocated, but the respective entry in this vector will
   *  be used to indicate that the profile is not active.
   */
  std::vector<bool> _isActive;
#endif
};
#endif // _MOLECULARDYNAMICS_PROFILEPLOTTER_H_
