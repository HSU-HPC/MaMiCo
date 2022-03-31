// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_MOLECULEMAPPINGS_WRITECHECKPOINTMAPPING_H_
#define _MOLECULARDYNAMICS_MOLECULEMAPPINGS_WRITECHECKPOINTMAPPING_H_

#include "simplemd/Molecule.h"
#include "simplemd/services/ParallelTopologyService.h"
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

namespace simplemd {
namespace moleculemappings { class WriteCheckPointMapping; }
}

/** writes checkpoint data for e.g. restarts.
 *
 *  @author Philipp Neumann
 */
class simplemd::moleculemappings::WriteCheckPointMapping {
public:
  WriteCheckPointMapping(const simplemd::services::ParallelTopologyService &
                             parallelTopologyService,
                         const std::string &filestem, const unsigned int &t)
      : _parallelTopologyService(parallelTopologyService), _file(NULL),
        _filestem(filestem), _t(t) {}
  ~WriteCheckPointMapping() {}

  void beginMoleculeIteration();
  void endMoleculeIteration();
  void handleMolecule(Molecule &molecule);

private:
  const simplemd::services::ParallelTopologyService &_parallelTopologyService;
  std::ofstream *_file;
  const std::string _filestem;
  const unsigned int _t;
  unsigned int _particleCounter;
  std::stringstream _moleculedata;
};

#endif // _MOLECULARDYNAMICS_MOLECULEMAPPINGS_WRITECHECKPOINTMAPPING_H_
