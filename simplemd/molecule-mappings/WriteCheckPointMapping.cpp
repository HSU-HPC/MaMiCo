// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#include "simplemd/molecule-mappings/WriteCheckPointMapping.h"

#include "simplemd/MolecularDynamicsDefinitions.h"

void
simplemd::moleculemappings::WriteCheckPointMapping::beginMoleculeIteration() {
  std::stringstream ss;
  ss << _filestem << "_" << _t << "_";
#if (MD_PARALLEL == MD_YES)
  ss << _parallelTopologyService.getRank();
#else
  ss << "0";
#endif
  ss << ".checkpoint";
  if (_file != NULL) {
    delete _file;
    _file = NULL;
  }
  _file = new std::ofstream(ss.str().c_str());
  if (!_file->is_open()) {
    std::cout << "ERROR WriteCheckPointMapping: Could not open file "
              << ss.str() << "!" << std::endl;
    exit(EXIT_FAILURE);
  }
  _file->precision(30);
  _moleculedata.precision(30);

  // reset particle counter
  _particleCounter = 0;
}

void
simplemd::moleculemappings::WriteCheckPointMapping::endMoleculeIteration() {
  *_file << _particleCounter << " " << MD_DIM << std::endl;
  *_file << _moleculedata.str() << std::endl;
  _moleculedata.clear();
  _moleculedata.flush();
  _file->close();
  delete _file;
  _file = NULL;
}

void simplemd::moleculemappings::WriteCheckPointMapping::handleMolecule(
    Molecule &molecule) {
  const tarch::la::Vector<MD_DIM, double> &forceOld =
      molecule.getConstForceOld();
  const tarch::la::Vector<MD_DIM, double> &vel = molecule.getConstVelocity();
  const tarch::la::Vector<MD_DIM, double> &pos = molecule.getConstPosition();
  for (unsigned int d = 0; d < MD_DIM; d++) {
    _moleculedata << std::fixed << pos[d] << " ";
  }
  for (unsigned int d = 0; d < MD_DIM; d++) {
    _moleculedata << std::fixed << vel[d] << " ";
  }
  for (unsigned int d = 0; d < MD_DIM - 1; d++) {
    _moleculedata << std::fixed << forceOld[d] << " ";
  }
  _moleculedata << std::fixed << forceOld[MD_DIM - 1] << std::endl;
  _particleCounter++;
}
