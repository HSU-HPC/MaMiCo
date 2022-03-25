// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_MOLECULEMAPPINGS_VELOCITYSTOERMERVERLETMAPPING_H_
#define _MOLECULARDYNAMICS_MOLECULEMAPPINGS_VELOCITYSTOERMERVERLETMAPPING_H_

#include "simplemd/MolecularDynamicsDefinitions.h"
#include "simplemd/Molecule.h"
#include "simplemd/services/MolecularPropertiesService.h"

namespace simplemd {
namespace moleculemappings {
class VelocityStoermerVerletMapping;
}
} // namespace simplemd

/** carries out the time integration for position and velocity of molecule.
 *  @author Philipp Neumann
 */
class simplemd::moleculemappings::VelocityStoermerVerletMapping {
public:
  VelocityStoermerVerletMapping(const double &kB, const double &dt, const double &mass,
                                const tarch::la::Vector<MD_LINKED_CELL_NEIGHBOURS, simplemd::BoundaryType> &boundary,
                                const tarch::la::Vector<MD_DIM, double> &domainOffset, const tarch::la::Vector<MD_DIM, double> &domainSize);
  ~VelocityStoermerVerletMapping();

  void beginMoleculeIteration();
  void endMoleculeIteration();

  void handleMolecule(Molecule &molecule);

private:
  tarch::la::Vector<2 * MD_DIM, bool> initReflectingBoundary(const tarch::la::Vector<MD_LINKED_CELL_NEIGHBOURS, simplemd::BoundaryType> &boundary) const;

  /** timestep */
  const double _dt;

  /** help parameter */
  const double _a;

  const tarch::la::Vector<MD_DIM, double> _zero;

  const tarch::la::Vector<2 * MD_DIM, bool> _boundary;
  const tarch::la::Vector<MD_DIM, double> _domainOffset;
  const tarch::la::Vector<MD_DIM, double> _domainSize;
};

#endif // _MOLECULARDYNAMICS_MOLECULEMAPPINGS_VELOCITYSTOERMERVERLETMAPPING_H_
