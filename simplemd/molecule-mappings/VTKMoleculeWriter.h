// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_MOLECULEMAPPINGS_VTKMOLECULEWRITER_H_
#define _MOLECULARDYNAMICS_MOLECULEMAPPINGS_VTKMOLECULEWRITER_H_

#include "simplemd/MolecularDynamicsDefinitions.h"
#include "simplemd/Molecule.h"
#include "simplemd/MoleculeContainer.h"
#include "simplemd/services/ParallelTopologyService.h"
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

namespace simplemd {
namespace moleculemappings {
class VTKMoleculeWriter;
}
} // namespace simplemd

/** writes molecule data to a vtk file.
 *  In case of parallel computations, the rank of the respective process is
 * added to the filename.
 *  @author Philipp Neumann
 */
class simplemd::moleculemappings::VTKMoleculeWriter {
public:
  /** Constructor
   * @param parallelTopologyService parallel topology service
   * @param moleculeContainer molecule container
   * @param filename filename for .vtk output file
   */
  VTKMoleculeWriter(const simplemd::services::ParallelTopologyService& parallelTopologyService, const simplemd::MoleculeContainer& moleculeContainer,
                    const std::string& filename);

  /** Destructor */
  ~VTKMoleculeWriter();

  /** updates local variable _timestep to global timestep
   * @param timestep global timestep
   */
  void setTimestep(const unsigned int& timestep);

  /** Opens new .vtk file for current timestep, writes meta-data, empty buffers */
  void beginMoleculeIteration();

  /** writes assembled data to file after iterating over all molecules */
  void endMoleculeIteration();

  /** adds position and velocity of the given molecule to the temporary stringstreams
   * @param molecule regarded molecule
   */
  void handleMolecule(Molecule& molecule);

  static const bool IsParallel = false;

private:
  const simplemd::services::ParallelTopologyService& _parallelTopologyService;
  const simplemd::MoleculeContainer& _moleculeContainer;
  /** filename */
  std::string _filename;
  /** current timestep */
  unsigned int _timestep;

  /** file stream */
  std::ofstream _file;

  /** stores the positions */
  std::stringstream _positions;

  /** stores the velocities */
  std::stringstream _velocities;

  /** stores the (old) forces */
  std::stringstream _forces;

  /** stores the information whether a particle is fixed in space */
  std::stringstream _fix;

  // stores the particleID
  std::stringstream _id;
};

#endif // _MOLECULARDYNAMICS_MOLECULEMAPPINGS_VTKMOLECULEWRITER_H_
