// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#if BUILD_WITH_ADIOS2

#ifndef _MOLECULARDYNAMICS_MOLECULEMAPPINGS_Adios2Writer_H_
#define _MOLECULARDYNAMICS_MOLECULEMAPPINGS_Adios2Writer_H_

#include <string>
#include "simplemd/Molecule.h"
#include "simplemd/services/MoleculeService.h"
#include "simplemd/MolecularDynamicsDefinitions.h"
#include "simplemd/services/ParallelTopologyService.h"
#include "simplemd/configurations/MolecularDynamicsConfiguration.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <cstdlib>

#include <errno.h>
#include <chrono>
#include <iomanip>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <variant>
#include <vector>
#include <numeric>
#include <adios2.h>


#if (MD_PARALLEL==MD_YES)
#include <mpi.h>
#endif


namespace simplemd {
  namespace moleculemappings {
    class Adios2Writer;
  }
}


/** writes molecule data to a vtk file.
 *  In case of parallel computations, the rank of the respective process is added to the filename.
 *  @author Philipp Neumann
 */
class simplemd::moleculemappings::Adios2Writer {
  public:
    Adios2Writer(const simplemd::services::ParallelTopologyService& parallelTopologyService,
      const simplemd::services::MoleculeService& moleculeService, const std::string& filename,
      tarch::la::Vector<MD_DIM,double> offset,
      int numberProcesses,
      const simplemd::configurations::MolecularDynamicsConfiguration &configuration
      #if (MD_PARALLEL==MD_YES)
      ,MPI_Comm communicator = MPI_COMM_WORLD
      #endif
      );
    ~Adios2Writer();
    void setTimestep(const unsigned int &timestep);

    void beginMoleculeIteration();
    void endMoleculeIteration();
    void handleMolecule(Molecule &molecule);

  private:
    const simplemd::services::ParallelTopologyService& _parallelTopologyService;
    const simplemd::services::MoleculeService& _moleculeService;
    /** filename */
    std::string _filename;
    /** current timestep */
    unsigned int _timestep;

    tarch::la::Vector<MD_DIM,double> _offset;
    int _numberProcesses;
    const simplemd::configurations::MolecularDynamicsConfiguration &_configuration;

    /** file stream */
    std::ofstream _file;

    adios2::Variable<float> rx_var;
    adios2::Variable<float> ry_var;
    adios2::Variable<float> rz_var;
    adios2::Variable<float> vx_var;
    adios2::Variable<float> vy_var;
    adios2::Variable<float> vz_var;
    adios2::Variable<float> global_box_var;
    adios2::Variable<uint64_t> component_id_var;
    adios2::Variable<double> time_var;

    /** stores the positions */
    std::vector<float> _positionsx;
    std::vector<float> _positionsy;
    std::vector<float> _positionsz;

    std::vector<uint64_t> _component_id;
    
    /** stores the velocities */
    std::vector<float> _velocitiesx;
    std::vector<float> _velocitiesy;
    std::vector<float> _velocitiesz;

    /** stores the information whether a particle is fixed in space */
    std::stringstream _fix;

    #if (MD_PARALLEL==MD_YES)
    void initAdios2(MPI_Comm communicator);
    #else
    void initAdios2();
    #endif
    std::shared_ptr<adios2::ADIOS> _inst;
	  std::shared_ptr<adios2::Engine> _engine;
	  std::shared_ptr<adios2::IO> _io;

    std::array<double, 3> _domainCenter;
};

#endif // _MOLECULARDYNAMICS_MOLECULEMAPPINGS_Adios2Writer_H_

#endif