// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#ifndef _MOLECULARDYNAMICS_MOLECULEMAPPINGS_Adios2Writer_H_
#define _MOLECULARDYNAMICS_MOLECULEMAPPINGS_Adios2Writer_H_

#include "simplemd/MolecularDynamicsDefinitions.h"
#include "simplemd/Molecule.h"
#include "simplemd/configurations/MolecularDynamicsConfiguration.h"
#include "simplemd/services/MoleculeService.h"
#include "simplemd/services/ParallelTopologyService.h"
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include <adios2.h>
#include <chrono>
#include <errno.h>
#include <iomanip>
#include <map>
#include <memory>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <variant>
#include <vector>

#if (MD_PARALLEL == MD_YES)
#include <mpi.h>
#endif

namespace simplemd {
  namespace moleculemappings {
    class Adios2Writer;
  }
} //namespace simplemd


/** writes molecule data to a adios2 output directory
 *  @author Vincent Fürst
 */
class simplemd::moleculemappings::Adios2Writer {
  public:
    /** Constructor 
   * @param parallelTopologyService parallel topology service
   * @param moleculeService molecule service
   * @param filename name for Adios2 output directory
   * @param offset local domain offset
   * @param configuration molecular dynamics configuration
   * @param communicator MPI communicator in parallel case
   */
    Adios2Writer(const simplemd::services::ParallelTopologyService& parallelTopologyService,
      const simplemd::services::MoleculeService& moleculeService, const std::string& filename,
      tarch::la::Vector<MD_DIM,double> offset,
      const simplemd::configurations::MolecularDynamicsConfiguration &configuration
      #if (MD_PARALLEL==MD_YES)
      ,MPI_Comm communicator = MPI_COMM_WORLD
      #endif
      );

    /** Destructor */
    ~Adios2Writer();

    /** updates local variable _timestep to global timestep
     * @param timestep global timestep
    */
    void setTimestep(const unsigned int &timestep);

    /** begins new timestep in Adios2 IO before iterating over molecules */
    void beginMoleculeIteration();

    /** writes assembled data to file after iterating over all molecules */
    void endMoleculeIteration();

    /** writes position and velocity of the given molecule to the output vector 
     * @param molecule regarded molecule
    */
    void handleMolecule(Molecule &molecule);

  private:
    const simplemd::services::ParallelTopologyService& _parallelTopologyService;
    const simplemd::services::MoleculeService& _moleculeService;
    /** filename */
    std::string _filename;
    /** current timestep */
    unsigned int _timestep;
    /** domain offset*/
    tarch::la::Vector<MD_DIM,double> _offset;
    /** molecular dynamics simulation*/
    const simplemd::configurations::MolecularDynamicsConfiguration &_configuration;

    #if(MD_PARALLEL==MD_YES)
      /** mpi-communicator*/
      MPI_Comm _communicator;
    #endif

    /** number of MPI-threads */
    int _numberProcesses;

    /** center of global domain */
    std::array<double, 3> _domainCenter;

    /** file stream */
    std::ofstream _file;

    //** adios 2 variables storing positions, velocities, domain size and position and current timestep */
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

    /** Does all initializations needed for Adios2 output and writes time-consistent data */
    void initAdios2();

    /** Adios2 instance */
    std::shared_ptr<adios2::ADIOS> _inst;
    /** Adios2 engine */
	  std::shared_ptr<adios2::Engine> _engine;
    /** Adios2 IO */
	  std::shared_ptr<adios2::IO> _io;
};

#endif // _MOLECULARDYNAMICS_MOLECULEMAPPINGS_Adios2Writer_H_