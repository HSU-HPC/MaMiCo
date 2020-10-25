// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_CONFIGURATIONS_DOMAINCONFIGURATION_H_
#define _MOLECULARDYNAMICS_CONFIGURATIONS_DOMAINCONFIGURATION_H_

#include "simplemd/MolecularDynamicsDefinitions.h"
#include "tarch/configuration/Configuration.h"
#include "tarch/la/Vector.h"
#include <iostream>
#include <fstream>

namespace simplemd {
  namespace configurations {
    class DomainConfiguration;
  }
}


/** configuration input for domain description, including domain offset, domain size, boundary types etc.
 *  @author Philipp Neumann
 */
class simplemd::configurations::DomainConfiguration: public tarch::configuration::Configuration {
  public:
    DomainConfiguration();
    virtual ~DomainConfiguration(){}

    void parseSubtag( tinyxml2::XMLElement* node );

    /**
     * Return name of xml tag that is associated to the configuration.
     */
    std::string getTag() const;

    /**
     * Is config valid?
     *
     * This operation usually fails, if
     *
     * - parseSubtag() hasn't been called, i.e. configuration has not been
     *   used, or
     * - parseSubtag() failed due to a wrong file.
     *
     * If a tag ain't optional and parseSubtag() was not called (first case)
     */
    bool isValid() const;

    /** getters for all parsed and computed quantities */
    const tarch::la::Vector<MD_DIM,unsigned int>& getMoleculesPerDirection() const { return _moleculesPerDirection;}
    const tarch::la::Vector<MD_DIM,double>& getGlobalDomainSize() const { return _domainSize; }
    const tarch::la::Vector<MD_DIM,double>& getGlobalDomainOffset() const {return _domainOffset;}
    const double& getCutoffRadius() const { return _cutoffRadius;}
    const tarch::la::Vector<MD_DIM,double>& getMeshWidth() const { return _meshWidth; }
    const double& getKB() const { return _kB;}
    const unsigned int& getBlockSize() const { return _blockSize;}
    const tarch::la::Vector<MD_LINKED_CELL_NEIGHBOURS,simplemd::BoundaryType>& getBoundary() const {
      return _boundary;
    }
    const std::string& getCheckpointFilestem() const { return _checkpointFilestem;}
    const bool& initFromCheckpoint() const { return _initFromCheckpoint;}
    const bool& initFromSequentialCheckpoint() const { return _initFromSequentialCheckpoint;}

    void setCheckpointFilestem(const std::string & filestem) {
      _checkpointFilestem = filestem;
    }
    void setInitFromCheckpoint(const bool & initFromCheckpoint) {
      _initFromCheckpoint = initFromCheckpoint;
    }
    void setInitFromSequentialCheckpoint (const bool & initFromSequentialCheckpoint) {
      _initFromSequentialCheckpoint = initFromSequentialCheckpoint;
    }

    unsigned int getNumberOfMolecules() const;

  private:
    static const std::string MOLECULES_PER_DIRECTION;
    static const std::string DOMAIN_SIZE;
    static const std::string DOMAIN_OFFSET;
    static const std::string CUTOFF_RADIUS;
    static const std::string LINKED_CELL_SIZE;
    static const std::string K_B;
    static const std::string BLOCK_SIZE;
    static const std::string BOUNDARY[MD_LINKED_CELL_NEIGHBOURS];
    static const std::string PERIODIC_BOUNDARY;
    static const std::string GEOMETRY_BOUNDARY;
    static const std::string OPEN_BOUNDARY;
    static const std::string REFLECTING_BOUNDARY;
    static const std::string RDF_FILENAME;
    static const std::string CELLS_PER_LINKED_CELL;
    static const std::string INIT_FROM_CHECKPOINT;
    static const std::string INIT_FROM_SEQUENTIAL_CHECKPOINT;
    static const std::string LINKED_CELLS_PER_NUMBER_DENSITY_EVALUATION;

    /** number of molecules in each direction */
    tarch::la::Vector<MD_DIM,unsigned int> _moleculesPerDirection;

    /** global domain size */
    tarch::la::Vector<MD_DIM,double> _domainSize;

    /** global domain offset */
    tarch::la::Vector<MD_DIM,double> _domainOffset;

    /** cut off radius for lennard jones potential. Determines the meshsize in the simulation. */
    double _cutoffRadius;

    /** size of linked cells. If not handed over in the config, this parameter will be determined automatically */
    tarch::la::Vector<MD_DIM,double> _meshWidth;

    /** dimensionless Boltzmann's constant */
    double _kB;

    /** blocksize to be used for molecule storage */
    unsigned int _blockSize;

    /** boundary types for all outer boundaries */
    tarch::la::Vector<MD_LINKED_CELL_NEIGHBOURS,simplemd::BoundaryType> _boundary;

    std::string _checkpointFilestem;
    bool _initFromCheckpoint;
    bool _initFromSequentialCheckpoint;

    /** isValid flag */
    bool _isValid;
};
#endif // _MOLECULARDYNAMICS_CONFIGURATIONS_DOMAINCONFIGURATION_H_
