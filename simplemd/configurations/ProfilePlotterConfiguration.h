// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_CONFIGURATIONS_PROFILEPLOTTERCONFIGURATION_H_
#define _MOLECULARDYNAMICS_CONFIGURATIONS_PROFILEPLOTTERCONFIGURATION_H_

#include "simplemd/MolecularDynamicsUserInput.h"
#include "tarch/la/Vector.h"
#include "tarch/configuration/Configuration.h"

namespace simplemd {
  namespace configurations {
    class ProfilePlotterConfiguration;
  }
}


/** reads information for profile measurements. Each profile is evaluated on a sub-set of the linked cell data structure.
 *  @author Philipp Neumann
 */
class simplemd::configurations::ProfilePlotterConfiguration: public tarch::configuration::Configuration {
  public:
    ProfilePlotterConfiguration();
    virtual ~ProfilePlotterConfiguration(){}

    void parseSubtag( tinyxml2::XMLElement* node );

    /**
     * Return name of xml tag that is associated to the configuration.
     */
    std::string getTag() const {return "profile-plotter";}

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
    bool isValid() const { return _isValid;}

    /** getters for all parsed and computed quantities */
    const tarch::la::Vector<MD_DIM,unsigned int>& getStartCell() const {return _startCell;}
    const tarch::la::Vector<MD_DIM,unsigned int>& getRange() const {return _range;}
    const unsigned int& getWriteEveryTimestep() const { return _writeEveryTimestep;}
    const unsigned int& getSampleEveryTimestep() const { return _sampleEveryTimestep;}
    const unsigned int& getStartAtTimestep() const {return _startAtTimestep;}

  private:
    static const std::string START_CELL;
    static const std::string RANGE;
    static const std::string WRITE_EVERY_TIMESTEP;
    static const std::string SAMPLE_EVERY_TIMESTEP;
    static const std::string START_AT_TIMESTEP;

    unsigned int _writeEveryTimestep;
    unsigned int _sampleEveryTimestep;
    unsigned int _startAtTimestep;
    tarch::la::Vector<MD_DIM,unsigned int> _startCell;
    tarch::la::Vector<MD_DIM,unsigned int> _range;

    /** isValid flag */
    bool _isValid;
};
#endif // _MOLECULARDYNAMICS_CONFIGURATIONS_PROFILEPLOTTERCONFIGURATION_H_
