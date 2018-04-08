// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_CONFIGURATIONS_MACROSCOPICCELLCONFIGURATION_H_
#define _MOLECULARDYNAMICS_COUPLING_CONFIGURATIONS_MACROSCOPICCELLCONFIGURATION_H_

#include "tarch/configuration/Configuration.h"
#include "tarch/la/Vector.h"
#include <iostream>

namespace coupling {
  namespace configurations {
    template<unsigned int dim>
    class MacroscopicCellConfiguration;
  }
}


/** configuration for output of macroscopic cell data to vtk files.
 *  @author Philipp Neumann
 */
template<unsigned int dim>
class coupling::configurations::MacroscopicCellConfiguration: public tarch::configuration::Configuration {
  public:
    MacroscopicCellConfiguration():
      _isValid(true),_macroscopicCellSize(0.0), _linkedCellsPerMacroscopicCell(0),
      _writeEveryMicroscopicTimestep(0), _microscopicFilename(""),
      _writeEveryMacroscopicTimestep(0), _macroscopicFilename(""){}

    ~MacroscopicCellConfiguration(){}

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
    bool isValid() const { return _isValid;}

    /** getters */
    tarch::la::Vector<dim,double> getMacroscopicCellSize() const { return _macroscopicCellSize;}

    tarch::la::Vector<dim,unsigned int> getNumberLinkedCellsPerMacroscopicCell() const {return _linkedCellsPerMacroscopicCell;}

    unsigned int getWriteEveryMicroscopicTimestep() const { return _writeEveryMicroscopicTimestep; }
    std::string getMicroscopicFilename() const { return _microscopicFilename;}

    unsigned int getWriteEveryMacroscopicTimestep() const { return _writeEveryMacroscopicTimestep; }
    std::string getMacroscopicFilename() const { return _macroscopicFilename;}

  private:
    static const std::string MACROSCOPIC_CELL_SIZE;
    static const std::string LINKED_CELLS_PER_MACROSCOPIC_CELL;
    static const std::string WRITE_EVERY_MICROSCOPIC_TIMESTEP;
    static const std::string WRITE_EVERY_MACROSCOPIC_TIMESTEP;
    static const std::string MICROSCOPIC_FILENAME;
    static const std::string MACROSCOPIC_FILENAME;

    bool _isValid;
    tarch::la::Vector<dim,double> _macroscopicCellSize;
    tarch::la::Vector<dim,unsigned int> _linkedCellsPerMacroscopicCell;
    unsigned int _writeEveryMicroscopicTimestep;
    std::string _microscopicFilename;
    unsigned int _writeEveryMacroscopicTimestep;
    std::string _macroscopicFilename;
};
#include "coupling/configurations/MacroscopicCellConfiguration.cpph"

#endif // _MOLECULARDYNAMICS_COUPLING_CONFIGURATIONS_MACROSCOPICCELLCONFIGURATION_H_
