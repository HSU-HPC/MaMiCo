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
template <unsigned int dim> class CouplingCellConfiguration;
}
} // namespace coupling

/** configuration for output of macroscopic cell data to vtk files. Derive from
 *the class tarch::configuration::Configuration
 *	@brief configuration for output of macroscopic cell data to vtk files.
 *	@tparam dim Number of dimensions; it can be 1, 2 or 3
 *  @author Philipp Neumann
 */
template <unsigned int dim> class coupling::configurations::CouplingCellConfiguration : public tarch::configuration::Configuration {
public:
  /** Constructor, initializes the class  */
  CouplingCellConfiguration()
      : _isValid(true), _macroscopicCellSize(0.0), _linkedCellsPerMacroscopicCell(0), _writeEveryMicroscopicTimestep(0), _microscopicFilename(""),
        _writeEveryMacroscopicTimestep(0), _macroscopicFilename("") {}

  /** Destructor */
  ~CouplingCellConfiguration() {}

  /** parseSubtag
   * 	@param node
   */
  void parseSubtag(tinyxml2::XMLElement* node);

  /** Returns name of xml tag that is associated to the configuration.
   * 	@return name of xml tag that is associated to the configuration
   */
  std::string getTag() const;

  /** checks if the configuration is valid. This operation usually fails, if
   *e.g.
   *	1. parseSubtag() hasn't been called, i.e. configuration has not been
   *used, or
   *  2. parseSubtag() failed due to a wrong file.
   * 	@return _isValid
   */
  bool isValid() const { return _isValid; }

  /** Returns the macroscopic cell size
   * 	@return _macroscopicCellSize
   */
  tarch::la::Vector<dim, double> getCouplingCellSize() const { return _macroscopicCellSize; }

  /** Returns the number of linked cell encapsulated within a macroscopic cell
   * 	@return _linkedCellsPerMacroscopicCell
   */
  tarch::la::Vector<dim, unsigned int> getNumberLinkedCellsPerCouplingCell() const { return _linkedCellsPerMacroscopicCell; }

  /** Returns step interval, at which the micro infos is to be written out
   * 	@return _writeEveryMicroscopicTimestep
   */
  unsigned int getWriteEveryMicroscopicTimestep() const { return _writeEveryMicroscopicTimestep; }
  /** Returns the microscopic file name
   * 	@return _microscopicFilename
   */
  std::string getMicroscopicFilename() const { return _microscopicFilename; }

  /** Returns step interval, at which the macro infos is to be written out
   * 	@return _writeEveryMacroscopicTimestep
   */
  unsigned int getWriteEveryMacroscopicTimestep() const { return _writeEveryMacroscopicTimestep; }
  /** Returns the macroscopic file name
   * 	@return _maroscopicFilename
   */
  std::string getMacroscopicFilename() const { return _macroscopicFilename; }

protected:
  CouplingCellConfiguration(tarch::la::Vector<dim, double> macroscopicCellSize, tarch::la::Vector<dim, unsigned int> linkedCellsPerMacroscopicCell,
                            unsigned int writeEveryMicroscopicTimestep, std::string microscopicFilename, unsigned int writeEveryMacroscopicTimestep,
                            std::string macroscopicFilename)
      : _isValid(true), _macroscopicCellSize(macroscopicCellSize), _linkedCellsPerMacroscopicCell(linkedCellsPerMacroscopicCell),
        _writeEveryMicroscopicTimestep(writeEveryMicroscopicTimestep), _microscopicFilename(microscopicFilename),
        _writeEveryMacroscopicTimestep(writeEveryMacroscopicTimestep), _macroscopicFilename(macroscopicFilename) {}

private:
  static const std::string MACROSCOPIC_CELL_SIZE;
  static const std::string LINKED_CELLS_PER_MACROSCOPIC_CELL;
  static const std::string WRITE_EVERY_MICROSCOPIC_TIMESTEP;
  static const std::string WRITE_EVERY_MACROSCOPIC_TIMESTEP;
  static const std::string MICROSCOPIC_FILENAME;
  static const std::string MACROSCOPIC_FILENAME;

  bool _isValid;
  tarch::la::Vector<dim, double> _macroscopicCellSize;
  tarch::la::Vector<dim, unsigned int> _linkedCellsPerMacroscopicCell;
  unsigned int _writeEveryMicroscopicTimestep;
  std::string _microscopicFilename;
  unsigned int _writeEveryMacroscopicTimestep;
  std::string _macroscopicFilename;
};
#include "coupling/configurations/MacroscopicCellConfiguration.cpph"

#endif // _MOLECULARDYNAMICS_COUPLING_CONFIGURATIONS_MACROSCOPICCELLCONFIGURATION_H_
