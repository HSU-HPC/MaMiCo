// Copyright (C) 2016 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _EVAPORATIONCONFIG_H_
#define _EVAPORATIONCONFIG_H_

#include "tarch/la/Vector.h"

namespace coupling {
namespace configurations {
struct EvaporationConfig;
}
}
/** Configuration parameters for Evaporation flow scenario
 *	@brief Configuration parameters for Evaporation flow scenario
 *	@author Helene Wittenberg
 */
struct coupling::configurations::EvaporationConfig {

public:
  /** Defines the type of md solver for the coupled simulation  */
  /*enum MicroSolverType {
    SIMPLEMD = 0, ///< the SimpleMD solver is used
    SYNTHETIC = 1 ///< the synthetic solver is used
  };*/

  /** @brief creates EvaporationConfig if all elements exist and can be read
   * 	@param filename file name of the configuation xml file
   */
  static EvaporationConfig parseConfiguration(const std::string& filename) {
    EvaporationConfig _EvapConfig;

    tinyxml2::XMLDocument conffile;
    tinyxml2::XMLElement* node = NULL;
    conffile.LoadFile("evaporation.xml");
    node = conffile.FirstChildElement("evaporation-test");
    if (node == NULL) {
      std::cout << "Could not read input file evaporation.xml: missing element <evaporation-test>" << std::endl;
      exit(EXIT_FAILURE);
    }
    tinyxml2::XMLElement* n_mamico = node->NextSiblingElement();
    if (n_mamico == NULL) {
      std::cout << "Could not read input file evaporation.xml: missing element <mamico>" << std::endl;
      exit(EXIT_FAILURE);
    }
    tinyxml2::XMLElement* n_md = n_mamico->NextSiblingElement();
    if (n_md == NULL) {
      std::cout << "Could not read input file evaporation.xml: missing element <molecular-dynamics>" << std::endl;
      exit(EXIT_FAILURE);
    }
    tinyxml2::XMLElement* n_fp = n_md->NextSiblingElement();
    if (n_fp == NULL) {
      std::cout << "Could not read input file evaporation.xml: missing element <filter-pipeline>" << std::endl;
      exit(EXIT_FAILURE);
    }
    tinyxml2::XMLElement* n_unexpected = n_fp->NextSiblingElement();
    if (n_unexpected != NULL) {
      std::cout << "Could not read input file evaporation.xml: unknown element " << n_unexpected->Name() << std::endl;
      exit(EXIT_FAILURE);
    }

    tinyxml2::XMLElement* subtag = node->FirstChildElement("coupling");
    if (subtag == NULL) {
      std::cout << "Could not read input file evaporation.xml: Missing subtag: coupling" << std::endl;
      exit(EXIT_FAILURE);
    }
    tarch::configuration::ParseConfiguration::readIntMandatory(_EvapConfig.couplingCycles, subtag, "coupling-cycles");
    tarch::configuration::ParseConfiguration::readIntMandatory(_EvapConfig.csvEveryTimestep, subtag, "write-csv-every-timestep");
    subtag = node->FirstChildElement("microscopic-solver");
    if (subtag == NULL) {
      std::cout << "Could not read input file evaporation.xml: Missing subtag: microscopic-solver" << std::endl;
      exit(EXIT_FAILURE);
    }
    std::string type;
    tarch::configuration::ParseConfiguration::readStringMandatory(type, subtag, "type");
    if (type == "md") {
      tarch::configuration::ParseConfiguration::readDoubleMandatory(_EvapConfig.temp, subtag, "temperature");
      tarch::configuration::ParseConfiguration::readIntMandatory(_EvapConfig.totalNumberMDSimulations, subtag, "number-md-simulations");
      if (_EvapConfig.totalNumberMDSimulations < 1) {
        std::cout << "Could not read input file evaporation.xml: number-md-simulations < 1" << std::endl;
        exit(EXIT_FAILURE);
      }
    } else {
      std::cout << "Could not read input file evaporation.xml: Unknown microscopic solver type!" << std::endl;
      exit(EXIT_FAILURE);
    }
    subtag = node->FirstChildElement("macroscopic-solver");
    if (subtag == NULL) {
      std::cout << "Could not read input file evaporation.xml: Missing subtag: macroscopic-solver" << std::endl;
      exit(EXIT_FAILURE);
    }

    return _EvapConfig;
  }

#if (BUILD_WITH_OPENFOAM)
  /** @brief holds the variables necessary for the interface to the IcoFoam
   * solver  */
  struct FoamConfig {
    /** @brief the path to the directory, where the folder with the IcoFoam
     * OpenFOAM setup files are */
    std::string directory;
    /** @brief the name of the OpenFOAM folder */
    std::string folder;
    /** @brief the cubic mesh with a cubic hole has 12 boundaries, this vector
     * tells which of them
     *         shall be coupled with the md; The order has to be the same as in
     * the blockMeshDict */
    tarch::la::Vector<12, unsigned int> boundariesWithMD;
  };
#endif

/** @brief channel is always expected to have origin at (0.0,0.0,0.0) and to have the same width as height */
double channelheight;
/** @brief the total length of the channel */
double channellength;
/** @brief number of coupling cycles, that is continuum time steps; MD/DPD: 1000 */
int couplingCycles;
/** @brief the time step interval for writing the md data to csv files  */
int csvEveryTimestep;
/** @brief the general density of the fluid under consideration */
double density;
/** @brief the start temperature for the fluid under consideration */
double temp;
/** @brief the number of md simulation instances in a multi-instance coupling  */
int totalNumberMDSimulations;

#if (BUILD_WITH_OPENFOAM)
  /** @brief the configurations for the OpenFoam solver */
  FoamConfig foam;
#endif
};

#endif
