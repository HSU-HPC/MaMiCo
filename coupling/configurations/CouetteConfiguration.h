// Copyright (C) 2016 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _COUETTECONFIG_H_
#define _COUETTECONFIG_H_

#include "tarch/configuration/ParseConfiguration.h"
#include "tarch/la/Vector.h"

using tarch::configuration::ParseConfiguration;

namespace coupling {
namespace configurations {

struct CouetteConfig;
}
} // namespace coupling
/** Configuration parameters for Couette flow scenario
 *	@brief Configuration parameters for Couette flow scenario
 *
 *
 *	@author Niklas Wittmer
 */
struct coupling::configurations::CouetteConfig {

private:
  // CouetteConfig(){}

public:
  /** Defines the type of continuum solver for the coupled simulation  */
  enum MacroSolverType {
    COUETTE_ANALYTICAL = 0, ///< the analytical couette solver is used
                            ///< (coupling::solvers::CouetteSolver)
    COUETTE_LB = 1,         ///< the Lattice-Bolzmann solver is used
                            ///< (coupling::solvers::LBCouetteSolver)
    COUETTE_FD = 2,         ///< the 1d finite-difference solver is used
                            ///< (coupling::solvers::FiniteDifferenceSolver)
    COUETTE_FOAM = 3        ///< the IcoFoam solver is used (coupling::solvers::IcoFoam)
  };
  /** Defines the type of md solver for the coupled simulation  */
  enum MicroSolverType {
    SIMPLEMD = 0,  ///< the SimpleMD solver is used
    SYNTHETIC = 1, ///< the synthetic solver is used
    LS1 = 2        ///< the LS1 solver is used
  };

  /** @brief creates CouetteConfig if all elements exist and can be read
   * 	@param filename
   */
  static CouetteConfig parseCouetteConfiguration(const std::string& filename) {
    CouetteConfig cfg;

    ParseConfiguration::XMLConfiguration xmlConfig = ParseConfiguration::XMLConfiguration::load(filename);
    tinyxml2::XMLElement* node = xmlConfig.root->FirstChildElement("couette-test");

    if (node == NULL) {
      std::cout << "Could not read input file " << filename
                << ": missing element "
                   "<couette-test>"
                << std::endl;
      exit(EXIT_FAILURE);
    }

    tinyxml2::XMLElement* n_mamico = node->NextSiblingElement();
    if (n_mamico == NULL) {
      std::cout << "Could not read input file " << filename << ": missing element <mamico>" << std::endl;
      exit(EXIT_FAILURE);
    }
    tinyxml2::XMLElement* n_md = n_mamico->NextSiblingElement();
    if (n_md == NULL) {
      std::cout << "Could not read input file " << filename
                << ": missing element "
                   "<molecular-dynamics>"
                << std::endl;
      exit(EXIT_FAILURE);
    }
    tinyxml2::XMLElement* n_fp = n_md->NextSiblingElement();
    if (n_fp == NULL) {
      std::cout << "Could not read input file " << filename
                << ": missing element "
                   "<filter-pipeline>"
                << std::endl;
      exit(EXIT_FAILURE);
    }
    tinyxml2::XMLElement* n_unexpected = n_fp->NextSiblingElement();
    if (n_unexpected != NULL) {
      std::cout << "Could not read input file " << filename << ": unknown element " << n_unexpected->Name() << std::endl;
      exit(EXIT_FAILURE);
    }

    tinyxml2::XMLElement* subtag = node->FirstChildElement("domain");
    if (subtag == NULL) {
      std::cout << "Could not read input file " << filename << ": Missing subtag: domain" << std::endl;
      exit(EXIT_FAILURE);
    }
    tarch::configuration::ParseConfiguration::readDoubleMandatory(cfg.channelheight, subtag, "channelheight");
    tarch::configuration::ParseConfiguration::readVectorMandatory<3, double>(cfg.wallVelocity, subtag, "wall-velocity");
    cfg.wallInitCycles = 0;
    tarch::configuration::ParseConfiguration::readIntOptional(cfg.wallInitCycles, subtag, "wall-init-cycles");
    if (cfg.wallInitCycles > 0)
      tarch::configuration::ParseConfiguration::readVectorMandatory<3, double>(cfg.wallInitVelocity, subtag, "wall-init-velocity");
    cfg.wallOscillations = 0;
    tarch::configuration::ParseConfiguration::readDoubleOptional(cfg.wallOscillations, subtag, "wall-oscillations");

    subtag = node->FirstChildElement("coupling");
    if (subtag == NULL) {
      std::cout << "Could not read input file couette.xml: Missing subtag: coupling" << std::endl;
      exit(EXIT_FAILURE);
    }
    tarch::configuration::ParseConfiguration::readIntMandatory(cfg.couplingCycles, subtag, "coupling-cycles");
    tarch::configuration::ParseConfiguration::readBoolMandatory(cfg.twoWayCoupling, subtag, "two-way-coupling");
    tarch::configuration::ParseConfiguration::readBoolMandatory(cfg.md2Macro, subtag, "send-from-md-to-macro");
    tarch::configuration::ParseConfiguration::readBoolMandatory(cfg.macro2Md, subtag, "send-from-macro-to-md");
    tarch::configuration::ParseConfiguration::readIntMandatory(cfg.filterInitCycles, subtag, "filter-init-cycles");
    tarch::configuration::ParseConfiguration::readIntMandatory(cfg.csvEveryTimestep, subtag, "write-csv-every-timestep");
    tarch::configuration::ParseConfiguration::readBoolMandatory(cfg.computeSNR, subtag, "compute-snr");

    subtag = node->FirstChildElement("microscopic-solver");
    if (subtag == NULL) {
      std::cout << "Could not read input file " << filename
                << ": Missing subtag: "
                   "microscopic-solver"
                << std::endl;
      exit(EXIT_FAILURE);
    }
    std::string type;
    tarch::configuration::ParseConfiguration::readStringMandatory(type, subtag, "type");
    if (type == "md") {
      cfg.miSolverType = SIMPLEMD;
      tarch::configuration::ParseConfiguration::readDoubleMandatory(cfg.temp, subtag, "temperature");
      tarch::configuration::ParseConfiguration::readIntMandatory(cfg.equSteps, subtag, "equilibration-steps");

      std::string num_MD;
      tarch::configuration::ParseConfiguration::readStringMandatory(num_MD, subtag, "number-md-simulations");
      if (num_MD == "dynamic") {
        cfg.totalNumberMDSimulations = -1;
        cfg.lowerBoundNumberMDSimulations = 1;
        tarch::configuration::ParseConfiguration::readIntOptional(cfg.lowerBoundNumberMDSimulations, subtag, "min-number-md");
        tarch::configuration::ParseConfiguration::readDoubleMandatory(cfg.absVelErrStart, subtag, "error-start");
        tarch::configuration::ParseConfiguration::readDoubleMandatory(cfg.absVelErrEnd, subtag, "error-end");
      } else {
        tarch::configuration::ParseConfiguration::readIntMandatory(cfg.totalNumberMDSimulations, subtag, "number-md-simulations");
        if (cfg.totalNumberMDSimulations < 1) {
          std::cout << "Could not read input file " << filename
                    << ": "
                       "number-md-simulations < 1"
                    << std::endl;
          exit(EXIT_FAILURE);
        }
      }
    } else if (type == "synthetic") {
      cfg.miSolverType = SYNTHETIC;
      tarch::configuration::ParseConfiguration::readDoubleMandatory(cfg.noiseSigma, subtag, "noise-sigma");
      cfg.totalNumberMDSimulations = 1;
      tarch::configuration::ParseConfiguration::readIntOptional(cfg.totalNumberMDSimulations, subtag, "number-md-simulations");
    } else if (type == "ls1") {
      cfg.miSolverType = LS1;

      cfg.totalNumberMDSimulations = 1;
      tarch::configuration::ParseConfiguration::readDoubleMandatory(cfg.temp, subtag, "temperature");
      tarch::configuration::ParseConfiguration::readIntMandatory(cfg.equSteps, subtag, "equilibration-steps");
      tarch::configuration::ParseConfiguration::readIntOptional(cfg.totalNumberMDSimulations, subtag, "number-md-simulations");
    } else {
      std::cout << "Could not read input file " << filename
                << ": Unknown microscopic "
                   "solver type!"
                << std::endl;
      exit(EXIT_FAILURE);
    }
    tarch::configuration::ParseConfiguration::readDoubleMandatory(cfg.density, subtag, "density");

    subtag = node->FirstChildElement("macroscopic-solver");
    if (subtag == NULL) {
      std::cout << "Could not read input file " << filename
                << ": Missing subtag: "
                   "macroscopic-solver"
                << std::endl;
      exit(EXIT_FAILURE);
    }
    cfg.lbNumberProcesses = tarch::la::Vector<3, unsigned int>(1);
    tarch::configuration::ParseConfiguration::readStringMandatory(type, subtag, "type");
    if (type == "lb") {
      cfg.maSolverType = COUETTE_LB;
      tarch::configuration::ParseConfiguration::readVectorMandatory<3, unsigned int>(cfg.lbNumberProcesses, subtag, "number-of-processes");
      tarch::configuration::ParseConfiguration::readIntMandatory(cfg.plotEveryTimestep, subtag, "plot-every-timestep");
    } else if (type == "fd") {
      cfg.maSolverType = COUETTE_FD;
      tarch::configuration::ParseConfiguration::readVectorMandatory<3, unsigned int>(cfg.lbNumberProcesses, subtag, "number-of-processes");
      tarch::configuration::ParseConfiguration::readIntMandatory(cfg.plotEveryTimestep, subtag, "plot-every-timestep");
    }
#if (BUILD_WITH_OPENFOAM)
    else if (type == "foam") {
      cfg.maSolverType = COUETTE_FOAM;
      tarch::configuration::ParseConfiguration::readIntMandatory(cfg.plotEveryTimestep, subtag, "plot-every-timestep");
      tarch::configuration::ParseConfiguration::readStringMandatory(cfg.foam.directory, subtag, "foam-setup-directory");
      tarch::configuration::ParseConfiguration::readStringMandatory(cfg.foam.folder, subtag, "foam-setup-folder");
      tarch::configuration::ParseConfiguration::readVectorMandatory<12, unsigned int>(cfg.foam.boundariesWithMD, subtag, "boundaries-with-MD");
      if (!cfg.twoWayCoupling && cfg.foam.boundariesWithMD != tarch::la::Vector<12, unsigned int>{0}) {
        std::cout << "ERROR: Two-way coupling is disabled, but boundaries with "
                     "MD for openfoam were defined"
                  << std::endl;
        exit(EXIT_FAILURE);
      }
    }
#endif
    else if (type == "analytical") {
      cfg.maSolverType = COUETTE_ANALYTICAL;
      if (!(cfg.wallVelocity[1] == 0.0 && cfg.wallVelocity[2] == 0.0)) {
        std::cout << "analytic solver only supports flow in x-direction" << std::endl;
        exit(EXIT_FAILURE);
      }
    } else {
      std::cout << "Could not read input file " << filename
                << ": Unknown macroscopic "
                   "solver type!"
                << std::endl;
      exit(EXIT_FAILURE);
    }
    double vis;
    tarch::configuration::ParseConfiguration::readDoubleMandatory(vis, subtag, "viscosity");
    cfg.kinVisc = vis / cfg.density;
    tarch::configuration::ParseConfiguration::readIntMandatory(cfg.initAdvanceCycles, subtag, "init-advance-cycles");

    subtag = node->FirstChildElement("tws-loop");
    if (subtag == NULL)
      cfg.twsLoop = false;
    else {
      cfg.twsLoop = true;
      tarch::configuration::ParseConfiguration::readIntMandatory(cfg.twsLoopMin, subtag, "min");
      tarch::configuration::ParseConfiguration::readIntMandatory(cfg.twsLoopMax, subtag, "max");
      cfg.twsLoopStep = 1;
      tarch::configuration::ParseConfiguration::readIntOptional(cfg.twsLoopStep, subtag, "step");
    }

    if (cfg.miSolverType == SYNTHETIC) {
      if (cfg.macro2Md || cfg.totalNumberMDSimulations > 1 || cfg.lbNumberProcesses[0] != 1 || cfg.lbNumberProcesses[1] != 1 ||
          cfg.lbNumberProcesses[2] != 1) {
        std::cout << "Invalid configuration: Synthetic MD runs sequentially on "
                     "rank 0 only. "
                  << "It does neither support parallel communication nor "
                     "multi-instance sampling"
                  << std::endl;
        exit(EXIT_FAILURE);
      }
    }
    if (cfg.maSolverType == COUETTE_ANALYTICAL) {
      if (cfg.twoWayCoupling) {
        std::cout << "Invalid configuration: COUETTE_ANALYTICAL does not "
                     "support twoWayCoupling"
                  << std::endl;
        exit(EXIT_FAILURE);
      }
    }

    return cfg;
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

  /** @brief channel is always expected to have origin at (0.0,0.0,0.0) and to
   * be cubic (MD 30: 50.0, MD 60: 100.0, MD 120: 200.0) */
  double channelheight;
  /** @brief velocity of moving wall (lower boundary moves) */
  tarch::la::Vector<3, double> wallVelocity;
  /** @brief number of coupling cycles that will be done, while the wall
   * velocity is the wallInitVelocity */
  int wallInitCycles;
  /** @brief the velocity at the moving wall (z=0) while the wallInitCycles*/
  tarch::la::Vector<3, double> wallInitVelocity;
  /** @brief total number of oscillation periods of the oscillating wall
   * velocity*/
  double wallOscillations;
  /** @brief number of coupling cycles, that is continuum time steps; MD/DPD:
   * 1000 */
  int couplingCycles;
  /** @brief true if data will be send & received from the md to the continuum
   * solver */
  bool md2Macro;
  /** @brief true if data will be send & received from the continuum to the md
   * solver  */
  bool macro2Md;
  /** @brief true if the signal to noise ratio shall be evaluated */
  bool computeSNR;
  /** @brief true if from the md solver will be applied as boundary condition
   * for the continuum solver */
  bool twoWayCoupling;
  /** @brief number of coupling cycles before a filter is applied*/
  int filterInitCycles;
  /** @brief the time step interval for writing the md data to csv files  */
  int csvEveryTimestep;
  /** @brief the type of continuum solver */
  MacroSolverType maSolverType;
  /** @brief the type of md solver */
  MicroSolverType miSolverType;
  /** @brief the general density of the fluid under consideration */
  double density;
  /** @brief the kinematic viscosity of the fluid under consideration */
  double kinVisc;
  /** @brief only for LB couette solver: number of processes */
  tarch::la::Vector<3, unsigned int> lbNumberProcesses;
  /** @brief  only for LB couette solver: VTK plotting per time step */
  int plotEveryTimestep;
  /** @brief number of cycles the continuum solver is advanced before the
   * coupling is enabled */
  int initAdvanceCycles;
  /** @brief the start temperature for the fluid under consideration */
  double temp;
  /** @brief number of equilibartion time steps = number of time steps that
   * the md will run before the coupling is enabled */
  int equSteps;
  /** @brief the number of md simulation instances in a multi-instance
   * coupling, -1 = dynamic */
  int totalNumberMDSimulations;
  /** @brief the minimum number of md simulation instances in dynamic MD coupling  */
  int lowerBoundNumberMDSimulations;
  /** @brief only for dynamic MD: target absolute velocity error at start of all coupling cycles */
  double absVelErrStart;
  /** @brief only for dynamic MD: target absolute velocity error at end of all coupling cycles */
  double absVelErrEnd;
  /** @brief the sigma for the random noise in the case of synthetic md solver
   */
  double noiseSigma;
  /** @todo piet */
  bool twsLoop;
  /** @todo piet*/
  int twsLoopMin;
  /** @todo piet */
  int twsLoopMax;
  /** @todo piet */
  int twsLoopStep;

#if (BUILD_WITH_OPENFOAM)
  /** @brief the configurations for the OpenFoam solver */
  FoamConfig foam;
#endif
};

#endif
