// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_CONFIGURATIONS_CHECKPOINTCONFIGURATION_H_
#define _MOLECULARDYNAMICS_CONFIGURATIONS_CHECKPOINTCONFIGURATION_H_

#include "simplemd/configurations/VTKConfiguration.h"

namespace simplemd {
namespace configurations { class CheckpointConfiguration; }
}

/** data for checkpoint writing */
class simplemd::configurations::CheckpointConfiguration
    : public simplemd::configurations::VTKConfiguration {
public:
  CheckpointConfiguration() : simplemd::configurations::VTKConfiguration() {}
  virtual ~CheckpointConfiguration() {}

  virtual std::string getTag() const { return "checkpoint-configuration"; }
};
#endif // _MOLECULARDYNAMICS_CONFIGURATIONS_CHECKPOINTCONFIGURATION_H_
