// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_MPIINDEXCONVERSION_H_
#define _MOLECULARDYNAMICS_COUPLING_MPIINDEXCONVERSION_H_

namespace coupling { template <unsigned int dim> class MPIIndexConversion; }

#include "coupling/IndexConversion.h"

template <unsigned int dim>
class coupling::MPIIndexConversion : public coupling::IndexConversion {}

#endif