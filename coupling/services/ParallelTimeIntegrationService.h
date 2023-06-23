// Copyright (C) 2023 Helmut Schmidt University
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder

#pragma once

namespace coupling {
namespace services {
template <class LinkedCell, unsigned int dim> class ParallelTimeIntegrationService;
} // namespace services
} // namespace coupling

//#include "coupling/configurations/MaMiCoConfiguration.h"
//#include "coupling/indexing/IndexingService.h"
#include "coupling/services/MultiMDCellService.h"

/** 
 * Service to manage timeloop of a coupled simulation scenario. Supports sequential or parallel-in-time integration using a Parareal variant,
 * as described in "Blumers, A. L., Li, Z., & Karniadakis, G. E. (2019). Supervised parallel-in-time algorithm for long-time Lagrangian 
 * simulations of stochastic dynamics: Application to hydrodynamics. Journal of Computational Physics, 393, 214-228".
 * 
 *  @author Piet Jarmatz
 */
template <class LinkedCell, unsigned int dim> class coupling::services::ParallelTimeIntegrationService {
public:
    ParallelTimeIntegrationService() {}

private:

};