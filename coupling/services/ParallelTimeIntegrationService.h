// Copyright (C) 2023 Helmut Schmidt University
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder

#pragma once

namespace coupling {
namespace services {
template <class LinkedCell, unsigned int dim> class ParallelTimeIntegrationService;
} // namespace services
} // namespace coupling

#include "coupling/configurations/MaMiCoConfiguration.h"
#include "coupling/interface/PintableMacroSolver.h"

/** 
 * Service to manage timeloop of a coupled simulation scenario. Supports sequential or parallel-in-time integration using a Parareal variant,
 * as described in "Blumers, A. L., Li, Z., & Karniadakis, G. E. (2019). Supervised parallel-in-time algorithm for long-time Lagrangian 
 * simulations of stochastic dynamics: Application to hydrodynamics. Journal of Computational Physics, 393, 214-228".
 * 
 *  @author Piet Jarmatz
 */
template <class LinkedCell, unsigned int dim> class coupling::services::ParallelTimeIntegrationService {
public:
    ParallelTimeIntegrationService(coupling::configurations::MaMiCoConfiguration<dim> mamicoConfig):
    _pint_domain(1),
    _cfg( mamicoConfig.getTimeIntegrationConfiguration() ) {
        #if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
        int world_size, world_rank;
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
        if(_cfg.isPinTEnabled()){
            if(world_size % _cfg.getPintDomains() != 0){
                if(world_rank == 0){
                    std::cout << "ERROR coupling::services::ParallelTimeIntegrationService: " <<
                        "MPI ranks not divisible by number of PinT subdomains!" << std::endl;
                    std::cout << "When PinT is used, the number of required MPI ranks increases by a factor of number-subdomains." << std::endl;
                    std::cout << "Check your configuration." << std::endl;
                }
                exit(EXIT_FAILURE);
            }
            int ranks_per_domain = world_size / _cfg.getPintDomains();
            _pint_domain = world_rank / ranks_per_domain;
            // This initializes _local_pint_comm by splitting MPI_COMM_WORLD into getPintDomains() disjoint communicators
            MPI_Comm_split(MPI_COMM_WORLD, _pint_domain, world_rank, &_local_pint_comm);
        } else {
            _local_pint_comm = MPI_COMM_WORLD;
        }
        MPI_Comm_rank(_local_pint_comm, &_rank);

        #ifdef PINT_DEBUG
        std::cout << "PINT_DEBUG: world_rank " << world_rank << " is rank " << _rank << " in pint domain " << _pint_domain << std::endl;
        #endif

        #endif
    }

    int getPintDomain() { return _pint_domain; }
    int getRank() { return _rank; }

    #if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
    MPI_Comm getPintComm() { return _local_pint_comm; }
    #endif

private:
    int _pint_domain;  // the index of the time domain to which this process belongs to
    int _rank;         // rank of current process in _local_pint_comm
    #if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
    MPI_Comm _local_pint_comm; // the communicator of the local time domain of this rank
    #endif
    coupling::configurations::TimeIntegrationConfiguration _cfg;
};