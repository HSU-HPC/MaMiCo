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
#include "coupling/scenario/Scenario.h"

/** 
 * Service to manage timeloop of a coupled simulation scenario. Supports sequential or parallel-in-time integration using a Parareal variant,
 * as described in "Blumers, A. L., Li, Z., & Karniadakis, G. E. (2019). Supervised parallel-in-time algorithm for long-time Lagrangian 
 * simulations of stochastic dynamics: Application to hydrodynamics. Journal of Computational Physics, 393, 214-228".
 * 
 *  @author Piet Jarmatz
 */
template <class LinkedCell, unsigned int dim> class coupling::services::ParallelTimeIntegrationService {
public:
    using State = coupling::interface::PintableMacroSolverState;
    using Solver = coupling::interface::PintableMacroSolver;

    ParallelTimeIntegrationService(coupling::configurations::MaMiCoConfiguration<dim> mamicoConfig, Scenario* scenario):
    _pint_domain(0), _rank(0), 
    _cfg( mamicoConfig.getTimeIntegrationConfiguration() ),
    _scenario(scenario) {
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

    void run(int num_cycles) {
        if(!_cfg.isPinTEnabled()){
            for (int cycle = 0; cycle < num_cycles; cycle++)
                _scenario->runOneCouplingCycle(cycle);
        } 
        else {
            // Domain setup //////////////////////////////////////////////////////////////////////////
            int pintDomainSize = (int)( num_cycles / _cfg.getPintDomains() );
            int minCycle{ _pint_domain * pintDomainSize };
            int maxCycle{ minCycle + pintDomainSize };
            // In case num_cycles is not divisible by _cfg.getPintDomains(), the last domain gets the remainder
            if(_pint_domain == _cfg.getPintDomains()-1) 
                maxCycle = num_cycles;
            #ifdef PINT_DEBUG
            std::cout << "PINT_DEBUG: _pint_domain " << _pint_domain << " has minCycle " << minCycle << " and maxCycle " << maxCycle << std::endl;
            #endif

            // Solver setup  ////////////////////////////////////////////////////////////////////////////////
            auto solver = dynamic_cast<Solver*>( _scenario->getSolver() );
            if(solver == nullptr){
                std::cout << "ERROR coupling::services::ParallelTimeIntegrationService: " <<
                    "macroscopic solver is not pintable (= not compatible with parallel in time coupling)" << std::endl;
                exit(EXIT_FAILURE);
            }
            _G = solver->getSupervisor(pintDomainSize, 2.0);
            auto F = [](State s&){
                solver->setState(s);
                // TODO enable momentumTransfer on inner cells for MD equilibration steps here
                // TODO run MD equilibration here
                for (int cycle = minCycle; cycle < maxCycle; cycle++)
                    _scenario->runOneCouplingCycle(cycle);
                // TODO double check that filter pipeline output ends up in solver.getState() here
                return solver.getState();
            };
            std::unique_ptr<State> u_0 = solver->getState();

            // Initialisation  /////////////////////////////////////////////////////////////////////
            if(_pint_domain == 0) 
                _u_last_past = u_0->clone();
            else
                receive(_u_last_past);
            _u_last_future = _G->operator()(*_u_last_past);
            if(_pint_domain != _cfg.getPintDomains()-1)
                send(_u_last_future);

            // Parareal iterations   //////////////////////////////////////////////////////////////////
            for(int it = 0; it < _cfg.getPintIterations(); it++){
                // correction
                std::unique_ptr<State> delta = F(_u_last_past) - _G->operator()(*_u_last_past);

                // prediction
                if(_pint_domain == 0)
                    _u_next_past = u_0->clone();
                else
                    receive(_u_next_past);
                std::unique_ptr<State> prediction = _G->operator()(*_u_next_past);

                // refinement
                _u_next_future = prediction + delta;
                if(_pint_domain != _cfg.getPintDomains()-1)
                    send(_u_next_future);
                
                // swap for next iteration
                _u_last_past = _u_next_past->clone();
                _last_future = _u_next_future->clone();
            }

            #ifdef PINT_DEBUG
            std::cout << "PINT_DEBUG: Finished all PinT iterations. " << std::endl;
            #endif
        }
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
    Scenario* _scenario;
    std::unique_ptr<Solver> _G;  // The supervisor, i.e. the coarse predictor

    // These objects represent coupled simulation states. There are used by the supervised parallel in time algorithm for operations
    // "last" and "next" describe two consecutive parareal iterations
    // "past" and "future" describe two points in simulation time, one pint time domain apart
    std::unique_ptr<State> _u_last_past;
    std::unique_ptr<State> _u_last_future;
    std::unique_ptr<State> _u_next_past;
    std::unique_ptr<State> _u_next_future;
};