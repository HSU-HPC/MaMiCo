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
#include <functional>

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
    _pint_domain(0), _rank(0), _ranks_per_domain(1), _world_rank(0),
    _cfg( mamicoConfig.getTimeIntegrationConfiguration() ),
    _scenario(scenario) {
        #if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
        int world_size;
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
        MPI_Comm_rank(MPI_COMM_WORLD, &_world_rank);
        if(_cfg.isPinTEnabled()){
            if(world_size % _cfg.getPintDomains() != 0){
                if(_world_rank == 0){
                    std::cout << "ERROR coupling::services::ParallelTimeIntegrationService: " <<
                        "MPI ranks not divisible by number of PinT subdomains!" << std::endl;
                    std::cout << "When PinT is used, the number of required MPI ranks increases by a factor of number-subdomains." << std::endl;
                    std::cout << "Check your configuration." << std::endl;
                }
                exit(EXIT_FAILURE);
            }
            _ranks_per_domain = world_size / _cfg.getPintDomains();
            _pint_domain = _world_rank / _ranks_per_domain;
            // This initializes _local_pint_comm by splitting MPI_COMM_WORLD into getPintDomains() disjoint communicators
            MPI_Comm_split(MPI_COMM_WORLD, _pint_domain, _world_rank, &_local_pint_comm);
        } else {
            _local_pint_comm = MPI_COMM_WORLD;
        }
        MPI_Comm_rank(_local_pint_comm, &_rank);

        #ifdef PINT_DEBUG
        std::cout << "PINT_DEBUG: world_rank " << _world_rank << " is rank " << _rank << " in pint domain " << _pint_domain << std::endl;
        #endif

        #endif
    }

    void run(int num_cycles) {
        if(!_cfg.isPinTEnabled())
            run_cycles(0, num_cycles);
        else {
            PintDomain domain = setup_domain(num_cycles);            
            setup_solvers(domain);
            init_parareal();
            run_parareal( _cfg.getPintIterations() );
        }
    }

    int getPintDomain() { return _pint_domain; }
    int getRank() { return _rank; }

    #if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
    MPI_Comm getPintComm() { return _local_pint_comm; }
    #endif

private:
    void run_cycles(int start, int end){
        for (int cycle = start; cycle < end; cycle++)
            _scenario->runOneCouplingCycle(cycle);
    }

    bool isFirst() { return _pint_domain == 0; }
    bool isLast()  { return _pint_domain == _cfg.getPintDomains()-1; }

    struct PintDomain {
        /** @brief number of time steps (coupling cycles) in this temporal domain */
        int size;
        /** @brief number of the first coupling cycle in this temporal domain (inclusive) */
        int minCycle;
        /** @brief number of the last coupling cycle of this temporal domain (exclusive) */
        int maxCycle;
    };

    PintDomain setup_domain(int num_cycles) {
        PintDomain res;
        res.size = (int)( num_cycles / _cfg.getPintDomains() );
        res.minCycle = _pint_domain * res.size;
        res.maxCycle = res.minCycle + res.size;
        // In case num_cycles is not divisible by _cfg.getPintDomains(), the last domain gets the remainder
        if( isLast() ) 
            res.maxCycle = num_cycles;
        #ifdef PINT_DEBUG
        std::cout << "PINT_DEBUG: _pint_domain " << _pint_domain << " has minCycle " << res.minCycle 
            << " and maxCycle " << res.maxCycle << std::endl;
        #endif
        return res;
    }

    void setup_solvers(PintDomain domain){
        auto solver = dynamic_cast<Solver*>( _scenario->getSolver() );
        if(solver == nullptr){
            std::cout << "ERROR coupling::services::ParallelTimeIntegrationService: " <<
                "macroscopic solver is not pintable (= not compatible with parallel in time coupling)" << std::endl;
            exit(EXIT_FAILURE);
        }
        _G = solver->getSupervisor(domain.size, 2.0);
        _F = [this, &solver, domain](const std::unique_ptr<State>& s){
            solver->setState(s);
            // TODO enable momentumTransfer on inner cells for MD equilibration steps here
            // TODO run MD equilibration here
            run_cycles(domain.minCycle, domain.maxCycle);
            // TODO double check that filter pipeline output ends up in solver.getState() here
            return solver->getState();
        };
        _u_0 = solver->getState();
    }

    void receive(std::unique_ptr<State>& state){
        int source_rank = _world_rank - _ranks_per_domain;
        if(!state) state = _u_0->clone();
        #if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
        MPI_Recv(state->getData(), state->getSizeBytes(), MPI_BYTE, source_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        #endif
    }

    void send(std::unique_ptr<State>& state) {
        int destination_rank = _world_rank + _ranks_per_domain;
        #if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
        MPI_Send(state->getData(), state->getSizeBytes(), MPI_BYTE, destination_rank, 0, MPI_COMM_WORLD);
        #endif
    }

    void get_past_state(std::unique_ptr<State>& state){
        if( isFirst() ) 
            state = _u_0->clone();
        else
            receive(state);
    }

    void init_parareal(){
        get_past_state(_u_last_past);
        auto& G = *_G;
        _u_last_future = G(_u_last_past);
        if( !isLast() )
            send(_u_last_future);
    }

    void run_parareal(int iterations){
        for(int it = 0; it < iterations; it++){
            // Correction step
            auto& G = *_G;
            auto delta = *_F(_u_last_past) - *G(_u_last_past);

            // Prediction step
            get_past_state(_u_next_past);
            auto prediction = G(_u_next_past);

            // Refinement step
            _u_next_future = *prediction + *delta;
            if( !isLast() )
                send(_u_next_future);
            
            // move for next iteration
            _u_last_past   = std::move(_u_next_past);
            _u_last_future = std::move(_u_next_future);
        }

        #ifdef PINT_DEBUG
        std::cout << "PINT_DEBUG: Finished all PinT iterations. " << std::endl;
        #endif
    }

    int _pint_domain;       // the index of the time domain to which this process belongs to
    int _rank;              // rank of current process in _local_pint_comm
    int _ranks_per_domain;  // number of MPI ranks in each time domain
    int _world_rank;        // rank of current process in MPI_COMM_WORLD
    #if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
    MPI_Comm _local_pint_comm; // the communicator of the local time domain of this rank
    #endif
    coupling::configurations::TimeIntegrationConfiguration _cfg;
    Scenario* _scenario;
    std::unique_ptr<Solver> _G;  // The supervisor, i.e. the coarse predictor
    std::function<std::unique_ptr<State>(const std::unique_ptr<State>&)> _F;

    // These objects represent coupled simulation states. There are used by the supervised parallel in time algorithm for operations
    // "last" and "next" describe two consecutive parareal iterations
    // "past" and "future" describe two points in simulation time, one pint time domain apart
    std::unique_ptr<State> _u_0;       // initial state
    std::unique_ptr<State> _u_last_past;
    std::unique_ptr<State> _u_last_future;
    std::unique_ptr<State> _u_next_past;
    std::unique_ptr<State> _u_next_future;
};