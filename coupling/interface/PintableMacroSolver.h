// Copyright (C) 2023 Helmut Schmidt University
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder

#pragma once

namespace coupling {
namespace interface {
class PintableMacroSolver;
class PintableMacroSolverState;
}
} 

/** 
 *  Interface for a macroscopic solver who supports parallel-in-time coupled simulations. This requires:
 *   - The macro solver can be switched between *supervising mode* and *coupling mode*
 *   - In supervising mode, the macro solver runs everywhere (incl inner MD cells)
 *   - In coupling mode, the macro solver runs only on the CFD domain (outer MD cells and outside MaMiCo)
 *   - Parameterization of the solver (e.g. viscosity, relaxation time) can differ between supervising and coupling mode
 *   - The macro solver can output / return a full state object, containing all data describing the CFD state on this ranks subdomain
 * 	 - In supervising mode, the macro solver can be non-zero-initialised to a given state (thus it can not be analytical)
 *   - Simulation state instances can perform the arithmetic operations addition, subtraction, copy assignment
 * 
 * 	All macroscopic solvers which want to support PinT have to implement this interface.
 * 
 *  @author Piet Jarmatz
 */
class coupling::interface::PintableMacroSolver {
	using State = PintableMacroSolverState;
	enum class Mode { supervising, coupling };
public:
	virtual ~PintableMacroSolver() = 0;
	virtual std::unique_ptr<State> getInitialState() = 0;
	virtual std::unique_ptr<State> operator()(const State&) = 0;
	virtual Mode getMode() = 0;
};

/**
 * Interface for state instances
 */
class coupling::interface::PintableMacroSolverState {
	using State = PintableMacroSolverState;
public:
	virtual ~PintableMacroSolverState() = 0;
	virtual long getSizeBytes() = 0;   // for MPI communication
	virtual std::unique_ptr<State> operator+(const State&) = 0;
	virtual std::unique_ptr<State> operator-(const State&) = 0;
};
