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

#include <memory>

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
public:
	using State = PintableMacroSolverState;
	enum class Mode { supervising, coupling };

	virtual ~PintableMacroSolver() {}
	virtual std::unique_ptr<State> getState() = 0;
	virtual void setState(const std::unique_ptr<State>&, int cycle = 0) = 0;
	virtual std::unique_ptr<State> operator()(const std::unique_ptr<State>&, int cycle = 0) = 0;
	virtual Mode getMode() const = 0;
	virtual std::unique_ptr<PintableMacroSolver> getSupervisor(int num_cycles, double visc_multiplier=1) const = 0;
	virtual void print(std::ostream& os) const = 0;
	virtual double get_avg_vel(const std::unique_ptr<State>& state) const { return 0; } 
	friend std::ostream& operator<<(std::ostream& os, const coupling::interface::PintableMacroSolver& s) {
		s.print(os);
		return os;
	}
};

/**
 * Interface for state instances
 */
class coupling::interface::PintableMacroSolverState {
public:
	using State = PintableMacroSolverState;
	virtual std::unique_ptr<State> clone() const = 0;
	virtual ~PintableMacroSolverState() {}
	virtual int getSizeBytes() const = 0;   // for MPI communication
	virtual std::unique_ptr<State> operator+(const State&) = 0;
	virtual std::unique_ptr<State> operator-(const State&) = 0;
	virtual bool operator==(const State&) const = 0;
	virtual double* getData() = 0;
	virtual const double* getData() const = 0;
	virtual void print(std::ostream& os) const = 0;
	friend std::ostream& operator<<(std::ostream& os, const coupling::interface::PintableMacroSolverState& s) {
		s.print(os);
		return os;
	}
};