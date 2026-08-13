// Copyright (C) 2023 Helmut Schmidt University
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder

#include "coupling/solvers/LBCouetteSolver.h"

using State = coupling::interface::PintableMacroSolverState;

std::unique_ptr<State> coupling::solvers::LBCouetteSolverState::operator+(const State& rhs) {
  const LBCouetteSolverState* other = dynamic_cast<const LBCouetteSolverState*>(&rhs);

#if (COUPLING_MD_ERROR == COUPLING_MD_YES)
  if (other == nullptr) {
    std::cout << "ERROR LBCouetteSolverState operator+ type mismatch" << std::endl;
    exit(EXIT_FAILURE);
  }
  if (other->_pdf.size() != _pdf.size()) {
    std::cout << "ERROR LBCouetteSolverState operator+ size mismatch" << std::endl;
    exit(EXIT_FAILURE);
  }
#endif

  std::unique_ptr<LBCouetteSolverState> res = std::make_unique<LBCouetteSolverState>(*this);
  for (std::vector<double>::size_type i = 0; i < _pdf.size(); ++i)
    res->_pdf[i] += other->_pdf[i];

  return res;
}

std::unique_ptr<State> coupling::solvers::LBCouetteSolverState::operator-(const State& rhs) {
  const LBCouetteSolverState* other = dynamic_cast<const LBCouetteSolverState*>(&rhs);

#if (COUPLING_MD_ERROR == COUPLING_MD_YES)
  if (other == nullptr) {
    std::cout << "ERROR LBCouetteSolverState operator- type mismatch" << std::endl;
    exit(EXIT_FAILURE);
  }
  if (other->_pdf.size() != _pdf.size()) {
    std::cout << "ERROR LBCouetteSolverState operator- size mismatch" << std::endl;
    exit(EXIT_FAILURE);
  }
#endif

  std::unique_ptr<LBCouetteSolverState> res = std::make_unique<LBCouetteSolverState>(*this);
  for (std::vector<double>::size_type i = 0; i < _pdf.size(); ++i)
    res->_pdf[i] -= other->_pdf[i];

  return res;
}
