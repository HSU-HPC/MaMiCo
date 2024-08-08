// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_TESTS_TESTCELL_H_
#define _MOLECULARDYNAMICS_COUPLING_TESTS_TESTCELL_H_

template <unsigned int dim> class TestCell {
public:
  TestCell() : _b1(tarch::la::Vector<dim, double>(0.0)), _b2(0.0) {}
  ~TestCell() {}

  void setBuffer1(tarch::la::Vector<dim, double> b1) { _b1 = b1; }
  tarch::la::Vector<dim, double> getBuffer1() const { return _b1; }
  void setBuffer2(double b2) { _b2 = b2; }
  double getBuffer2() const { return _b2; }

private:
  tarch::la::Vector<dim, double> _b1;
  double _b2;
};

#endif // _MOLECULARDYNAMICS_COUPLING_TESTS_TESTCELL_H_
