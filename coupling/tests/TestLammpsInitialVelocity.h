// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _TESTLAMMPSINITIALVELOCITY_H_
#define _TESTLAMMPSINITIALVELOCITY_H_

#include "TestLammps.h"
#include "tarch/la/Vector.h"

template <unsigned int dim>
class TestLammpsInitialVelocity : public TestLammps<dim> {
public:
  TestLammpsInitialVelocity(int argc, char **argv, std::string name)
      : TestLammps<dim>(argc, argv, name) {}
  virtual ~TestLammpsInitialVelocity() {}

  virtual void run() {
    // initialise all interfaces and simulation parts
    if (dim == 2) {
      TestLammps<dim>::loadLammpsTestConfiguration(
          "inputpositions2D_moleculeiterator.xyz", 24);
    } else {
      TestLammps<dim>::loadLammpsTestConfiguration(
          "inputpositions3D_moleculeiterator.xyz", 96);
    }
    TestLammps<dim>::loadMacroscopicSolverConfiguration();
    TestLammps<dim>::loadMamicoTestConfiguration();
    coupling::interface::MDSolverInterface<LAMMPS_NS::MamicoCell, dim> *
        mdInterface = coupling::interface::MamicoInterfaceProvider<
            LAMMPS_NS::MamicoCell, dim>::getInstance().getMDSolverInterface();

    const double kB = mdInterface->getKB(); // value for kB
    const double temperature = 3.5;         // value for temperature
    const tarch::la::Vector<dim, double> meanVelocity(
        1.0);                               // value for avg velocity
    const double stdDev =
        std::sqrt(dim * kB * temperature / mdInterface->getMoleculeMass());
    const int numberExperiments = 100000;
    const double tol = 1.0e-2;

    std::vector<tarch::la::Vector<dim, double> > velocities;
    tarch::la::Vector<dim, double> expAvg(0.0);
    double expStdDev = 0.0;

    // sample velocities
    for (int i = 0; i < numberExperiments; i++) {
      tarch::la::Vector<dim, double> initVel(0.0);
      mdInterface->getInitialVelocity(meanVelocity, kB, temperature, initVel);
      velocities.push_back(initVel);
    }

    // compute avg
    for (int i = 0; i < numberExperiments; i++) {
      expAvg = expAvg + velocities[i];
    }
    expAvg = (1.0 / numberExperiments) * expAvg;
    // compute std dev
    for (int i = 0; i < numberExperiments; i++) {
      expStdDev = expStdDev + tarch::la::dot(velocities[i] - expAvg,
                                             velocities[i] - expAvg);
    }
    expStdDev = std::sqrt(expStdDev / numberExperiments);

    std::cout << "Average: " << meanVelocity << " measured: " << expAvg
              << ", ||. - .||_2: " << tarch::la::norm2(meanVelocity - expAvg)
              << std::endl;
    std::cout << "Std-dev: " << stdDev << " measured: " << expStdDev
              << ", ||. - .||_2: " << fabs(stdDev - expStdDev) << std::endl;
    if (tarch::la::norm2(meanVelocity - expAvg) > tol) {
      std::cout
          << "ERROR TestLammpsInitialVelocity: Avg. velocity deviation exceeds "
          << tol << "!" << std::endl;
      exit(EXIT_FAILURE);
    }
    if (fabs(stdDev - expStdDev) > tol) {
      std::cout
          << "ERROR TestLammpsInitialVelocity: Std. dev. deviation exceeds "
          << tol << "!" << std::endl;
      exit(EXIT_FAILURE);
    }
  }

};
#endif // _TESTLAMMPSINITIALVELOCITY_H_
