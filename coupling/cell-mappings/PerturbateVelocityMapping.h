// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _COUPLING_CELLMAPPINGS_PERTURBATEVELOCITYMAPPING_H_
#define _COUPLING_CELLMAPPINGS_PERTURBATEVELOCITYMAPPING_H_

#include "coupling/interface/MDSolverInterface.h"

namespace coupling {
  namespace cellmappings {
    template<class LinkedCell, unsigned int dim> class PerturbateVelocityMapping;
  }
}

template<class LinkedCell, unsigned int dim>
class coupling::cellmappings::PerturbateVelocityMapping {
public:
  PerturbateVelocityMapping(coupling::interface::MDSolverInterface<LinkedCell, dim> * const mdSolverInterface, const tarch::la::Vector<dim, double> & velocity, const double & temperature)
    : _mdSolverInterface(mdSolverInterface), _molecularMass(mdSolverInterface->getMoleculeMass()), _kB(mdSolverInterface->getKB()), 
      _temperature(temperature), _velocity(velocity), _sigma(_mdSolverInterface->getMoleculeSigma()) {}

  ~PerturbateVelocityMapping() {}

  void beginCellIteration() {}
  void endCellIteration() {}

  void handleCell(LinkedCell &cell,const unsigned int &cellIndex) {
    
    double stdDeviation = std::sqrt(MD_DIM * _kB * _temperature / _molecularMass);
    
    tarch::la::Vector<MD_DIM, double> randomNumbers(0.0);

    for(coupling::interface::MoleculeIterator<LinkedCell,dim> *molecule = _mdSolverInterface->getMoleculeIterator(cell);
          molecule->continueIteration(); 
          molecule->next()
    ) {
      randomNumbers[0] = tarch::utils::RandomNumberService::getInstance().getGaussianRandomNumber();
      for(unsigned int d=1;d<MD_DIM;++d) {
        randomNumbers[d] = tarch::utils::RandomNumberService::getInstance().getGaussianRandomNumber();
      }

      tarch::la::Vector<MD_DIM, double> mVelocity = molecule->get().getVelocity();
#if (MD_DIM==1)
      mVelocity = _velocity + stdDeviation*randomNumbers;
#elif (MD_DIM==2)
      mVelocity[0] = _velocity[0] + stdDeviation*(randomNumbers[0]*std::cos(randomNumbers[1]));
      mVelocity[1] = _velocity[1] + stdDeviation*(randomNumbers[0]*std::sin(randomNumbers[1]));
#elif (MD_DIM==3)
      mVelocity[0] = _velocity[0] + stdDeviation*(randomNumbers[0]*std::sin(randomNumbers[1])*std::cos(randomNumbers[2]));
      mVelocity[1] = _velocity[1] + stdDeviation*(randomNumbers[0]*std::sin(randomNumbers[1])*std::sin(randomNumbers[2]));
      mVelocity[2] = _velocity[2] + stdDeviation*(randomNumbers[0]*std::cos(randomNumbers[1]));
#endif
    molecule->get().setVelocity(mVelocity);
    } 
  }


private:
  coupling::interface::MDSolverInterface<LinkedCell,dim> * const _mdSolverInterface;
  const double _molecularMass;
  const double _kB;
  const double _temperature;
  const tarch::la::Vector<dim, double> _velocity;
  const double _sigma;
};


#endif //_COUPLING_CELLMAPPINGS_VARYCHECKPOINTMAPPING_H_