// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _COUPLING_CELLMAPPINGS_VARYCHECKPOINTMAPPING_H_
#define _COUPLING_CELLMAPPINGS_VARYCHECKPOINTMAPPING_H_

#include "coupling/interface/MDSolverInterface.h"

namespace coupling {
  namespace cellmappings {
    template<class LinkedCell, unsigned int dim> class VaryCheckpointMapping;
  }
}

template<class LinkedCell, unsigned int dim>
class coupling::cellmappings::VaryCheckpointMapping {
public:
  VaryCheckpointMapping(const coupling::interface::MDSolverInterface<LinkedCell, dim> * const mdSolverInterface, double& molecularMass, const double & kB, const double & temperature, const double & sigma)
    : _mdSolverInterface(mdSolverInterface), _molecularMass(molecularMass), _kB(kB), _temperature(temperature), _sigma(sigma) {}

  ~VaryCheckpointMapping() {}

  void beginCellIteration() {}
  void endCellIteration() {}

  void handleCell(LinkedCell &cell,const unsigned int &cellIndex) {
    tarch::la::Vector<MD_DIM, double> meanVelocityForCell(0.0);
    double stdDeviation = std::sqrt(MD_DIM * _kB * _temperature / _molecularMass);
    
    tarch::la::Vector<MD_DIM, double> randomNumbers(0.0);
    

    for(coupling::interface::MoleculeIterator<LinkedCell,dim> *molecule = _mdSolverInterface->getMoleculeIterator(cell);
          molecule.continueIteration(); 
          molecule->next()
    ) {
      meanVelocityForCell += molecule->get().getVelocity();
    }
    meanVelocityForCell = meanVelocityForCell / (double)cell.getConstList().size() / _molecularMass;

    for(coupling::interface::MoleculeIterator<LinkedCell,dim> *molecule = _mdSolverInterface->getMoleculeIterator(cell);
          molecule.continueIteration(); 
          molecule->next()
    ) {
      randomNumbers[0] = tarch::utils::RandomNumberService::getInstance().getGaussianRandomNumber();
      for(unsigned int d=1;d<MD_DIM;++d) {
        randomNumbers[d] = tarch::utils::RandomNumberService::getInstance().getGaussianRandomNumber();
      }

      tarch::la::Vector<MD_DIM, double> & mVelocity = molecule->get().getVelocity();
#if (MD_DIM==1)
      mVelocity = meanVelocityForCell + stdDeviation*randomNumbers;
#elif (MD_DIM==2)
      mVelocity[0] = meanVelocityForCell[0] + stdDeviation*(randomNumbers[0]*std::cos(randomNumbers[1]));
      mVelocity[1] = meanVelocityForCell[1] + stdDeviation*(randomNumbers[0]*std::sin(randomNumbers[1]));
#elif (MD_DIM==3)
      mVelocity[0] = meanVelocityForCell[0] + stdDeviation*(randomNumbers[0]*std::sin(randomNumbers[1])*std::cos(randomNumbers[2]));
      mVelocity[1] = meanVelocityForCell[1] + stdDeviation*(randomNumbers[0]*std::sin(randomNumbers[1])*std::sin(randomNumbers[2]));
      mVelocity[2] = meanVelocityForCell[2] + stdDeviation*(randomNumbers[0]*std::cos(randomNumbers[1]));
#endif
    } 
  }


private:
  coupling::interface::MDSolverInterface<LinkedCell,dim> * const _mdSolverInterface;
  const double _molecularMass;
  const double _kB;
  const double _temperature;
  const double _sigma;
};


#endif //_COUPLING_CELLMAPPINGS_VARYCHECKPOINTMAPPING_H_