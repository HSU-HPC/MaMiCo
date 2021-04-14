// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _COUPLING_DATASTRUCTURES_DUMMYCELL_H
#define _COUPLING_DATASTRUCTURES_DUMMYCELL_H

#include<vector>
#include<iostream>
#include "coupling/datastructures/Molecule.h"

namespace coupling {
  namespace datastructures {
    class DummyCell;
  }
}

/** describes a dummy cell.
 *  The dimension is hardcoded to 3d
 *  @author Philipp Neumann
 */
class coupling::datastructures::DummyCell {
  public:
    /** initialises linked cell list with numberMolecules empty entries */
    DummyCell(const tarch::la::Vector<3,double> cellOffset=tarch::la::Vector<3,double>(0,0,0)): _molecules(12) {
      _molecules[0].setPosition(tarch::la::Vector<3,double>(0.09438, 2.04958, 0.61647)+cellOffset);
      _molecules[0].setVelocity(tarch::la::Vector<3,double>(-0.70105, -1.03295, 0.65857));
      _molecules[1].setPosition(tarch::la::Vector<3,double>(1.76113, 1.20659, 1.41395)+cellOffset);
      _molecules[1].setVelocity(tarch::la::Vector<3,double>(0.28982, 0.69947, -0.59525));
      _molecules[2].setPosition(tarch::la::Vector<3,double>(1.09300, 0.36040, 1.45549)+cellOffset);
      _molecules[2].setVelocity(tarch::la::Vector<3,double>(0.10770, 1.72253, 0.18966));
      _molecules[3].setPosition(tarch::la::Vector<3,double>(0.09611, 0.83016, 2.13748)+cellOffset);
      _molecules[3].setVelocity(tarch::la::Vector<3,double>(-1.76393, -1.01643, 1.18727));
      _molecules[4].setPosition(tarch::la::Vector<3,double>(0.99965, 1.16250, 0.58097)+cellOffset);
      _molecules[4].setVelocity(tarch::la::Vector<3,double>(-0.19865, 0.28953, -0.58370));
      _molecules[5].setPosition(tarch::la::Vector<3,double>(1.51687, 2.01781, 0.87839)+cellOffset);
      _molecules[5].setVelocity(tarch::la::Vector<3,double>(2.02154, -0.15842, -0.32287));
      _molecules[6].setPosition(tarch::la::Vector<3,double>(2.43428, 2.32386, 1.81244)+cellOffset);
      _molecules[6].setVelocity(tarch::la::Vector<3,double>(1.49034, -1.10203, 0.83891));
      _molecules[7].setPosition(tarch::la::Vector<3,double>(0.81611, 0.14770, 0.19385)+cellOffset);
      _molecules[7].setVelocity(tarch::la::Vector<3,double>(-1.44460, -1.73636, 0.52284));
      _molecules[8].setPosition(tarch::la::Vector<3,double>(2.48619, 0.68831, 0.73015)+cellOffset);
      _molecules[8].setVelocity(tarch::la::Vector<3,double>(-1.59348, 0.39840, -1.14072));
      _molecules[9].setPosition(tarch::la::Vector<3,double>(1.81832, 0.31171, 2.47495)+cellOffset);
      _molecules[9].setVelocity(tarch::la::Vector<3,double>(1.16062, 1.80032, -1.62150));
      _molecules[10].setPosition(tarch::la::Vector<3,double>(1.91300, 1.59869, 2.32975)+cellOffset);
      _molecules[10].setVelocity(tarch::la::Vector<3,double>(0.47142, -0.27172, -0.92108));
      _molecules[11].setPosition(tarch::la::Vector<3,double>(0.92002, 1.67267, 2.01992)+cellOffset);
      _molecules[11].setVelocity(tarch::la::Vector<3,double>(-0.05393, 1.38331, -0.04296));
    }
    ~DummyCell(){
    }

    coupling::datastructures::Molecule<3>& operator[](unsigned int i){return _molecules[i];}

    unsigned int size()const{return 12;}

  private:
    std::vector<coupling::datastructures::Molecule<3>> _molecules;
};

#endif // _COUPLING_DATASTRUCTURES_DUMMYCELL_H
