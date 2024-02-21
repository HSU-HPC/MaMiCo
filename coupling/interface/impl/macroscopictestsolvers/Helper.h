// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_INTERFACE_MACROSCOPICTESTSOLVERSHELPER_H_
#define _MOLECULARDYNAMICS_COUPLING_INTERFACE_MACROSCOPICTESTSOLVERSHELPER_H_

#include <cmath>

namespace coupling {
namespace interface {
class Helper;
}
} // namespace coupling

class coupling::interface::Helper {
public:
  template <unsigned int dim>
  static tarch::la::Vector<dim, unsigned int> getGlobalNumberCouplingCells(tarch::la::Vector<dim, double> mdDomainSize,
                                                                           tarch::la::Vector<dim, double> couplingCellSize) {
    tarch::la::Vector<dim, unsigned int> numberCouplingCells(0);
    for (unsigned int d = 0; d < dim; d++) {
      numberCouplingCells[d] = (unsigned int)floor(mdDomainSize[d] / couplingCellSize[d] + 0.5);
    }
    return numberCouplingCells;
  }
};
#endif // _MOLECULARDYNAMICS_COUPLING_INTERFACE_MACROSCOPICTESTSOLVERSHELPER_H_
