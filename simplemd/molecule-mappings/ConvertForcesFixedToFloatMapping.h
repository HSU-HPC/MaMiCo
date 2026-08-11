#ifndef _MOLECULARDYNAMICS_MOLECULEMAPPINGS_CONVERTFORCES_MAPPING_H_
#define _MOLECULARDYNAMICS_MOLECULEMAPPINGS_CONVERTFORCES_MAPPING_H_

#include "tarch/utils/Utils.h"

namespace simplemd {
namespace moleculemappings {
class ConvertForcesFixedToFloatMapping;
}
} // namespace simplemd

/*
 * fixed-point math for force accumulation
 * only active in debug mode, useful for verification of simulation results
 * because results do not depend on order of force summation
 * this expects forces to contain long int and converts them back into double
 */
class simplemd::moleculemappings::ConvertForcesFixedToFloatMapping {
public:
  void beginMoleculeIteration() const {}
  KOKKOS_FUNCTION void handleMolecule(simplemd::Molecule& molecule) const {
    DEFINE_DECIMAL_FP_LIMITS(6);

    tarch::la::Vector<MD_DIM, double>& force = molecule.getForce();
    for (unsigned int d = 0; d < MD_DIM; d++) {
      force[d] = *(long long*)(&force[d]) * minFP6;
    }
  }
  void endMoleculeIteration() const {}
  static const bool IsParallel = true;
};

#endif // _MOLECULARDYNAMICS_MOLECULEMAPPINGS_CONVERTFORCES_MAPPING_H_
