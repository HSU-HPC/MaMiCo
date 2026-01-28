#ifndef _MOLECULARDYNAMICS_MOLECULEMAPPINGS_CONVERTFORCES_MAPPING_H_
#define _MOLECULARDYNAMICS_MOLECULEMAPPINGS_CONVERTFORCES_MAPPING_H_

namespace simplemd {
namespace moleculemappings {
class ConvertForcesMapping;
}
} // namespace simplemd

/*
 * fixed-point math for force accumulation
 * only active in debug mode, useful for verification of simulation results
 * because results do not depend on order of force summation
 * this expects forces to contain long int and converts them back into double
 */
class simplemd::moleculemappings::ConvertForcesMapping {
public:
  void beginMoleculeIteration() const {}
  KOKKOS_FUNCTION void handleMolecule(simplemd::Molecule& molecule) const {
    constexpr double maxF = 1e6;
    constexpr double stepF = (double)(std::numeric_limits<long long>::max()) / maxF;
    constexpr double minF = 1 / stepF;

    tarch::la::Vector<MD_DIM, double>& force = molecule.getForce();
    for (unsigned int d = 0; d < MD_DIM; d++) {
      force[d] = *(long long*)(&force[d]) * minF;
    }
  }
  void endMoleculeIteration() const {}
  static const bool IsParallel = true;
};

#endif // _MOLECULARDYNAMICS_MOLECULEMAPPINGS_CONVERTFORCES_MAPPING_H_
