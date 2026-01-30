#pragma once

namespace simplemd {
namespace moleculemappings {
class ConvertForcesFloatToFixedMapping;
}
} // namespace simplemd

/*
 * fixed-point math for force accumulation
 * only active in debug mode, useful for verification of simulation results
 * because results do not depend on order of force summation
 * this expects forces to contain long int and converts them back into double
 */
class simplemd::moleculemappings::ConvertForcesFloatToFixedMapping {
public:
  void beginMoleculeIteration() const {}

  KOKKOS_FUNCTION void handleMolecule(simplemd::Molecule& molecule) const {
    constexpr double maxF = 1e6;
    constexpr double stepF = (double)(std::numeric_limits<long long>::max()) / maxF;

    tarch::la::Vector<MD_DIM, double> force = stepF * molecule.getForce();
    long long fb0{(long long)(force[0])};
    long long fb1{(long long)(force[1])};
    long long fb2{(long long)(force[2])};
    long long& fb0r = fb0;
    long long& fb1r = fb1;
    long long& fb2r = fb2;
    molecule.getForce()[0] = *(double*)(&fb0r);
    molecule.getForce()[1] = *(double*)(&fb1r);
    molecule.getForce()[2] = *(double*)(&fb2r);
  }

  void endMoleculeIteration() const {}

  static const bool IsParallel = true;
};
