#pragma once

#include "simplemd/cell-mappings/LennardJonesForceMapping.h"

#include <Kokkos_Core.hpp>

namespace simplemd {
namespace moleculewithcellmappings {
class LennardJonesForceMapping;
}
} // namespace simplemd

class simplemd::moleculewithcellmappings::LennardJonesForceMapping : public simplemd::cellmappings::LennardJonesForceMapping {
public:
  LennardJonesForceMapping(simplemd::services::ExternalForceService& externalForceService,
                           const simplemd::services::MolecularPropertiesService& molecularPropertiesService);
  KOKKOS_FUNCTION virtual ~LennardJonesForceMapping() {}

  void beginMoleculeIteration(const Kokkos::View<double**[3], Kokkos::LayoutRight>& posData,
    const Kokkos::View<int*, Kokkos::LayoutRight>& linkedCellNumMolecules);
  void endMoleculeIteration() { }

  KOKKOS_FUNCTION void handleMolecule(Molecule& molecule, const LinkedCell& cell) const;
  KOKKOS_FUNCTION void handleMoleculeVeryFast(Molecule& molecule, int cellIndex) const;
  KOKKOS_FUNCTION void handleMoleculeVeryFastES1(Molecule& molecule, int cellIndex) const;

  static const bool IsParallel = true;
  static const bool IsReadonly = false;

private:
  Kokkos::View<const double**[3], Kokkos::LayoutRight> _posData;
  Kokkos::View<const int*, Kokkos::LayoutRight> _linkedCellNumMolecules;
};
