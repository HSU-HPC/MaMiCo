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

  KOKKOS_FUNCTION void beginMoleculeIteration();

  KOKKOS_FUNCTION void endMoleculeIteration() { Kokkos::Profiling::popRegion(); }
  KOKKOS_FUNCTION void handleMolecule(Molecule& molecule, LinkedCell& cell) const;

  static const bool IsParallel = true;

};
