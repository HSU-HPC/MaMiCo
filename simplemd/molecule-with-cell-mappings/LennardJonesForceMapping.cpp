#include "simplemd/molecule-with-cell-mappings/LennardJonesForceMapping.h"

simplemd::moleculewithcellmappings::LennardJonesForceMapping::LennardJonesForceMapping(simplemd::services::ExternalForceService& externalForceService,
                                                                           const simplemd::services::MolecularPropertiesService& molecularPropertiesService)
    : simplemd::cellmappings::LennardJonesForceMapping(externalForceService,molecularPropertiesService) { }

void simplemd::moleculewithcellmappings::LennardJonesForceMapping::beginMoleculeIteration() {
#if (MD_DEBUG == MD_YES)
  Kokkos::printf("simplemd::moleculewithcellmappings::LennardJonesForceMapping::beginMoleculeIteration()\n");
#endif
  Kokkos::Profiling::pushRegion("simplemd::moleculewithcellmappings::LennardJonesForceMapping");
}

/*
 * uni-directional version, i.e. Newton3 OFF
 * fixed-point math for force accumulation
 * only active in debug mode, useful for verification of simulation results
 * because results do not depend on order of force summation
 * this expects target and forceBuffer to contain correctly formatted long long data already, not double
 */
KOKKOS_INLINE_FUNCTION void addForce(tarch::la::Vector<MD_DIM, double>& target,
                                     tarch::la::Vector<MD_DIM, double>& forceBuffer) {
#if (TARCH_DEBUG == TARCH_YES)
  *(long long*)(&target[0]) += *(long long*)(&forceBuffer[0]);
  *(long long*)(&target[1]) += *(long long*)(&forceBuffer[1]);
  *(long long*)(&target[2]) += *(long long*)(&forceBuffer[2]);
#else
  target += forceBuffer;
#endif
}

void simplemd::moleculewithcellmappings::LennardJonesForceMapping::handleMolecule(Molecule& m1, const LinkedCell& cell) const {

	tarch::la::Vector<MD_DIM, double> forceBuffer(0.0); 
	tarch::la::Vector<MD_DIM, double>& force1 = m1.getForce();
    const tarch::la::Vector<MD_DIM, double>& position1 = m1.getConstPosition();

#if (TARCH_DEBUG == TARCH_YES)
    if (_externalForce != tarch::la::Vector<MD_DIM, double>{0.0}) {
      Kokkos::abort("ERROR simplemd::moleculewithcellmappings::LennardJonesForceMapping::handleCell(): externalForce not implemented in fixed point math debug mode!\n");
    }
#else
    force1 += _externalForce;
#endif

    const auto end = cell.end();
    const auto begin = cell.begin();
    for (auto m2 = begin; m2 != end; m2++) {
    	forceBuffer = getLennardJonesForce(position1, m2->getConstPosition());
#if (MD_DEBUG == MD_YES)
      if (tarch::la::dot(forceBuffer, forceBuffer) > 1e12) {
        const auto position2 = m2->getConstPosition();
        Kokkos::printf("Force: %lf %lf %lf; "
                       "Position1: %lf %lf %lf; "
                       "Position2: %lf %lf %lf; "
                       "ID1: %u;"
                       "ID2: %u"
                       "\n",
                       forceBuffer[0], forceBuffer[1], MD_DIM3_OR0(forceBuffer[2]), position1[0], position1[1], MD_DIM3_OR0(position1[2]), position2[0],
                       position2[1], MD_DIM3_OR0(position2[2]), m1.getID(), m2->getID());
        Kokkos::abort("ERROR simplemd::moleculewithcellmappings::LennardJonesForceMapping::handleCellPair: Force out of range!\n");
      }
#endif
      addForce(force1, forceBuffer);
    }
}