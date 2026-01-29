// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#include "simplemd/cell-mappings/LennardJonesForceMapping.h"
#include <limits>

simplemd::cellmappings::LennardJonesForceMapping::LennardJonesForceMapping(simplemd::services::ExternalForceService& externalForceService,
                                                                           const simplemd::services::MolecularPropertiesService& molecularPropertiesService)
    : _epsilon(molecularPropertiesService.getMolecularProperties().getEpsilon()),
      _sigma6(molecularPropertiesService.getMolecularProperties().getSigma() * molecularPropertiesService.getMolecularProperties().getSigma() *
              molecularPropertiesService.getMolecularProperties().getSigma() * molecularPropertiesService.getMolecularProperties().getSigma() *
              molecularPropertiesService.getMolecularProperties().getSigma() * molecularPropertiesService.getMolecularProperties().getSigma()),
      _cutOffRadiusSquared(molecularPropertiesService.getMolecularProperties().getCutOffRadius() *
                           molecularPropertiesService.getMolecularProperties().getCutOffRadius()),
      _externalForce(0) {
  // Note: This copies the (static) external force from the configuration
  externalForceService.addExternalForce(_externalForce);
}

void simplemd::cellmappings::LennardJonesForceMapping::beginCellIteration() {
#if (MD_DEBUG == MD_YES)
  Kokkos::printf("simplemd::cellmappings::LennardJonesForceMapping::beginCellIteration()\n");
#endif
}

/*
 * fixed-point math for force accumulation
 * only active in debug mode, useful for verification of simulation results
 * because results do not depend on order of force summation
 * this expects force1, force2 and forceBuffer to contain correctly formatted long long data already, not double
 */
KOKKOS_INLINE_FUNCTION void addForce(tarch::la::Vector<MD_DIM, double>& force1, tarch::la::Vector<MD_DIM, double>& force2,
                                     tarch::la::Vector<MD_DIM, double>& forceBuffer) {
#if (TARCH_DEBUG == TARCH_YES)
  *(long long*)(&force1[0]) += *(long long*)(&forceBuffer[0]);
  *(long long*)(&force1[1]) += *(long long*)(&forceBuffer[1]);
  *(long long*)(&force1[2]) += *(long long*)(&forceBuffer[2]);
  *(long long*)(&force2[0]) -= *(long long*)(&forceBuffer[0]);
  *(long long*)(&force2[1]) -= *(long long*)(&forceBuffer[1]);
  *(long long*)(&force2[2]) -= *(long long*)(&forceBuffer[2]);
#else
  force1 += forceBuffer;
  force2 -= forceBuffer;
#endif
}

void simplemd::cellmappings::LennardJonesForceMapping::handleCell(LinkedCell& cell) const {
  // force buffer
  tarch::la::Vector<MD_DIM, double> forceBuffer(0.0);

  // iterate over all molecules
  const auto end = cell.end();
  const auto begin = cell.begin();
  for (auto m1 = begin; m1 != end; m1++) {
    auto m2 = m1;
    tarch::la::Vector<MD_DIM, double>& force1 = m1->getForce();
    const tarch::la::Vector<MD_DIM, double>& position1 = m1->getConstPosition();

#if (TARCH_DEBUG == TARCH_YES)
    if (_externalForce != tarch::la::Vector<MD_DIM, double>{0.0}) {
      Kokkos::abort("ERROR simplemd::cellmappings::LennardJonesForceMapping::handleCell(): externalForce not implemented in fixed point math debug mode!");
    }
#else
    force1 += _externalForce;
#endif

    // iterate over all other molecules not touched so far
    m2++;
    while (m2 != end) {
      tarch::la::Vector<MD_DIM, double>& force2 = m2->getForce();
      forceBuffer = getLennardJonesForce(position1, m2->getConstPosition());
#if (MD_DEBUG == MD_YES)
      if (tarch::la::dot(forceBuffer, forceBuffer) > 1e12) {
        const auto position2 = m2->getConstPosition();
        Kokkos::printf(
          "Force: %lf %lf %lf; "
          "Position1: %lf %lf %lf; "
          "Position2: %lf %lf %lf; "
          "ID1: %u;"
          "ID2: %u"
          "\n",
          forceBuffer[0], forceBuffer[1], MD_DIM > 2 ? forceBuffer[2] : 0,
          position1[0], position1[1], MD_DIM > 2 ? position1[2] : 0,
          position2[0], position2[1], MD_DIM > 2 ? position2[2] : 0,
          m1->getID(),
          m2->getID()
        );
        Kokkos::abort("ERROR simplemd::cellmappings::LennardJonesForceMapping::handleCellPair: Force out of range!");
      }
#endif
      addForce(force1, force2, forceBuffer);
      m2++;
    }
  }
}

void simplemd::cellmappings::LennardJonesForceMapping::handleCellPair(const LinkedCell& cell1, const LinkedCell& cell2, const unsigned int& cellIndex1,
                                                                      const unsigned int& cellIndex2) const {

  // force buffer
  tarch::la::Vector<MD_DIM, double> forceBuffer(0.0);

  // iterate over pairs of molecules
  const auto endCell1 = cell1.end();
  const auto endCell2 = cell2.end();
  const auto beginCell1 = cell1.begin();
  const auto beginCell2 = cell2.begin();
  for (auto m1 = beginCell1; m1 != endCell1; m1++) {
    tarch::la::Vector<MD_DIM, double>& force1 = m1->getForce();
    const tarch::la::Vector<MD_DIM, double>& position1 = m1->getConstPosition();

    for (auto m2 = beginCell2; m2 != endCell2; m2++) {
      tarch::la::Vector<MD_DIM, double>& force2 = m2->getForce();
      forceBuffer = getLennardJonesForce(position1, m2->getConstPosition());
#if (MD_DEBUG == MD_YES)
      if (tarch::la::dot(forceBuffer, forceBuffer) > 1e12) {
        const auto position2 = m2->getConstPosition();
        Kokkos::printf(
          "Force: %lf %lf %lf; "
          "Position1: %lf %lf %lf; "
          "Position2: %lf %lf %lf; "
          "ID1: %u;"
          "ID2: %u"
          "\n",
          forceBuffer[0], forceBuffer[1], MD_DIM > 2 ? forceBuffer[2] : 0,
          position1[0], position1[1], MD_DIM > 2 ? position1[2] : 0,
          position2[0], position2[1], MD_DIM > 2 ? position2[2] : 0,
          m1->getID(),
          m2->getID()
        );
        Kokkos::abort("ERROR simplemd::cellmappings::LennardJonesForceMapping::handleCellPair: Force out of range!");
      }
#endif
      addForce(force1, force2, forceBuffer);
    }
  }
}

tarch::la::Vector<MD_DIM, double>
simplemd::cellmappings::LennardJonesForceMapping::getLennardJonesForce(const tarch::la::Vector<MD_DIM, double>& position1,
                                                                       const tarch::la::Vector<MD_DIM, double>& position2) const {
  const tarch::la::Vector<MD_DIM, double> rij(position2 - position1);
  const double rij2 = tarch::la::dot(rij, rij);
#if (MD_ERROR == MD_YES)
  if (tarch::la::equals(rij2, 0.0, 1e-4)) {
    Kokkos::printf("Position: %lf %lf %lf, %lf %lf %lf\n", position1[0], position1[1], MD_DIM > 2 ? position1[2] : 0, position2[0], position2[1],
                   MD_DIM > 2 ? position2[2] : 0);
    Kokkos::abort("ERROR simplemd::cellmappings::LennardJonesForceMapping::getLennardJonesForce(): Particle positions are identical!");
  }
#endif

  if (rij2 <= _cutOffRadiusSquared) {
    const double rij6 = rij2 * rij2 * rij2;
#if (TARCH_DEBUG == TARCH_YES)
    tarch::la::Vector<MD_DIM, double> res{24.0 * _epsilon / rij2 * (_sigma6 / rij6) * (1.0 - 2.0 * (_sigma6 / rij6)) * rij};
    constexpr double maxF = 1e6;
    constexpr double stepF = (double)(std::numeric_limits<long long>::max()) / maxF;
    res = stepF * res;
    long long fb0{(long long)(res[0])};
    long long fb1{(long long)(res[1])};
    long long fb2{(long long)(res[2])};
    long long& fb0r = fb0;
    long long& fb1r = fb1;
    long long& fb2r = fb2;
    res[0] = *(double*)(&fb0r);
    res[1] = *(double*)(&fb1r);
    res[2] = *(double*)(&fb2r);
    return res;
#else
    return 24.0 * _epsilon / rij2 * (_sigma6 / rij6) * (1.0 - 2.0 * (_sigma6 / rij6)) * rij;
#endif
  } else {
    return tarch::la::Vector<MD_DIM, double>(0.0);
  }
}
