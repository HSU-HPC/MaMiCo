// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#include "simplemd/cell-mappings/LennardJonesForceMapping.h"

simplemd::cellmappings::LennardJonesForceMapping::LennardJonesForceMapping(simplemd::services::ExternalForceService& externalForceService,
                                                                           const simplemd::services::MolecularPropertiesService& molecularPropertiesService)
    : _epsilon(molecularPropertiesService.getMolecularProperties().getEpsilon()),
      _sigma6(molecularPropertiesService.getMolecularProperties().getSigma() * molecularPropertiesService.getMolecularProperties().getSigma() *
              molecularPropertiesService.getMolecularProperties().getSigma() * molecularPropertiesService.getMolecularProperties().getSigma() *
              molecularPropertiesService.getMolecularProperties().getSigma() * molecularPropertiesService.getMolecularProperties().getSigma()),
      _cutOffRadiusSquared(molecularPropertiesService.getMolecularProperties().getCutOffRadius() *
                           molecularPropertiesService.getMolecularProperties().getCutOffRadius()),
      _externalForceService(externalForceService) {}

void simplemd::cellmappings::LennardJonesForceMapping::beginCellIteration() {
#if (MD_DEBUG == MD_YES)
  std::cout << "simplemd::cellmappings::LennardJonesForceMapping::beginCellIteration() " << std::endl;
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

    // add external force
    _externalForceService.addExternalForce(force1);

    // iterate over all other molecules not touched so far
    m2++;
    while (m2 != end) {
      tarch::la::Vector<MD_DIM, double>& force2 = m2->getForce();
      forceBuffer = getLennardJonesForce(position1, m2->getConstPosition());
#if (MD_DEBUG == MD_YES)
      if (tarch::la::dot(forceBuffer, forceBuffer) > 1e12) {
        std::cout << "ERROR simplemd::cellmappings::LennardJonesForceMapping::handleCell: Force " << forceBuffer << " out of range!" << std::endl;
        std::cout << "Position1: " << position1 << std::endl;
        std::cout << "Position2: " << m2->getConstPosition() << std::endl;
        std::cout << "ID1: " << m1->getID() << std::endl;
        std::cout << "ID2: " << m2->getID() << std::endl;
        exit(EXIT_FAILURE);
      }
#endif
      force1 += forceBuffer;
      force2 -= forceBuffer;

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
        std::cout << "ERROR simplemd::cellmappings::LennardJonesForceMapping::handleCellPair: Force " << forceBuffer << " out of range!" << std::endl;
        std::cout << "Position1: " << position1 << std::endl;
        std::cout << "Position2: " << m2->getConstPosition() << std::endl;
        std::cout << "ID1: " << m1->getID() << std::endl;
        std::cout << "ID2: " << m2->getID() << std::endl;
        exit(EXIT_FAILURE);
      }
#endif
      force1 += forceBuffer;
      force2 -= forceBuffer;
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
    Kokkos::printf("Position:");
    for (int d = 0; d < MD_DIM; d++)
      Kokkos::printf(" %f", position1[d]);
    Kokkos::printf(",");
    for (int d = 0; d < MD_DIM; d++)
      Kokkos::printf(" %f", position2[d]);
    Kokkos::printf("\n");
    Kokkos::abort("ERROR simplemd::cellmappings::LennardJonesForceMapping::getLennardJonesForce(): Particle positions are identical!");
  }
#endif

  if (rij2 <= _cutOffRadiusSquared) {
    const double rij6 = rij2 * rij2 * rij2;
    return 24.0 * _epsilon / rij2 * (_sigma6 / rij6) * (1.0 - 2.0 * (_sigma6 / rij6)) * rij;
  } else {
    return tarch::la::Vector<MD_DIM, double>(0.0);
  }
}
