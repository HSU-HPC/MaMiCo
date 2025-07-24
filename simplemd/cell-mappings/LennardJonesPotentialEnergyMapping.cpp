// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#include "simplemd/cell-mappings/LennardJonesPotentialEnergyMapping.h"

simplemd::cellmappings::LennardJonesPotentialEnergyMapping::LennardJonesPotentialEnergyMapping(
    const simplemd::services::MolecularPropertiesService& molecularPropertiesService)
    : _epsilon(molecularPropertiesService.getMolecularProperties().getEpsilon()),
      _sigma6(molecularPropertiesService.getMolecularProperties().getSigma() * molecularPropertiesService.getMolecularProperties().getSigma() *
              molecularPropertiesService.getMolecularProperties().getSigma() * molecularPropertiesService.getMolecularProperties().getSigma() *
              molecularPropertiesService.getMolecularProperties().getSigma() * molecularPropertiesService.getMolecularProperties().getSigma()),
      _cutOffRadiusSquared(molecularPropertiesService.getMolecularProperties().getCutOffRadius() *
                           molecularPropertiesService.getMolecularProperties().getCutOffRadius()),
      _cutOffEnergy(4.0 * _epsilon * _sigma6 / (_cutOffRadiusSquared * _cutOffRadiusSquared * _cutOffRadiusSquared) *
                    (_sigma6 / (_cutOffRadiusSquared * _cutOffRadiusSquared * _cutOffRadiusSquared) - 1.0)) {}

void simplemd::cellmappings::LennardJonesPotentialEnergyMapping::beginCellIteration() {}

void simplemd::cellmappings::LennardJonesPotentialEnergyMapping::handleCell(LinkedCell& cell, const unsigned int& cellIndex) {

  simplemd::services::MoleculeService* _moleculeService = nullptr;
  throw "Not yet implemented: init _moleculeService";

  // iterate over all molecules
  auto itEnd = cell.end();
  auto itBegin = cell.begin(*_moleculeService);
  for (auto m1 = itBegin; m1 != itEnd; ++m1) {
    auto m2 = m1;
    double& potentialEnergy1 = (*m1)->getPotentialEnergy();

    // iterate over all other molecules not touched so far
    ++m2;
    while (m2 != itEnd) {
#if (MD_DEBUG == MD_YES)
      std::cout << "Compute potential energy " << (*m1)->getID() << " <-> " << (*m2)->getID() << std::endl;
#endif

      double& potentialEnergy2 = (*m2)->getPotentialEnergy();
      const double rij2 = tarch::la::dot(((*m2)->getConstPosition() - (*m1)->getConstPosition()), ((*m2)->getConstPosition() - (*m1)->getConstPosition()));
#if (MD_ERROR == MD_YES)
      if (rij2 == 0.0) {
        std::cout << "ERROR simplemd::cellmappings::LennardJonesPotentialEnergyMapping::handleCell(): Particle positions are identical!" << std::endl;
        exit(EXIT_FAILURE);
      }
#endif

      if (rij2 < _cutOffRadiusSquared) {
        const double rij6 = rij2 * rij2 * rij2;
        const double energyBuffer = 4.0 * _epsilon * (_sigma6 / rij6) * ((_sigma6 / rij6) - 1.0) - _cutOffEnergy;
        potentialEnergy1 += 0.5 * energyBuffer;
        potentialEnergy2 += 0.5 * energyBuffer;
      }

      ++m2;
    }
  }
}

void simplemd::cellmappings::LennardJonesPotentialEnergyMapping::handleCellPair(LinkedCell& cell1, LinkedCell& cell2, const unsigned int& cellIndex1,
                                                                                const unsigned int& cellIndex2) {

  simplemd::services::MoleculeService* _moleculeService = nullptr;
  throw "Not yet implemented: init _moleculeService";

  // iterate over pairs of molecules
  auto m1End = cell1.end();
  auto m2End = cell2.end();
  auto m1Begin = cell1.begin(*_moleculeService);
  auto m2Begin = cell2.begin(*_moleculeService);
  for (auto m1 = m1Begin; m1 != m1End; ++m1) {
    double& potentialEnergy1 = (*m1)->getPotentialEnergy();

    for (auto m2 = m2Begin; m2 != m2End; ++m2) {
#if (MD_DEBUG == MD_YES)
      std::cout << "Compute potential energy " << (*m1)->getID() << " <-> " << (*m2)->getID() << std::endl;
#endif

      double& potentialEnergy2 = (*m2)->getPotentialEnergy();
      const double rij2 = tarch::la::dot(((*m2)->getConstPosition() - (*m1)->getConstPosition()), ((*m2)->getConstPosition() - (*m1)->getConstPosition()));
#if (MD_ERROR == MD_YES)
      if (rij2 == 0.0) {
        std::cout << "ERROR simplemd::cellmappings::LennardJonesPotentialEnergyMapping::handleCellPairs(): Particle positions are identical!" << std::endl;
        exit(EXIT_FAILURE);
      }
#endif

      if (rij2 < _cutOffRadiusSquared) {
        const double rij6 = rij2 * rij2 * rij2;
        const double energyBuffer = 4.0 * _epsilon * (_sigma6 / rij6) * ((_sigma6 / rij6) - 1.0) - _cutOffEnergy;
        potentialEnergy1 += 0.5 * energyBuffer;
        potentialEnergy2 += 0.5 * energyBuffer;
      }
    }
  }
}
