// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#include "simplemd/cell-mappings/AxilrodTellerPotentialEnergyMapping.h"

simplemd::cellmappings::AxilrodTellerPotentialEnergyMapping::AxilrodTellerPotentialEnergyMapping(
    const simplemd::services::MolecularPropertiesService& molecularPropertiesService)
    : _v(molecularPropertiesService.getMolecularProperties().getV()),
      _cutOffRadiusSquared(molecularPropertiesService.getMolecularProperties().getCutOffRadius() *
                           molecularPropertiesService.getMolecularProperties().getCutOffRadius()),
      _cutOffEnergy(_v * 1.375 / (_cutOffRadiusSquared * _cutOffRadiusSquared * _cutOffRadiusSquared * // TODO
                    _cutOffRadiusSquared * molecularPropertiesService.getMolecularProperties().getCutOffRadius())) {}

void simplemd::cellmappings::AxilrodTellerPotentialEnergyMapping::beginCellIteration() {}

void simplemd::cellmappings::AxilrodTellerPotentialEnergyMapping::handleCell(LinkedCell& cell, const unsigned int& cellIndex) {
  
  // potential energy buffer
  tarch::la::Vector<3, double> potentialEnergyBuffer(0.0);

  // iterate over all molecules
  const std::list<Molecule*>::const_iterator itEnd = cell.end();
  const std::list<Molecule*>::const_iterator itBegin = cell.begin();
  for (std::list<Molecule*>::const_iterator m1 = itBegin; m1 != itEnd; m1++) {
    std::list<Molecule*>::const_iterator m2 = m1;
#if (AD_RES == MD_NO)
    double& potentialEnergy1 = (*m1)->getPotentialEnergy();
#else
    double& potentialEnergy1 = (*m1)->getThreeBodyPotentialEnergy();
#endif

    // iterate over all other molecules not touched so far
    m2++;
    while (m2 != itEnd) {
      std::list<Molecule*>::const_iterator m3 = m2;
#if (AD_RES == MD_NO)
      double& potentialEnergy2 = (*m2)->getPotentialEnergy();
#else
      double& potentialEnergy2 = (*m2)->getThreeBodyPotentialEnergy();
#endif

      while (m3 != itEnd) {
#if (MD_DEBUG == MD_YES)
        std::cout << "Compute potential energy " << (*m1)->getID() << " <-> " << (*m2)->getID() << " <-> " << (*m3)->getID() << std::endl;
#endif

#if (AD_RES == MD_NO)
        double& potentialEnergy3 = (*m3)->getPotentialEnergy();
#else
        double& potentialEnergy3 = (*m3)->getThreeBodyPotentialEnergy();
#endif
        potentialEnergyBuffer = getAxilrodTellerPotentialEnergy((*m1)->getConstPosition(), (*m2)->getConstPosition(), (*m3)->getConstPosition());

        potentialEnergy1 += potentialEnergyBuffer[0];
        potentialEnergy2 += potentialEnergyBuffer[1];
        potentialEnergy3 += potentialEnergyBuffer[2];
        m3++;
      }

      m2++;
    }
  }
}

void simplemd::cellmappings::AxilrodTellerPotentialEnergyMapping::handleCellPair(LinkedCell& cell1, LinkedCell& cell2, const unsigned int& cellIndex1,
                                                                                const unsigned int& cellIndex2) {
  
  // potential energy buffer
  tarch::la::Vector<3, double> potentialEnergyBuffer(0.0);

  // iterate over triplets of molecules
  const std::list<Molecule*>::const_iterator m1End = cell1.end();
  const std::list<Molecule*>::const_iterator m2End = cell2.end();
  const std::list<Molecule*>::const_iterator m1Begin = cell1.begin();
  const std::list<Molecule*>::const_iterator m2Begin = cell2.begin();
  // iterate over triplets of molecules with two molecules from cell1
  for (std::list<Molecule*>::const_iterator m1 = m1Begin; m1 != m1End; m1++) {
    std::list<Molecule*>::const_iterator m2 = m1;
#if (AD_RES == MD_NO)
    double& potentialEnergy1 = (*m1)->getPotentialEnergy();
#else
    double& potentialEnergy1 = (*m1)->getThreeBodyPotentialEnergy();
#endif

    // iterate over all molecules not touched so far from cell1
    m2++;
    while (m2 != m1End) {
#if (AD_RES == MD_NO)
      double& potentialEnergy2 = (*m2)->getPotentialEnergy();
#else
      double& potentialEnergy2 = (*m2)->getThreeBodyPotentialEnergy();
#endif

      for (std::list<Molecule*>::const_iterator m3 = m2Begin; m3 != m2End; m3++) {
#if (MD_DEBUG == MD_YES)
        std::cout << "Compute potential energy " << (*m1)->getID() << " <-> " << (*m2)->getID() << " <-> " << (*m3)->getID() << std::endl;
#endif
#if (AD_RES == MD_NO)
        double& potentialEnergy3 = (*m3)->getPotentialEnergy();
#else
        double& potentialEnergy3 = (*m3)->getThreeBodyPotentialEnergy();
#endif

        potentialEnergyBuffer = getAxilrodTellerPotentialEnergy((*m1)->getConstPosition(), (*m2)->getConstPosition(), (*m3)->getConstPosition());

        potentialEnergy1 += potentialEnergyBuffer[0];
        potentialEnergy2 += potentialEnergyBuffer[1];
        potentialEnergy3 += potentialEnergyBuffer[2];
      }
      m2++;
    }
  }

  // iterate over triplets of molecules with two molecules from cell2
  for (std::list<Molecule*>::const_iterator m1 = m1Begin; m1 != m1End; m1++) {
#if (AD_RES == MD_NO)
    double& potentialEnergy1 = (*m1)->getPotentialEnergy();
#else
    double& potentialEnergy1 = (*m1)->getThreeBodyPotentialEnergy();
#endif

    for (std::list<Molecule*>::const_iterator m2 = m2Begin; m2 != m2End; m2++) {
      std::list<Molecule*>::const_iterator m3 = m2;
#if (AD_RES == MD_NO)
      double& potentialEnergy2 = (*m2)->getPotentialEnergy();
#else
      double& potentialEnergy2 = (*m2)->getThreeBodyPotentialEnergy();
#endif

      // iterate over all molecoles not touched so far from cell2
      m3++;
      while (m3 != m2End) {  
#if (MD_DEBUG == MD_YES)
        std::cout << "Compute potential energy " << (*m1)->getID() << " <-> " << (*m2)->getID() << " <-> " << (*m3)->getID() << std::endl;
#endif
#if (AD_RES == MD_NO)
        double& potentialEnergy3 = (*m3)->getPotentialEnergy();
#else
        double& potentialEnergy3 = (*m3)->getThreeBodyPotentialEnergy();
#endif

        potentialEnergyBuffer = getAxilrodTellerPotentialEnergy((*m1)->getConstPosition(), (*m2)->getConstPosition(), (*m3)->getConstPosition());

        potentialEnergy1 += potentialEnergyBuffer[0];
        potentialEnergy2 += potentialEnergyBuffer[1];
        potentialEnergy3 += potentialEnergyBuffer[2];
        m3++;
      }
    }
  }
}

void simplemd::cellmappings::AxilrodTellerPotentialEnergyMapping::handleCellTriplet(LinkedCell& cell1, LinkedCell& cell2, LinkedCell& cell3,
                                                                      const unsigned int& cellIndex1, const unsigned int& cellIndex2, const unsigned int& cellIndex3) {

  // potential energy buffer
  tarch::la::Vector<3, double> potentialEnergyBuffer(0.0);

  // iterate over triplets of molecules
  const std::list<Molecule*>::const_iterator m1End = cell1.constEnd();
  const std::list<Molecule*>::const_iterator m2End = cell2.constEnd();
  const std::list<Molecule*>::const_iterator m3End = cell3.constEnd();
  const std::list<Molecule*>::const_iterator m1Begin = cell1.constBegin();
  const std::list<Molecule*>::const_iterator m2Begin = cell2.constBegin();
  const std::list<Molecule*>::const_iterator m3Begin = cell3.constBegin();
  for (std::list<Molecule*>::const_iterator m1 = m1Begin; m1 != m1End; m1++) {
#if (AD_RES == MD_NO)
    double& potentialEnergy1 = (*m1)->getPotentialEnergy();
#else
    double& potentialEnergy1 = (*m1)->getThreeBodyPotentialEnergy();
#endif

    for (std::list<Molecule*>::const_iterator m2 = m2Begin; m2 != m2End; m2++) {
#if (AD_RES == MD_NO)
      double& potentialEnergy2 = (*m2)->getPotentialEnergy();
#else
      double& potentialEnergy2 = (*m2)->getThreeBodyPotentialEnergy();
#endif

      for (std::list<Molecule*>::const_iterator m3 = m3Begin; m3 != m3End; m3++) {
#if (MD_DEBUG == MD_YES)
        std::cout << "Compute potential energy " << (*m1)->getID() << " <-> " << (*m2)->getID() << " <-> " << (*m3)->getID() << std::endl;
#endif
#if (AD_RES == MD_NO)
        double& potentialEnergy3 = (*m3)->getPotentialEnergy();
#else
        double& potentialEnergy3 = (*m3)->getThreeBodyPotentialEnergy();
#endif

        potentialEnergyBuffer = getAxilrodTellerPotentialEnergy((*m1)->getConstPosition(), (*m2)->getConstPosition(), (*m3)->getConstPosition());

        potentialEnergy1 += potentialEnergyBuffer[0];
        potentialEnergy2 += potentialEnergyBuffer[1];
        potentialEnergy3 += potentialEnergyBuffer[2];
      }
    }
  }
}

tarch::la::Vector<3, double>
simplemd::cellmappings::AxilrodTellerPotentialEnergyMapping::getAxilrodTellerPotentialEnergy(const tarch::la::Vector<MD_DIM, double>& position1,
                                                                       const tarch::la::Vector<MD_DIM, double>& position2,
                                                                       const tarch::la::Vector<MD_DIM, double>& position3) {

  const double rij2 = tarch::la::dot(position2 - position1, position2 - position1);
  const double rik2 = tarch::la::dot(position3 - position1, position3 - position1);
  const double rjk2 = tarch::la::dot(position3 - position2, position3 - position2);

#if (MD_ERROR == MD_YES)
  if (tarch::la::equals(rij2, 0.0, 1e-4)) {
    std::cout << "ERROR simplemd::cellmappings::AxilrodTellerForceMapping::getAxilrodTellerForce(): Particle positions are identical!" << std::endl;
    std::cout << "Position1: " << position1 << ", "
              << "Position2: " << position2 << std::endl;
  }
  if (tarch::la::equals(rik2, 0.0, 1e-4)) {
    std::cout << "ERROR simplemd::cellmappings::AxilrodTellerForceMapping::getAxilrodTellerForce(): Particle positions are identical!" << std::endl;
    std::cout << "Position1: " << position1 << ", "
              << "Position3: " << position3 << std::endl;
  }
  if (tarch::la::equals(rjk2, 0.0, 1e-4)) {
    std::cout << "ERROR simplemd::cellmappings::AxilrodTellerForceMapping::getAxilrodTellerForce(): Particle positions are identical!" << std::endl;
    std::cout << "Position2: " << position2 << ", "
              << "Position3: " << position3 << std::endl;
  }
#endif

  if (rij2 <= _cutOffRadiusSquared && rik2 <= _cutOffRadiusSquared && rjk2 <= _cutOffRadiusSquared) {

    // Distances between molecules and powers of those
    const double rij1 = std::sqrt(rij2);
    const double rik1 = std::sqrt(rik2);
    const double rjk1 = std::sqrt(rjk2);

    const double rij3 = rij1 * rij2;
    const double rik3 = rik1 * rik2;
    const double rjk3 = rjk1 * rjk2;

    const double rij5 = rij3 * rij2;
    const double rik5 = rik3 * rik2;
    const double rjk5 = rjk3 * rjk2;

    // Potential Energy
    const double totalEnergy = _v * (1.0 / (rij3 * rik3 * rjk3) 
                    + (3.0 * (-rij2 + rik2 + rjk2) * (rij2 - rik2 + rjk2) * (rij2 + rik2 - rjk2)) / (8.0 * rij5 * rik5 * rjk5));

    // Factors
    // For now all factors are equal (potential energies are never used individually)
    // Can be changed later
    const double fi = 1.0 / 3.0;
    const double fj = 1.0 / 3.0;
    const double fk = 1.0 / 3.0;

    tarch::la::Vector<3, double> factors;
    factors[0] = fi;
    factors[1] = fj;
    factors[2] = fk;

    return factors * totalEnergy;
  } else {
    return tarch::la::Vector<3, double>(0.0);
  }
}