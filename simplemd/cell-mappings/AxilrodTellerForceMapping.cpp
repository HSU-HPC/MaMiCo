// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#include "simplemd/cell-mappings/AxilrodTellerForceMapping.h"

simplemd::cellmappings::AxilrodTellerForceMapping::AxilrodTellerForceMapping(simplemd::services::ExternalForceService& externalForceService,
                                                                           const simplemd::services::MolecularPropertiesService& molecularPropertiesService)
    : _v(molecularPropertiesService.getMolecularProperties().getV()),
      _cutOffRadiusSquared(molecularPropertiesService.getMolecularProperties().getCutOffRadius() *
                           molecularPropertiesService.getMolecularProperties().getCutOffRadius()),
      _externalForceService(externalForceService) {}

void simplemd::cellmappings::AxilrodTellerForceMapping::beginCellIteration() {}

void simplemd::cellmappings::AxilrodTellerForceMapping::handleCell(const LinkedCell& cell, const unsigned int& cellIndex) const {
  // force buffer
  std::vector<tarch::la::Vector<MD_DIM, double>> forceBuffer(3, tarch::la::Vector<MD_DIM, double> (0.0));

  // iterate over all molecules
  const std::list<Molecule*>::const_iterator end = cell.constEnd();
  const std::list<Molecule*>::const_iterator begin = cell.constBegin();
  for (std::list<Molecule*>::const_iterator m1 = begin; m1 != end; m1++) {
    std::list<Molecule*>::const_iterator m2 = m1;
    tarch::la::Vector<MD_DIM, double>& force1 = (*m1)->getForce();
    const tarch::la::Vector<MD_DIM, double>& position1 = (*m1)->getConstPosition();

    // add external force
    _externalForceService.addExternalForce(force1);

    // iterate over all other molecules not touched so far
    m2++;
    while (m2 != end) {
      tarch::la::Vector<MD_DIM, double>& force2 = (*m2)->getForce();
      const tarch::la::Vector<MD_DIM, double>& position2 = (*m2)->getConstPosition();  

      std::list<Molecule*>::const_iterator m3 = m2;
      m3++;
      while (m3 != end) {  
#if (MD_DEBUG == MD_YES)
        std::cout << "Compute force " << (*m1)->getID() << " <-> " << (*m2)->getID() << " <-> " << (*m3)->getID() << std::endl;
#endif

        tarch::la::Vector<MD_DIM, double>& force3 = (*m3)->getForce();
        forceBuffer = getAxilrodTellerForce(position1, position2, (*m3)->getConstPosition());
        force1 -= forceBuffer[0];
        force2 -= forceBuffer[1];
        force3 -= forceBuffer[2];

        m3++;
      }
      m2++;
    }
  }
}

void simplemd::cellmappings::AxilrodTellerForceMapping::handleCellPair(const LinkedCell& cell1, const LinkedCell& cell2, const unsigned int& cellIndex1,
                                                                      const unsigned int& cellIndex2) const {

  // force buffer
  std::vector<tarch::la::Vector<MD_DIM, double>> forceBuffer(3, tarch::la::Vector<MD_DIM, double> (0.0));

  // iterate over triplets of molecules
  const std::list<Molecule*>::const_iterator endCell1 = cell1.constEnd();
  const std::list<Molecule*>::const_iterator endCell2 = cell2.constEnd();
  const std::list<Molecule*>::const_iterator beginCell1 = cell1.constBegin();
  const std::list<Molecule*>::const_iterator beginCell2 = cell2.constBegin();

  // iterate over triplets of molecules with two molecules from cell1
  for (std::list<Molecule*>::const_iterator m1 = beginCell1; m1 != endCell1; m1++) {
    std::list<Molecule*>::const_iterator m2 = m1;
    tarch::la::Vector<MD_DIM, double>& force1 = (*m1)->getForce();
    const tarch::la::Vector<MD_DIM, double>& position1 = (*m1)->getConstPosition();

    // iterate over all molecules not touched so far from cell1
    m2++;
    while (m2 != endCell1) {
      tarch::la::Vector<MD_DIM, double>& force2 = (*m2)->getForce();
      const tarch::la::Vector<MD_DIM, double>& position2 = (*m2)->getConstPosition();

      for (std::list<Molecule*>::const_iterator m3 = beginCell2; m3 != endCell2; m3++) {
#if (MD_DEBUG == MD_YES)
        std::cout << "Compute force " << (*m1)->getID() << " <-> " << (*m2)->getID() << " <-> " << (*m3)->getID() << std::endl;
#endif
        tarch::la::Vector<MD_DIM, double>& force3 = (*m3)->getForce();

        forceBuffer = getAxilrodTellerForce(position1, position2, (*m3)->getConstPosition());
        force1 -= forceBuffer[0];
        force2 -= forceBuffer[1];
        force3 -= forceBuffer[2];
      }
      m2++;
    }
  }

  // iterate over triplets of molecules with two molecules from cell2
  for (std::list<Molecule*>::const_iterator m1 = beginCell1; m1 != endCell1; m1++) {
    tarch::la::Vector<MD_DIM, double>& force1 = (*m1)->getForce();
    const tarch::la::Vector<MD_DIM, double>& position1 = (*m1)->getConstPosition();

    for (std::list<Molecule*>::const_iterator m2 = beginCell2; m2 != endCell2; m2++) {
      std::list<Molecule*>::const_iterator m3 = m2;
      tarch::la::Vector<MD_DIM, double>& force2 = (*m2)->getForce();
      const tarch::la::Vector<MD_DIM, double>& position2 = (*m2)->getConstPosition();

      // iterate over all molecoles not touched so far from cell2
      m3++;
      while (m3 != endCell2) {  
#if (MD_DEBUG == MD_YES)
        std::cout << "Compute force " << (*m1)->getID() << " <-> " << (*m2)->getID() << " <-> " << (*m3)->getID() << std::endl;
#endif
        tarch::la::Vector<MD_DIM, double>& force3 = (*m3)->getForce();

        forceBuffer = getAxilrodTellerForce(position1, position2, (*m3)->getConstPosition());
        force1 -= forceBuffer[0];
        force2 -= forceBuffer[1];
        force3 -= forceBuffer[2];
        m3++;
      }
    }
  }
}

void simplemd::cellmappings::AxilrodTellerForceMapping::handleCellTriplet(const LinkedCell& cell1, const LinkedCell& cell2, const LinkedCell& cell3,
                                                                      const unsigned int& cellIndex1, const unsigned int& cellIndex2, const unsigned int& cellIndex3) const {

  // force buffer
  std::vector<tarch::la::Vector<MD_DIM, double>> forceBuffer(3, tarch::la::Vector<MD_DIM, double> (0.0));

  // iterate over triplets of molecules
  const std::list<Molecule*>::const_iterator endCell1 = cell1.constEnd();
  const std::list<Molecule*>::const_iterator endCell2 = cell2.constEnd();
  const std::list<Molecule*>::const_iterator endCell3 = cell3.constEnd();
  const std::list<Molecule*>::const_iterator beginCell1 = cell1.constBegin();
  const std::list<Molecule*>::const_iterator beginCell2 = cell2.constBegin();
  const std::list<Molecule*>::const_iterator beginCell3 = cell3.constBegin();
  for (std::list<Molecule*>::const_iterator m1 = beginCell1; m1 != endCell1; m1++) {
    tarch::la::Vector<MD_DIM, double>& force1 = (*m1)->getForce();
    const tarch::la::Vector<MD_DIM, double>& position1 = (*m1)->getConstPosition();

    for (std::list<Molecule*>::const_iterator m2 = beginCell2; m2 != endCell2; m2++) {
      tarch::la::Vector<MD_DIM, double>& force2 = (*m2)->getForce();
      const tarch::la::Vector<MD_DIM, double>& position2 = (*m2)->getConstPosition();  

      for (std::list<Molecule*>::const_iterator m3 = beginCell3; m3 != endCell3; m3++) {
#if (MD_DEBUG == MD_YES)
        std::cout << "Compute force " << (*m1)->getID() << " <-> " << (*m2)->getID() << " <-> " << (*m3)->getID() << std::endl;
#endif
        tarch::la::Vector<MD_DIM, double>& force3 = (*m3)->getForce();

        forceBuffer = getAxilrodTellerForce(position1, position2, (*m3)->getConstPosition());
        force1 -= forceBuffer[0];
        force2 -= forceBuffer[1];
        force3 -= forceBuffer[2];
      }
    }
  }
}

std::vector<tarch::la::Vector<MD_DIM, double>>
simplemd::cellmappings::AxilrodTellerForceMapping::getAxilrodTellerForce(const tarch::la::Vector<MD_DIM, double>& position1,
                                                                       const tarch::la::Vector<MD_DIM, double>& position2,
                                                                       const tarch::la::Vector<MD_DIM, double>& position3) const {
  // Vectors between molecules
  const tarch::la::Vector<MD_DIM, double> rij(position2 - position1);
  const tarch::la::Vector<MD_DIM, double> rik(position3 - position1);
  const tarch::la::Vector<MD_DIM, double> rjk(position3 - position2);

  const double rij2 = tarch::la::dot(rij, rij);
  const double rik2 = tarch::la::dot(rik, rik);
  const double rjk2 = tarch::la::dot(rjk, rjk);

#if (MD_ERROR == MD_YES)
  if (tarch::la::equals(rij2, 0.0, 1e-4)) {
    std::cout << "ERROR simplemd::cellmappings::AxilrodTellerForceMapping::getAxilrodTellerForce(): Particle positions are identical!" << std::endl;
    std::cout << "Position1: " << position1 << ","
              << "Position2: " << position2 << std::endl;
  }
  if (tarch::la::equals(rik2, 0.0, 1e-4)) {
    std::cout << "ERROR simplemd::cellmappings::AxilrodTellerForceMapping::getAxilrodTellerForce(): Particle positions are identical!" << std::endl;
    std::cout << "Position1: " << position1 << ","
              << "Position3: " << position3 << std::endl;
  }
  if (tarch::la::equals(rjk2, 0.0, 1e-4)) {
    std::cout << "ERROR simplemd::cellmappings::AxilrodTellerForceMapping::getAxilrodTellerForce(): Particle positions are identical!" << std::endl;
    std::cout << "Position2: " << position2 << ","
              << "Position3: " << position3 << std::endl;
  }
#endif

  if (rij2 <= _cutOffRadiusSquared && rik2 <= _cutOffRadiusSquared && rjk2 <= _cutOffRadiusSquared) {

    // Distances between molecules and powers of those
    const double rij1 = tarch::la::norm2(rij);
    const double rik1 = tarch::la::norm2(rik);
    const double rjk1 = tarch::la::norm2(rjk);

    const double rij3 = rij1 * rij2;
    const double rik3 = rik1 * rik2;
    const double rjk3 = rjk1 * rjk2;

    const double rij4 = rij2 * rij2;
    const double rik4 = rik2 * rik2;
    const double rjk4 = rjk2 * rjk2;

    const double rij5 = rij4 * rij1;
    const double rik5 = rik4 * rik1;
    const double rjk5 = rjk4 * rjk1;

    const double rij6 = rij4 * rij2;
    const double rik6 = rik4 * rik2;
    const double rjk6 = rjk4 * rjk2;

    // Gradient
    const double dVdRij = ((3.0 * _v / (8.0 * rij1)) * (-8.0 / (rij4 * rik3 * rjk3) - 1.0 / (rik5 * rjk5)
                        + 5.0 * rik1 / (rij6 * rjk5) + 5.0 * rjk1 / (rij6 * rik5) - 1.0 / (rij2 * rik3 * rjk5)
                        - 1.0 / (rij2 * rik5 * rjk3) - 3.0 / (rij4 * rik1 * rjk5) - 3.0 / (rij4 * rik5 * rjk1)
                        - 5.0 / (rij6 * rik1 * rjk3) - 5.0 / (rij6 * rik3 * rjk1) + 6.0 / (rij4 * rik3 * rjk3)));
    const double dVdRik = ((3.0 * _v / (8.0 * rik1)) * (-8.0 / (rik4 * rij3 * rjk3) - 1.0 / (rij5 * rjk5)
                        + 5.0 * rij1 / (rik6 * rjk5) + 5.0 * rjk1 / (rik6 * rij5) - 1.0 / (rik2 * rij3 * rjk5)
                        - 1.0 / (rik2 * rij5 * rjk3) - 3.0 / (rik4 * rij1 * rjk5) - 3.0 / (rik4 * rij5 * rjk1)
                        - 5.0 / (rik6 * rij1 * rjk3) - 5.0 / (rik6 * rij3 * rjk1) + 6.0 / (rik4 * rij3 * rjk3)));
    const double dVdRjk = ((3.0 * _v / (8.0 * rjk1)) * (-8.0 / (rjk4 * rik3 * rij3) - 1.0 / (rik5 * rij5)
                        + 5.0 * rik1 / (rjk6 * rij5) + 5.0 * rij1 / (rjk6 * rik5) - 1.0 / (rjk2 * rik3 * rij5)
                        - 1.0 / (rjk2 * rik5 * rij3) - 3.0 / (rjk4 * rik1 * rij5) - 3.0 / (rjk4 * rik5 * rij1)
                        - 5.0 / (rjk6 * rik1 * rij3) - 5.0 / (rjk6 * rik3 * rij1) + 6.0 / (rjk4 * rik3 * rij3)));

    // Forces
    const tarch::la::Vector<MD_DIM, double> Fi = rij * dVdRij + rik * dVdRik;
    const tarch::la::Vector<MD_DIM, double> Fj = rij * (-dVdRij) + rjk * dVdRjk;
    const tarch::la::Vector<MD_DIM, double> Fk = rik * (-dVdRik) + rjk * (-dVdRjk);

    std::vector<tarch::la::Vector<MD_DIM, double>> forces = {Fi, Fj, Fk};

    return forces;
  } else {
    std::vector<tarch::la::Vector<MD_DIM, double>> forces(3, tarch::la::Vector<MD_DIM, double>(0.0));
    return forces;
  }
}