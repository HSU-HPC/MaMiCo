// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#include "simplemd/molecule-mappings/VelocityStoermerVerletMapping.h"

simplemd::moleculemappings::VelocityStoermerVerletMapping::VelocityStoermerVerletMapping(
const double& kB,const double& dt,const double& mass,
const tarch::la::Vector<MD_LINKED_CELL_NEIGHBOURS,simplemd::BoundaryType>& boundary,
const tarch::la::Vector<MD_DIM,double> &domainOffset, const tarch::la::Vector<MD_DIM,double> &domainSize
):
_dt(dt),
_a(_dt/(2.0*mass)),_zero(0.0),
_boundary(initReflectingBoundary(boundary)),
_domainOffset(domainOffset),
_domainSize(domainSize)
{}

simplemd::moleculemappings::VelocityStoermerVerletMapping::~VelocityStoermerVerletMapping(){}


void simplemd::moleculemappings::VelocityStoermerVerletMapping::beginMoleculeIteration(){}


void simplemd::moleculemappings::VelocityStoermerVerletMapping::endMoleculeIteration(){}



void simplemd::moleculemappings::VelocityStoermerVerletMapping::handleMolecule(Molecule &molecule){
  // if the molecule is fixed in space, return immediately:
  if(molecule.isFixed()) return;


  // do time integration update (position and velocity) ------------------------------
  tarch::la::Vector<MD_DIM,double> &position = molecule.getPosition();
  const tarch::la::Vector<MD_DIM,double> oldPosition(molecule.getConstPosition());
  #if (MD_ERROR==MD_YES)
  const tarch::la::Vector<MD_DIM,double> oldVelocity(molecule.getConstVelocity());
  #endif
  tarch::la::Vector<MD_DIM,double> &velocity = molecule.getVelocity();

  // v_n = v_(n-1) + a*(f_n + f_(n-1))
  velocity += _a*(molecule.getConstForce() + molecule.getConstForceOld());
  // x_(n+1) = x_n + dt*(v_n + a*f_n)
  position += _dt*( molecule.getConstVelocity() + _a*molecule.getConstForce() );
  #if (MD_ERROR == MD_YES)
  for (unsigned int d = 0; d < MD_DIM; d++){
    if ( isnan(position[d]) || isinf(position[d]) ){
      std::cout << "ERROR simplemd::moleculemappings::VelocityStoermerVerletMapping::handleMolecule: Position ";
      std::cout << d << " is out of range" << std::endl;
      std::cout << "Position: " << position << ", molecule: " << molecule.getID() << std::endl;
      std::cout << "Velocity: " << velocity << ", molecule: " << molecule.getID() << std::endl;
      std::cout << "OldVelocity: " << oldVelocity << ", molecule: " << molecule.getID() << std::endl;
      std::cout << "Force: " << molecule.getConstForce() << ", old: " << molecule.getConstForceOld() << std::endl;
      std::cout << "Old position: " << oldPosition << std::endl;
      exit(EXIT_FAILURE);
    }
    if (isnan(velocity[d]) || isinf(velocity[d])){
      std::cout << "ERROR simplemd::moleculemappings::VelocityStoermerVerletMapping::handleMolecule: Velocity ";
      std::cout << d << " is NaN or Inf" << std::endl;
      std::cout << velocity << std::endl;
      std::cout << molecule.getConstForce() << ", " << molecule.getConstForceOld() << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  #endif

  // apply reflection
  for (unsigned int d = 0; d < MD_DIM; d++){
    // left/front/bottom reflecting boundary
    if (_boundary[2*d] && (position[d]<_domainOffset[d]) ){
      //std::cout << "Reflect particle " << position << " d=" << d << ", " << 2*d << std::endl;
      position[d] = position[d] + 2.0*(_domainOffset[d]-position[d]);
      velocity[d] = -velocity[d];
    }
    // right/back/top reflecting boundary
    if (_boundary[2*d+1] && (position[d]>_domainOffset[d]+_domainSize[d]) ){
      //std::cout << "Reflect particle " << position << " d=" << d << ", " << 2*d+1 << std::endl;
      position[d] = position[d] + 2.0*(_domainOffset[d]+_domainSize[d]-position[d]);
      velocity[d] = -velocity[d];
    }
  }

  // store force in force_old
  molecule.setForceOld(molecule.getConstForce());
  // reset force
  molecule.setForce(_zero);
}


tarch::la::Vector<2*MD_DIM,bool> simplemd::moleculemappings::VelocityStoermerVerletMapping::
initReflectingBoundary(const tarch::la::Vector<MD_LINKED_CELL_NEIGHBOURS,simplemd::BoundaryType>& boundary) const{
  tarch::la::Vector<2*MD_DIM,bool> reflect;
  for (unsigned int d = 0; d < 2*MD_DIM; d++){reflect[d] = false;}

  // check 6 sides for reflecting boundaries
  #if (MD_DIM==1)
  if (boundary[0] == simplemd::REFLECTING_BOUNDARY){ reflect[0] = true; }
  if (boundary[1] == simplemd::REFLECTING_BOUNDARY){ reflect[1] = true; }
  #elif (MD_DIM==2)
  // bottom
  if (boundary[1] == simplemd::REFLECTING_BOUNDARY){
    reflect[2] = true;
    if ( (boundary[0] != simplemd::REFLECTING_BOUNDARY) || (boundary[2] != simplemd::REFLECTING_BOUNDARY) ){
      std::cout << "ERROR simplemd::moleculemappings::VelocityStoermerVerletMapping::initReflectingBoundary: Boundaries 0,2 are not reflecting!" << std::endl; exit(EXIT_FAILURE);
    }
  }
  // left
  if (boundary[3] == simplemd::REFLECTING_BOUNDARY){
    reflect[0] = true;
    if ( (boundary[0] != simplemd::REFLECTING_BOUNDARY) || (boundary[5] != simplemd::REFLECTING_BOUNDARY) ){
      std::cout << "ERROR simplemd::moleculemappings::VelocityStoermerVerletMapping::initReflectingBoundary: Boundaries 0,5 are not reflecting!" << std::endl; exit(EXIT_FAILURE);
    }
  }
  // right
  if (boundary[4] == simplemd::REFLECTING_BOUNDARY){
    reflect[1] = true;
    if ( (boundary[2] != simplemd::REFLECTING_BOUNDARY) || (boundary[7] != simplemd::REFLECTING_BOUNDARY) ){
      std::cout << "ERROR simplemd::moleculemappings::VelocityStoermerVerletMapping::initReflectingBoundary: Boundaries 2,7 are not reflecting!" << std::endl; exit(EXIT_FAILURE);
    }
  }
  // top
  if (boundary[6] == simplemd::REFLECTING_BOUNDARY){
    reflect[3] = true;
    if ( (boundary[5] != simplemd::REFLECTING_BOUNDARY) || (boundary[7] != simplemd::REFLECTING_BOUNDARY) ){
      std::cout << "ERROR simplemd::moleculemappings::VelocityStoermerVerletMapping::initReflectingBoundary: Boundaries 5,7 are not reflecting!" << std::endl; exit(EXIT_FAILURE);
    }
  }
  #elif (MD_DIM==3)
  // bottom
  if (boundary[4] == simplemd::REFLECTING_BOUNDARY){
    reflect[4] = true;
    if ( (boundary[0] != simplemd::REFLECTING_BOUNDARY) || (boundary[1] != simplemd::REFLECTING_BOUNDARY) || (boundary[2] != simplemd::REFLECTING_BOUNDARY)
       ||(boundary[3] != simplemd::REFLECTING_BOUNDARY) || (boundary[5] != simplemd::REFLECTING_BOUNDARY) || (boundary[6] != simplemd::REFLECTING_BOUNDARY)
       ||(boundary[7] != simplemd::REFLECTING_BOUNDARY) || (boundary[8] != simplemd::REFLECTING_BOUNDARY) ){
      std::cout << "ERROR simplemd::moleculemappings::VelocityStoermerVerletMapping::initReflectingBoundary: Boundaries 0,1,2,3,5,6,7,8 are not reflecting!" << std::endl; exit(EXIT_FAILURE);
    }
  }
  // front
  if (boundary[10] == simplemd::REFLECTING_BOUNDARY){
    reflect[2] = true;
    if ( (boundary[0] != simplemd::REFLECTING_BOUNDARY) || (boundary[1] != simplemd::REFLECTING_BOUNDARY) || (boundary[2] != simplemd::REFLECTING_BOUNDARY)
       ||(boundary[9] != simplemd::REFLECTING_BOUNDARY) || (boundary[11] != simplemd::REFLECTING_BOUNDARY) || (boundary[17] != simplemd::REFLECTING_BOUNDARY)
       ||(boundary[18] != simplemd::REFLECTING_BOUNDARY) || (boundary[19] != simplemd::REFLECTING_BOUNDARY) ){
      std::cout << "ERROR simplemd::moleculemappings::VelocityStoermerVerletMapping::initReflectingBoundary: Boundaries 0,1,2,9,11,17,18,19 are not reflecting!" << std::endl; exit(EXIT_FAILURE);
    }
  }
  // left
  if (boundary[12] == simplemd::REFLECTING_BOUNDARY){
    reflect[0] = true;
    if ( (boundary[0] != simplemd::REFLECTING_BOUNDARY) || (boundary[3] != simplemd::REFLECTING_BOUNDARY) || (boundary[6] != simplemd::REFLECTING_BOUNDARY)
       ||(boundary[9] != simplemd::REFLECTING_BOUNDARY) || (boundary[14] != simplemd::REFLECTING_BOUNDARY) || (boundary[17] != simplemd::REFLECTING_BOUNDARY)
       ||(boundary[20] != simplemd::REFLECTING_BOUNDARY) || (boundary[23] != simplemd::REFLECTING_BOUNDARY) ){
      std::cout << "ERROR simplemd::moleculemappings::VelocityStoermerVerletMapping::initReflectingBoundary: Boundaries 0,3,6,9,14,17,20,23 are not reflecting!" << std::endl; exit(EXIT_FAILURE);
    }
  }
  // right
  if (boundary[13] == simplemd::REFLECTING_BOUNDARY){
    reflect[1] = true;
    if ( (boundary[2] != simplemd::REFLECTING_BOUNDARY) || (boundary[5] != simplemd::REFLECTING_BOUNDARY) || (boundary[8] != simplemd::REFLECTING_BOUNDARY)
       ||(boundary[11] != simplemd::REFLECTING_BOUNDARY) || (boundary[16] != simplemd::REFLECTING_BOUNDARY) || (boundary[19] != simplemd::REFLECTING_BOUNDARY)
       ||(boundary[22] != simplemd::REFLECTING_BOUNDARY) || (boundary[25] != simplemd::REFLECTING_BOUNDARY) ){
      std::cout << "ERROR simplemd::moleculemappings::VelocityStoermerVerletMapping::initReflectingBoundary: Boundaries 2,5,8,11,16,19,22,25 are not reflecting!" << std::endl; exit(EXIT_FAILURE);
    }
  }
  // back
  if (boundary[15] == simplemd::REFLECTING_BOUNDARY){
    reflect[3] = true;
    if ( (boundary[6] != simplemd::REFLECTING_BOUNDARY) || (boundary[7] != simplemd::REFLECTING_BOUNDARY) || (boundary[8] != simplemd::REFLECTING_BOUNDARY)
       ||(boundary[14] != simplemd::REFLECTING_BOUNDARY) || (boundary[16] != simplemd::REFLECTING_BOUNDARY) || (boundary[23] != simplemd::REFLECTING_BOUNDARY)
       ||(boundary[24] != simplemd::REFLECTING_BOUNDARY) || (boundary[25] != simplemd::REFLECTING_BOUNDARY) ){
      std::cout << "ERROR simplemd::moleculemappings::VelocityStoermerVerletMapping::initReflectingBoundary: Boundaries 6,7,8,14,16,23,24,25 are not reflecting!" << std::endl; exit(EXIT_FAILURE);
    }
  }
  // top
  if (boundary[21] == simplemd::REFLECTING_BOUNDARY){
    reflect[5] = true;
    if ( (boundary[17] != simplemd::REFLECTING_BOUNDARY) || (boundary[18] != simplemd::REFLECTING_BOUNDARY) || (boundary[19] != simplemd::REFLECTING_BOUNDARY)
       ||(boundary[20] != simplemd::REFLECTING_BOUNDARY) || (boundary[22] != simplemd::REFLECTING_BOUNDARY) || (boundary[23] != simplemd::REFLECTING_BOUNDARY)
       ||(boundary[24] != simplemd::REFLECTING_BOUNDARY) || (boundary[25] != simplemd::REFLECTING_BOUNDARY) ){
      std::cout << "ERROR simplemd::moleculemappings::VelocityStoermerVerletMapping::initReflectingBoundary: Boundaries 17,18,19,20,22,23,24,25 are not reflecting!" << std::endl; exit(EXIT_FAILURE);
    }
  }
  #endif

  return reflect;
}
