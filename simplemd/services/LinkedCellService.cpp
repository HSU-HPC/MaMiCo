// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#include "simplemd/services/LinkedCellService.h"
#include "simplemd/services/MoleculeService.h"

simplemd::services::LinkedCellService::LinkedCellService(
  const tarch::la::Vector<MD_DIM,double>& domainSize,
  const tarch::la::Vector<MD_DIM,double>& domainOffset,
  const simplemd::services::ParallelTopologyService& parallelTopologyService,
  simplemd::services::MoleculeService& moleculeService
): _cells(NULL), _domainSize(domainSize), _domainOffset(domainOffset),
   _meshWidth(getMeshwidth(domainSize,parallelTopologyService.getLocalNumberOfCells())),
   _numberOfCells(parallelTopologyService.getLocalNumberOfCells()),_indexOffset(tarch::la::Vector<MD_DIM,unsigned int>(1)),
   _totalNumberOfCells(_numberOfCells+2u*_indexOffset)
   #if (MD_DIM>2)
   ,_totalNumberOfCells_X_By_totalNumberOfCells_Y(_totalNumberOfCells[0]*_totalNumberOfCells[1])
   #endif
{
  // initialise memory for cell structure according to number of cells and index of first cell on this process
  initCellStructure();

  // sort molecules into cells. This is done as an initial step.
  simplemd::moleculemappings::UpdateLinkedCellListsMapping updateLinkedCellLists(parallelTopologyService,*this);
  moleculeService.iterateMolecules(updateLinkedCellLists,false);
}


const tarch::la::Vector<MD_DIM,unsigned int>&
simplemd::services::LinkedCellService::getLocalIndexOfFirstCell() const {
  return _indexOffset;
}


const tarch::la::Vector<MD_DIM,unsigned int>&
simplemd::services::LinkedCellService::getLocalNumberOfCells() const {
  return _numberOfCells;
}



void simplemd::services::LinkedCellService::initCellStructure(){
  unsigned int numberCells = 1;
  if (_cells != NULL){
    delete [] _cells;
    _cells = NULL;
  }

  for (unsigned int d = 0; d < MD_DIM; d++){
    // increase number of cells by two to account for ghost layer
    numberCells = numberCells * _totalNumberOfCells[d];
  }

  _cells = new LinkedCell[numberCells];
  if (_cells ==NULL){
    std::cout << "ERROR simplemd::services::LinkedCellService::initCellStructure(): _cells==NULL!" << std::endl; exit(EXIT_FAILURE);
  }
}


void simplemd::services::LinkedCellService::shutdown(){
  if (_cells != NULL){
    delete [] _cells;
    _cells = NULL;
  }
}


unsigned int simplemd::services::LinkedCellService::getLocalIndexFromLocalVector(
  const tarch::la::Vector<MD_DIM,unsigned int>& coords
) const {
  return coords[0]
#if (MD_DIM > 1)
       + coords[1]*_totalNumberOfCells[0]
#endif
#if (MD_DIM > 2)
       + coords[2]*_totalNumberOfCells_X_By_totalNumberOfCells_Y
#endif
  ;
}


const tarch::la::Vector<MD_DIM,double>& simplemd::services::LinkedCellService::getMeshWidth() const {
  return _meshWidth;
}


const tarch::la::Vector<MD_DIM,double>& simplemd::services::LinkedCellService::getLocalDomainOffset() const {
  return _domainOffset;
}


const tarch::la::Vector<MD_DIM,double>& simplemd::services::LinkedCellService::getLocalDomainSize() const {
  return _domainSize;
}


void simplemd::services::LinkedCellService::addMoleculeToLinkedCell(Molecule &molecule, const tarch::la::Vector<MD_DIM,unsigned int> &localCellIndex){
#if (MD_DEBUG == MD_YES)
  for (unsigned int d = 0; d < MD_DIM; d++){
    if (localCellIndex[d] >= _totalNumberOfCells[d]){
      std::cout << "ERROR simplemd::services::LinkedCellService::addMoleculeToLinkedCell: localCellIndex out of range!" << std::endl;
      std::cout << "Cell index=" << localCellIndex << "; " << molecule.getConstPosition() << "; " << molecule.getConstVelocity() << ";" << molecule.getConstForce() << ";" << molecule.getConstForceOld() << std::endl;
      exit(EXIT_FAILURE);
    }
  }
#endif
  unsigned int index = localCellIndex[0]
#if (MD_DIM > 1)
                     + _totalNumberOfCells[0]*localCellIndex[1]
#endif
#if (MD_DIM > 2)
                     +   _totalNumberOfCells_X_By_totalNumberOfCells_Y*localCellIndex[2]
#endif
                     ;

  _cells[index].addMolecule(&molecule);
}


void simplemd::services::LinkedCellService::addMoleculeToLinkedCell(Molecule &molecule, const unsigned int &localCellIndex){
  _cells[localCellIndex].addMolecule(&molecule);
}


simplemd::LinkedCell& simplemd::services::LinkedCellService::getLinkedCell(const tarch::la::Vector<MD_DIM,unsigned int> &localCellIndex){
#if (MD_DEBUG==MD_YES)
  for (unsigned int d = 0; d < MD_DIM; d++){
    if (localCellIndex[d] >= _totalNumberOfCells[d]){
      std::cout << "ERROR simplemd::services::LinkedCellService::getLinkedCell: localCellIndex out of range!" << std::endl;
      exit(EXIT_FAILURE);
    }
  }
#endif
  unsigned int index = localCellIndex[0]
#if (MD_DIM>1)
                     + localCellIndex[1]*_totalNumberOfCells[0]
#endif
#if (MD_DIM > 2)
                     + localCellIndex[2]*_totalNumberOfCells_X_By_totalNumberOfCells_Y
#endif
                     ;
  return _cells[index];
}


void simplemd::services::LinkedCellService::deleteMoleculeFromLinkedCell(
Molecule& molecule, const tarch::la::Vector<MD_DIM,unsigned int> &localCellIndex){
  // compute index of respective cell
  unsigned int index = localCellIndex[0]
#if (MD_DIM>1)
                     + localCellIndex[1]*_totalNumberOfCells[0]
#endif
#if (MD_DIM > 2)
                     + localCellIndex[2]*_totalNumberOfCells_X_By_totalNumberOfCells_Y
#endif
                     ;
  _cells[index].deleteMolecule(&molecule);
}


bool simplemd::services::LinkedCellService::isGhostCell(const unsigned int &cellIndex) const {
  unsigned int help = cellIndex;
  for (int d = MD_DIM-1; d>-1; d--){
    unsigned int div = 1;
    for (int e = 0; e < d; e++){
      div = div* _totalNumberOfCells[e];
    }
    const unsigned int coord = help/div;
    // if the coordinate is at the beginning or end, return true; otherwise: consider next coordinate
    if ((coord==0) || (coord==_totalNumberOfCells[d]-1) ){
      return true;
    }
    help = help%div;
  }
  // return false if no ghost cell coordinate could be detected
  return false;
}


tarch::la::Vector<MD_DIM,unsigned int> simplemd::services::LinkedCellService::
getLocalCellIndexVector(const unsigned int cellIndex) const {
  unsigned int help = cellIndex;
  tarch::la::Vector<MD_DIM,unsigned int> localCellIndexVector(0);

  for (int d = MD_DIM-1; d>-1; d--){
    unsigned int div = 1;
    for (int e = 0; e < d; e++){
      div = div* _totalNumberOfCells[e];
    }
    localCellIndexVector[d] = help/div;

    help = help%div;
  }
  return localCellIndexVector;
}


unsigned int simplemd::services::LinkedCellService::
getLocalCellIndex(const tarch::la::Vector<MD_DIM,unsigned int> &cellIndexVector) const {
  unsigned int cellIndex = 0;

  for (unsigned int d = 0; d < MD_DIM; d++){
    unsigned int help = 1;
    for (unsigned int e = 0; e < d; e++){
      help = help*_totalNumberOfCells[e];
    }
    cellIndex += help*cellIndexVector[d];
  }
  return cellIndex;
}
