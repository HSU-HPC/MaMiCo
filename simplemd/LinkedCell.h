// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_LINKEDCELL_H_
#define _MOLECULARDYNAMICS_LINKEDCELL_H_

#include<list>
// #include<vector>
#include<iostream>
#include<fstream>
#include<string>
#include<cstdlib>
#include "simplemd/Molecule.h"

namespace simplemd {
  class LinkedCell;
}

/** describes a linked cell.
 *  @author Philipp Neumann
 */
class simplemd::LinkedCell {
  public:
    /** initialises linked cell list with numberMolecules empty entries */
    LinkedCell(const unsigned int numberMolecules=0): _molecules(numberMolecules), _ghostMolecules(0){}
    ~LinkedCell(){
      _molecules.clear();
      _ghostMolecules.clear();
    }

    /** iterators to begin and end position */
    std::list<Molecule*>::iterator begin(){ return _molecules.begin();}
    std::list<Molecule*>::const_iterator constBegin() const{ return _molecules.begin();}
    std::list<Molecule*>::iterator end(){ return _molecules.end();}
    std::list<Molecule*>::const_iterator constEnd() const{ return _molecules.end();}

    std::list<Molecule*>::iterator beginGhost(){ return _ghostMolecules.begin();}
    std::list<Molecule*>::const_iterator constBeginGhost() const{ return _ghostMolecules.begin();}
    std::list<Molecule*>::iterator endGhost(){ return _ghostMolecules.end();}
    std::list<Molecule*>::const_iterator constEndGhost() const{ return _ghostMolecules.end();}
    std::list<Molecule*>::iterator beginGhost(const unsigned int dim){
      return _ghostIterators[dim];
    }

    std::list<Molecule*>::iterator endGhost(const unsigned int dim){
      return _ghostIterators[dim+3];
    }

    std::list<Molecule*>& getList(){return _molecules;}
    const std::list<Molecule*>& getConstList(){return _molecules;}

    /** adds a molecule pointer */
    void addMolecule(Molecule *molecule){
      _molecules.push_back(molecule);
    }

    void addGhostMolecule(Molecule *molecule, const unsigned int ghostIteratorIndex = 0){
      std::string filename = "ghostParticles.txt";
      std::ofstream file(filename, std::ios_base::app );
      if (!file.is_open()){std::cout << "ERROR CouetteTest::write2CSV(): Could not open file " << filename << "!" << std::endl; exit(EXIT_FAILURE);}
      auto position = molecule->getPosition();
      auto velocity = molecule->getVelocity();
      file << position[0] << " " << position[1] << " " << position[2] << " " << ghostIteratorIndex<< std::endl;
      file.close();
      auto myMolecule = new Molecule(position,velocity);
      _ghostMolecules.push_back(myMolecule);
      if(ghostIteratorIndex>0){
        _ghostIterators[ghostIteratorIndex-1] = _ghostMolecules.begin();
        std::advance(_ghostIterators[ghostIteratorIndex-1],  _ghostMolecules.size() - 1);
      }
    }

    /** deletes the molecule 'molecule' from the list, if it is contained */
    void deleteMolecule(Molecule *molecule){
      #if (MD_ERROR == MD_YES)
      if (molecule==NULL){
        std::cout << "ERROR simplemd::LinkedCell::deleteMolecule: molecule==NULL!" << std::endl;
        exit(EXIT_FAILURE);
      }
      #endif
      _molecules.remove(molecule);
    }

    /** initialises the molecule pointer at it with the value of molecule */
    void setMolecule(std::list<Molecule*>::iterator& it, Molecule* molecule){
      *it = molecule;
    }

    bool hasGhostMolecules()const{
      return !_ghostMolecules.empty();
    }

  private:
    std::list<Molecule*> _molecules;
    std::list<Molecule*> _ghostMolecules;
    tarch::la::Vector<6, std::list<Molecule*>::iterator> _ghostIterators;
};

#endif // _MOLECULARDYNAMICS_LINKEDCELL_H_
