#ifndef LS1_MOLECULE_ITERATOR_H_
#define LS1_MOLECULE_ITERATOR_H_

//MD library
#include "ls1/src/molecules/Molecule.h"
#include "ls1/src/particleContainer/ParticleContainer.h"
#include "ls1/src/Simulation.h"
#include "ls1/src/particleContainer/RegionParticleIterator.h"

//utils
#include "tarch/la/Vector.h"

//MaMiCo
#include "coupling/interface/MoleculeIterator.h"
#include "LS1RegionWrapper.h"
#include "LS1Molecule.h"

namespace coupling
{
    namespace interface
    {
        class LS1MoleculeIterator;
    }
}

class coupling::interface::LS1MoleculeIterator : public coupling::interface::MoleculeIterator<ls1::LS1RegionWrapper, 3>
{
public:
    LS1MoleculeIterator(ls1::LS1RegionWrapper& cell) : coupling::interface::MoleculeIterator<ls1::LS1RegionWrapper, 3>(cell) { }

    /** sets the iterator to the first element */
    virtual void begin() { _cell.iteratorReset(); }

    /** returns true, if the iterator should continue the molecule traversal */
    virtual bool continueIteration() const { return _cell.iteratorValid(); }

    /** sets the iterator to the next molecule */
    virtual void next() { _cell.iteratorNext(); }

    /** returns a reference to the molecule that this iterator currently points to */
    virtual coupling::interface::Molecule<3>& get() { _tempMolecule.setMolecule(_cell.getParticleAtIterator()); return _tempMolecule; }

    /** returns a const reference to the current molecule for pure reading purposes */
    virtual const coupling::interface::Molecule<3>& getConst() { _tempMolecule.setMolecule(_cell.getParticleAtIterator()); return _tempMolecule; }

private:
    coupling::interface::LS1Molecule _tempMolecule;
};

#endif