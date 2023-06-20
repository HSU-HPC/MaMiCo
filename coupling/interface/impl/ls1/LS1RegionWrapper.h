// Maps regions in the LS1 MD domain to mamico cells
// Causes behavior identical to fixed regular cells
#ifndef LS1_REGION_WRAPPER_H_
#define LS1_REGION_WRAPPER_H_


#include "ls1/src/particleContainer/RegionParticleIterator.h"
#include "ls1/src/particleContainer/ParticleContainer.h"
#include "ls1/src/molecules/Molecule.h"
#include "ls1/src/Simulation.h"
#include "ls1/src/parallel/DomainDecompMPIBase.h"
//#include "ls1/ensemble/EnsembleBase.h"


#include "coupling/interface/Molecule.h"
#include "coupling/interface/impl/ls1/LS1StaticCommData.h"

namespace ls1 {

class LS1RegionWrapper
{
public:
    LS1RegionWrapper(double startRegion[3], double endRegion[3], Simulation* simulation)
    {
        _curParticleID = 0;
        _IDinited = false;
        for(int i = 0; i < 3; i++)
        {
            _startRegion[i] = startRegion[i];
            _endRegion[i] = endRegion[i]; //boundingboxmin
        }
        _locSimulation = simulation;
        _particleContainer = _locSimulation->getMoleculeContainer();
        _iterator = _particleContainer->regionIterator(_startRegion, _endRegion, ParticleIterator::ONLY_INNER_AND_BOUNDARY);
    }
    LS1RegionWrapper() : _startRegion{0,0,0}, _endRegion{0,0,0}, _curParticleID(0), _IDinited(false) {}

    void setRegion(double startRegion[3], double endRegion[3])
    {
        for(int i = 0; i < 3; i++)
        {
            _startRegion[i] = startRegion[i];
            _endRegion[i] = endRegion[i];
        }
        _iterator = _particleContainer->regionIterator(_startRegion, _endRegion, ParticleIterator::ONLY_INNER_AND_BOUNDARY);
    }

    void iteratorReset()
    {
        _iterator = _particleContainer->regionIterator(_startRegion, _endRegion, ParticleIterator::ONLY_INNER_AND_BOUNDARY);
    }

    void iteratorNext() { ++_iterator; }

    bool iteratorValid() { return _iterator.isValid(); }

    ::Molecule* getParticleAtIterator()
    {
        ::Molecule* temp = &(*_iterator);
        return temp;
    }

    bool isInRegion(const double point[3])
    {
        bool isInRegion = true;
        for (int i = 0; i < 3; i++) { isInRegion &= ((point[i] >= _startRegion[i]) & (point[i] < _endRegion[i])); }
        return isInRegion;
    }
    bool isInRegion(const tarch::la::Vector<3, double> point)
    {
        return isInRegion(new double[3]{point[0], point[1], point[2]});
    }

    bool isInRegion(const double startBox[3], const double endBox[3])
    {
        bool isInRegion = true;
        for (int i = 0; i < 3; i++) { isInRegion &= ((startBox[i] >= _startRegion[i]) & (endBox[i] < _endRegion[i])); }
        return isInRegion;
    }

    void addMolecule(::Molecule &molecule)
    {
        _particleContainer->addParticle(molecule);
    }

    void deleteMolecule(ParticleIterator &iterator)
    {
        _particleContainer->deleteMolecule(iterator, false);
    }

    void addMolecule(const coupling::interface::Molecule<3> &molecule)
    {
        ::Molecule temp;

        auto position = molecule.getPosition();
        auto velocity = molecule.getVelocity();
        auto force = molecule.getForce();

        temp.setr(0, position[0] - coupling::interface::LS1StaticCommData::getInstance().getBoxOffsetAtDim(0)); //temporary till ls1 offset is natively supported
        temp.setr(1, position[1] - coupling::interface::LS1StaticCommData::getInstance().getBoxOffsetAtDim(1));
        temp.setr(2, position[2] - coupling::interface::LS1StaticCommData::getInstance().getBoxOffsetAtDim(2));
        
        temp.setv(0, velocity[0]);
        temp.setv(1, velocity[1]);
        temp.setv(2, velocity[2]);

        temp.setF(0, force[0]);  //more important for usher particles
        temp.setF(1, force[1]);
        temp.setF(2, force[2]);
        
        if(!_IDinited)
        {
            _curParticleID = _locSimulation->getTotalNumberOfMolecules() + 1;
            _IDIncrementor = 1;
            
            #ifdef ENABLE_MPI
                int curRank = global_simulation->domainDecomposition().getRank();
                _curParticleID += curRank + 1;
                _IDIncrementor = global_simulation->domainDecomposition().getNumProcs();
            #endif

            _IDinited = true;
        }

        temp.setid(_curParticleID);
        _curParticleID += _IDIncrementor;

        addMolecule(temp);
        //_particleContainer->addParticle(temp);
    }

    void deleteMolecule(const coupling::interface::Molecule<3> &molecule)
    {
        //check if coords in region
        auto molPosition = molecule.getPosition();
        for(int i = 0; i < 3; i++) {molPosition[i] = molPosition[i] - coupling::interface::LS1StaticCommData::getInstance().getBoxOffsetAtDim(i);} //temporary till ls1 offset is natively supported
        if(!isInRegion(molPosition))
            return;
        //check if molecule at location specified
        double cutoff = _locSimulation->getcutoffRadius();

        double startBox[] = {molPosition[0]-cutoff/10, molPosition[1]-cutoff/10, molPosition[2]-cutoff/10};
        double endBox[] = {molPosition[0]+cutoff/10, molPosition[1]+cutoff/10, molPosition[2]+cutoff/10};

        RegionParticleIterator temp = _particleContainer->regionIterator(startBox, endBox, ParticleIterator::ONLY_INNER_AND_BOUNDARY);
        bool found = false;
        while(temp.isValid())
        {
            ::Molecule* temp = &(*_iterator);
            if(temp->r(0) == molPosition[0] && temp->r(1) == molPosition[1] && temp->r(2) == molPosition[2])
            {
                found = true;
                break;
            }
            temp++;
        }
        //delete molecule
        if(!found)
            return;

        deleteMolecule(temp);
        //_particleContainer->deleteMolecule(temp, false);
    }

private:
    double _startRegion[3], _endRegion[3];
    unsigned long int _curParticleID;
    int _IDIncrementor;
    bool _IDinited;
    Simulation* _locSimulation;
    RegionParticleIterator _iterator;
    ParticleContainer* _particleContainer;
};
}

#endif