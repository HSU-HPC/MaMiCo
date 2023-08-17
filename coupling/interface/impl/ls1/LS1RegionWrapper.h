// Maps regions in the LS1 MD domain to mamico cells
// Causes behavior identical to fixed regular cells
#ifndef LS1_REGION_WRAPPER_H_
#define LS1_REGION_WRAPPER_H_

#include <tuple>

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
    LS1RegionWrapper(double startRegion[3], double endRegion[3], Simulation* simulation) :
        _locSimulation(simulation),
        _particleContainer(simulation->getMoleculeContainer()),
        _cutoff(simulation->getcutoffRadius()),
        _cutoff2(_cutoff*_cutoff),
        _sigma(simulation->getEnsemble()->getComponent(0)->getSigma(0)),
        _sigma6(_sigma*_sigma*_sigma*_sigma*_sigma*_sigma),
        _epsilon(simulation->getEnsemble()->getComponent(0)->ljcenter(0).eps()),
        _cutoffEnergy(2.0 * _epsilon * _sigma6 / (_cutoff2 * _cutoff2 * _cutoff2) * (_sigma6 / (_cutoff2* _cutoff2 * _cutoff2) - 1.0))
    {
        _curParticleID = 0;
        _IDinited = false;
        for(int i = 0; i < 3; i++)
        {
            _startRegion[i] = startRegion[i];
            _endRegion[i] = endRegion[i]; //boundingboxmin
        }
        _iterator = simulation->getMoleculeContainer()->regionIterator(_startRegion, _endRegion, ParticleIterator::ALL_CELLS);
    }

    void setRegion(double startRegion[3], double endRegion[3])
    {
        for(int i = 0; i < 3; i++)
        {
            _startRegion[i] = startRegion[i];
            _endRegion[i] = endRegion[i];
        }
        _iterator = _particleContainer->regionIterator(_startRegion, _endRegion, ParticleIterator::ALL_CELLS);
    }

    void iteratorReset()
    {
        _iterator = _particleContainer->regionIterator(_startRegion, _endRegion, ParticleIterator::ALL_CELLS);
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

        temp.setComponent(_locSimulation->getEnsemble()->getComponent(0));

        addMolecule(temp);
        //_particleContainer->addParticle(temp);
    }

    void deleteMolecule(const coupling::interface::Molecule<3> &molecule)
    {
        //check if coords in region
        auto molPosition = molecule.getPosition();
        for(int i = 0; i < 3; i++) {molPosition[i] = molPosition[i] - coupling::interface::LS1StaticCommData::getInstance().getBoxOffsetAtDim(i);} //temporary till ls1 offset is natively supported;
        if(!isInRegion(molPosition))
            return;
        //check if molecule at location specified
        double margin = _cutoff/10;
        double startBox[] = {molPosition[0]-margin, molPosition[1]-margin, molPosition[2]-margin};
        double endBox[] = {molPosition[0]+margin, molPosition[1]+margin, molPosition[2]+margin};

        RegionParticleIterator _curIterator = _particleContainer->regionIterator(startBox, endBox, ParticleIterator::ALL_CELLS);
        bool found = false;
        while(_curIterator.isValid())
        {
            ::Molecule* temp = &(*_curIterator);
            if(temp->r(0) == molPosition[0] && temp->r(1) == molPosition[1] && temp->r(2) == molPosition[2])
            {
                found = true;
                break;
            }
            ++_curIterator;
        }
        //delete molecule
        if(!found)
            return;

        deleteMolecule(_curIterator);
        //_particleContainer->deleteMolecule(temp, false);
    }

    std::tuple<tarch::la::Vector<3,double>, double> calculateForceAndPotentialAtPoint(const tarch::la::Vector<3,double> position, double adjustCutoff)
    {
        tarch::la::Vector<3,double> force (0.0);
        double potentialEnergy = 0.0;

        //molecule position
        tarch::la::Vector<3,double> moleculePosition = position;
        for(int i = 0; i < 3; i++)
        {
            moleculePosition[i] = moleculePosition[i] - coupling::interface::LS1StaticCommData::getInstance().getBoxOffsetAtDim(i); //temporary till ls1 offset is natively supported
        }

        tarch::la::Vector<3,double> tempMoleculePosition;

        //calculate force
        //find all molecules within cutoff

        double startRegion[] = {moleculePosition[0] - _cutoff, moleculePosition[1] - _cutoff, moleculePosition[2] - _cutoff};
        double endRegion[] = {moleculePosition[0] + _cutoff, moleculePosition[1] + _cutoff, moleculePosition[2] + _cutoff};

        ls1::LS1RegionWrapper region(startRegion, endRegion, global_simulation);

        //calculate lennard jones energy
        while(region.iteratorValid())
        {
            ::Molecule* temp = region.getParticleAtIterator();
            tempMoleculePosition = { temp->r(0), temp->r(1), temp->r(2) };
            const auto r = tempMoleculePosition - moleculePosition;
            const double r2 = tarch::la::dot(r, r);
            if(r2 < _cutoff2)
            {
                const double r6 = r2 * r2 * r2;
                const auto forceContrib =  (24.0 * _epsilon / r2 * (_sigma6 / r6)) * (1.0 - 2.0 * (_sigma6 / r6)) * r;
                const double uContrib =  2.0* _epsilon * (_sigma6 / r6) * ((_sigma6 / r6) - 1.0) - (adjustCutoff?_cutoffEnergy:0);
                potentialEnergy += uContrib;
                force += forceContrib;
            }

            region.iteratorNext();
        }
        return std::make_tuple(force, potentialEnergy);
    }


private:
    double _startRegion[3], _endRegion[3];
    unsigned long int _curParticleID;
    int _IDIncrementor;
    bool _IDinited;
    const double _cutoff, _cutoff2;
    const double _sigma, _sigma6;
    const double _epsilon;
    const double _cutoffEnergy;
    Simulation* _locSimulation;
    RegionParticleIterator _iterator;
    ParticleContainer* _particleContainer;
};
}

#endif