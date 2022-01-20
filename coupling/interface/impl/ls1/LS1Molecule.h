#ifndef LS1_MOLECULE_H_
#define LS1_MOLECULE_H_

//MD library
#include "ls1/src/molecules/Molecule.h"
#include "ls1/src/Simulation.h"
#include "ls1/src/ensemble/EnsembleBase.h"

//utils
#include "tarch/la/Vector.h"

//MaMiCo
#include "coupling/interface/Molecule.h"
#include "coupling/interface/impl/ls1/LS1RegionWrapper.h"
#include "coupling/interface/impl/ls1/LS1StaticCommData.h"

namespace coupling
{
    namespace interface
    {
        class LS1Molecule;
    }
}

class coupling::interface::LS1Molecule: public coupling::interface::Molecule<3>
{
public:
    LS1Molecule(::Molecule *myMolecule): _myMolecule(myMolecule){}
    LS1Molecule(): _myMolecule(){}
    virtual ~LS1Molecule(){}

    void setMolecule(::Molecule *newMolecule) { _myMolecule=newMolecule;}

    /** returns/ sets the velocity of the molecule */
    virtual tarch::la::Vector<3,double> getVelocity() const
    {
        tarch::la::Vector<3, double> velocity(0.0);
      	for(int i = 0; i < 3; i++)
        {	
			velocity[i] = _myMolecule->v(i);
		}
      	return velocity;
    }
    virtual void setVelocity(const tarch::la::Vector<3,double>& velocity)
    {
        for(int i = 0; i < 3; i++)
        {
            _myMolecule->setv(i, velocity[i]);
        }
    }

    /** returns/ sets the position of the molecule */
    virtual tarch::la::Vector<3,double> getPosition() const
    {
        tarch::la::Vector<3, double> position(0.0);
        for(int i = 0; i < 3; i++)
        {
            position[i] = _myMolecule->r(i) + coupling::interface::LS1StaticCommData::getInstance().getBoxOffsetAtDim(i); //temporary till ls1 offset is natively supported
        }
        return position;
    }
    virtual void setPosition(const tarch::la::Vector<3,double>& position)
    {
        for(int i = 0; i < 3; i++)
        {
            _myMolecule->setr(i, position[i] - coupling::interface::LS1StaticCommData::getInstance().getBoxOffsetAtDim(i)); //temporary till ls1 offset is natively supported
        }        
    }

    /** sets the force acting on this molecule. This function is only called in the USHER
     *  scheme so far If you want to set the force of a newly created molecule,
     *  you need to implement this function.
     */
    virtual void setForce(const tarch::la::Vector<3,double>& force)
    {
        for(int i = 0; i < 3; i++)
        {
            _myMolecule->setF(i, force[i]);
        }
    }
    virtual tarch::la::Vector<3,double> getForce() const
    {
        tarch::la::Vector<3, double> force(0.0);
        for(int i = 0; i < 3; i++)
        {
            force[i] = _myMolecule->F(i);
        }
        return force;
    }

    /** returns/ sets the potential energy of the molecule */
    virtual double getPotentialEnergy() const
    {
        double u = 0.0;

        tarch::la::Vector<3,double> location1, location2;
        //if(_myMolecule == NULL) return 0.0;
        
        location1 = { _myMolecule->r(0), _myMolecule->r(1), _myMolecule->r(2)};
        
        //find all molecules within cutoff
        double cutoff = global_simulation->getcutoffRadius();
        Ensemble* ensemble = global_simulation->getEnsemble();
        const double sigma = ensemble->getComponent(0)->getSigma(0);
        const double epsilon = ensemble->getComponent(0)->ljcenter(0).eps();

        const double sigma2 = sigma * sigma;
        const double sigma6 = sigma2 * sigma2 * sigma2;

        double startRegion[] = {location1[0] - cutoff, location1[1] - cutoff, location1[2] - cutoff};
        double endRegion[] = {location1[0] + cutoff, location1[1] + cutoff, location1[2] + cutoff};

        ls1::LS1RegionWrapper region(startRegion, endRegion);
        double cutoff2 = cutoff * cutoff;

        //calculate lennard jones energy
        while(region.iteratorValid())
        {
            ::Molecule* temp = region.getParticleAtIterator();
            location2 = { temp->r(0), temp->r(1), temp->r(2) };
            const double r2 = tarch::la::dot(location2 - location1, location2 - location1);
            if(r2 < cutoff2)
            {
                const double r6 = r2 * r2 * r2;
                const double contrib =  2.0* epsilon * (sigma6/r6) * ((sigma6/r6) - 1.0);
                u += contrib;
            }

            region.iteratorNext();
        }

        //sum and return
        return u;
    }
    virtual void setPotentialEnergy(const double& potentialEnergy) {};

private:
    ::Molecule* _myMolecule;
};
#endif