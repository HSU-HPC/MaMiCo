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
    LS1Molecule(::Molecule *myMolecule): _myMolecule(myMolecule), _u(0.0){}
    LS1Molecule(): _myMolecule(), _u(0.0){}
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
        return _u;
    }
    virtual void setPotentialEnergy(const double& potentialEnergy) 
    {
        _u = potentialEnergy;
    }

private:
    ::Molecule* _myMolecule;
    double _u;
};
#endif