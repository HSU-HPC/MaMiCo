// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_ZHOUBOUNDARYFORCE_H_
#define _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_ZHOUBOUNDARYFORCE_H_

#include <cmath>
#include <iostream>
#include "coupling/interface/MDSolverInterface.h"
#include "coupling/interface/Molecule.h"

namespace coupling {
  namespace cellmappings {
    template<class LinkedCell,unsigned int dim>
    class ZhouBoundaryForce;
  }
}

/** applies the Zhou boundary force to all molecules assuming a cut-off radius r_c=2.5.
 *  We only consider molecules in the outermost macroscopic cell. Thus, the current method will only be
 *  employed to molecules within a distance=min(r_c,size of macroscopic cell).
 *  The force is based on the research paper:
 *  W.J. Zhou, H.B. Luan, Y.L. He, J. Sun, W.Q. Tao.
 *  A study on boundary force model used in multiscale simulations with non-periodic boundary condition
 *  Microfluid Nanofluid 16(3): 587-595, 2014
 *
 *  The force is computed based on an interpolation formula which has been obtained from various MD simulation
 *  runs at different temperature and density values. The investigated range is:
 *  0.4*m*sigma-3<density<0.9*m*sigma-3, 1.3*eps/kB<temperature<3.9*eps/kB.
 *  Moreover, only 3D simulations have been considered in this study.
 *	@brief This class applies the Zhou boundary force to all molecules assuming a cut-off radius r_c=2.5.
 *	@tparam LinkedCell cell type
 *	@tparam dim Number of dimensions; it can be 1, 2 or 3
 *  @author Philipp Neumann
 *	\todo Philipp please take a look on this class
 */
template<class LinkedCell,unsigned int dim>
class coupling::cellmappings::ZhouBoundaryForce {
  public:
    /** currently, we support constant temperature and density and 3D only. To be extended
	 *	@param density
	 *	@param temperature
	 *	@param epsilon
	 *	@param sigma
	 *	@param boundary
	 *	@param domainOffset
	 *	@param domainSize
	 *	@param mdSolverInterface
	 */
    ZhouBoundaryForce(
      const double &density,
      const double &temperature,
      const double &epsilon,
      const double &sigma,
      const tarch::la::Vector<2*dim,bool>& boundary,
      const tarch::la::Vector<dim,double>& domainOffset,
      const tarch::la::Vector<dim,double>& domainSize,
      coupling::interface::MDSolverInterface<LinkedCell,dim> * const mdSolverInterface
    ): _boundary(boundary),_domainLowerLeft(domainOffset),_domainUpperRight(domainSize+domainOffset),
      _mdSolverInterface(mdSolverInterface),
      _p1(getP1(density,temperature)),
      _p2(getP2(density,temperature)),
      _p3(getP3(density,temperature)),
      _q1(getQ1(density)),
      _q2(getQ2(density)),
      _q3(getQ3(density)),
      _forceFactor(epsilon/sigma),
      _sigma(sigma),
      /** Fixed step width for numerical integration of boundary force term */
      _energyResolution(1000),
      _energyTable(new double[_energyResolution]){
        // simple trapezoidal integration scheme
        double force(0);
        double energy(0);
        double dx(2.5 * sigma / (_energyResolution-1));
        for(unsigned int i=0;i<_energyResolution;i++){
          force = getScalarBoundaryForce(2.5 - i*dx);
          _energyTable[i] = energy + force/2*dx;
          energy += force * dx;
        }

        // for(unsigned int i=0;i<_energyResolution;i++){
        //   std::cout << "_energyTable[" << i << "] = " << _energyTable[i] << std::endl;
        // }
        // for(double pos = 0; pos <= 2.5; pos += .25/_energyResolution){
        //   std::cout << "getPotentialEnergy(" << pos << ") = " << getPotentialEnergy(pos) << std::endl;
        // }
    }

    /** Destructor */
	~ZhouBoundaryForce(){delete [] _energyTable;}

    /** empty function
	 */
	void beginCellIteration(){}

    /** empty function
	 */
	void endCellIteration(){}

    /** extracts position and force of each molecule, add boundary force and applies it to Molecule
	 *	@param cell
	 *	@param cellIndex
	 */
	void handleCell(LinkedCell& cell,const unsigned int &cellIndex){
      coupling::interface::MoleculeIterator<LinkedCell,dim> *it = _mdSolverInterface->getMoleculeIterator(cell);
      it->begin();
      while(it->continueIteration()){
        coupling::interface::Molecule<dim> &wrapper(it->get());
        //std::cout << "apply force to Molecule " << wrapper.getPosition() << std::endl;

        // extract position and force of each molecule
        const tarch::la::Vector<dim,double> position(wrapper.getPosition());
        tarch::la::Vector<dim,double> force(wrapper.getForce());
        // add boundary force
        force = force + getBoundaryForces(position);
        wrapper.setForce(force);

        it->next();
      }
      delete it;
    }

    /** returns the potential energy
	 *	@param position
	 *	@return energy
     */
	double getPotentialEnergy(const tarch::la::Vector<dim,double>& position) const {
      double energy(0);
      const double distance = 2.5*_sigma;
      for (unsigned int d = 0; d < dim; d++){
        if (_boundary[2*d]   && (position[d]-_domainLowerLeft[d]<distance)) { energy += getPotentialEnergy((position[d]-_domainLowerLeft[d])/_sigma);}
        if (_boundary[2*d+1] && (_domainUpperRight[d]-position[d]<distance)){ energy += getPotentialEnergy((_domainUpperRight[d]-position[d])/_sigma);}
      }
      return energy;
    }

    /** checks the boundary flags for open boundaries. If an open boundary is encountered and this molecule is close to that boundary (distance smaller than 2.5, cf. Zhou paper),
     *  the respective force contribution in that dimension is added/subtracted.
	 *	@param position
	 *	@return force
     */
    tarch::la::Vector<dim,double> getBoundaryForces(const tarch::la::Vector<dim,double>& position) const {
      tarch::la::Vector<dim,double> force(0.0);
      const double distance = 2.5*_sigma;

      for (unsigned int d = 0; d < dim; d++){
        // for left/lower/back boundary: add force; for right/upper/front boundary: subtract force value
        if (_boundary[2*d]   && (position[d]-_domainLowerLeft[d]<distance)) { force[d] += _forceFactor*getScalarBoundaryForce( (position[d]-_domainLowerLeft[d])/_sigma);   }
        if (_boundary[2*d+1] && (_domainUpperRight[d]-position[d]<distance)){ force[d] -= _forceFactor*getScalarBoundaryForce((_domainUpperRight[d]-position[d])/_sigma);   }
      }
      return force;
    }

  private:
    
    /** perform interpolation between nearest energy lookup table values
	 *	@param rw
	 *	@return nearest energy to the ????
	 */
    double getPotentialEnergy(double rw) const {
      double dx(2.5 * _sigma / (_energyResolution-1));
      int upper = (int)((2.5-rw)/2.5*(_energyResolution-1)-.5);
      int lower = (int)((2.5-rw)/2.5*(_energyResolution-1)+.5);
      return _energyTable[lower] * ((2.5 - upper * dx - rw) / dx) 
           + _energyTable[upper] * ((rw - 2.5 + lower * dx) / dx);
    }

    /** evaluates the force expression for a scalar component of the force vector. We thus add up dimensional contributions if several boundaries are located beside each other.
	 *	@param rw
	 *	@return Scalar boundary force
	 */
    double getScalarBoundaryForce(double rw) const {
      if (rw < 1.04){ return _p1 + _p2*exp(pow((rw+0.25),3.4))*cos(_p3*rw) ;   }
      else { return -1.0/(_q1+_q2*(2.5-rw)*(2.5-rw) + _q3/((2.5-rw)*(2.5-rw))); }
    }

    /** computes and returns the coefficient p1 according to Zhou-paper
	 *	@param density
	 *	@param temperature
	 *	@return p1
	 */
    double getP1(const double& density, const double& temperature) const {
      const double logrho = log(density);
      const double T2     = temperature*temperature;
      const double T3     = T2*temperature;

      // an additional factor 10**-5 is in the expression 4.599/(...) -> private communication with Wenjing Zhou (typo in paper)
      const double p1 = (-18.953 + 53.369*temperature - 1.253*T2 + 4.599/100000.0*T3
                         +59.871*logrho + 19.737*logrho*logrho)/
                        (1.0 + 2.592*temperature - 0.557*T2 + 0.049*T3
                        -13.912*logrho + 18.657*logrho*logrho );
      return p1;
    }
	
	/** computes and returns the coefficient p2 according to Zhou-paper
	 *	@param density
	 *	@param temperature
	 *	@return p2
	 */
    double getP2(const double& density,const double &temperature) const {
      const double logrho = log(density);
      const double T2     = temperature*temperature;
      const double T3     = T2*temperature;

      const double p2 = (-0.094 + 2.808*temperature - 0.019*T2 - 0.001*T3
                         +2.823*logrho + 2.071*logrho*logrho)/ (1.0 + 0.168*temperature - 0.013*T2
                         -4.323*logrho + 2.557*logrho*logrho - 2.155*logrho*logrho*logrho);
      return p2;
    }
	
	/** computes and returns the coefficient p3 according to Zhou-paper
	 *	@param density
	 *	@param temperature
	 *	@return p3
	 */
    double getP3(const double& density,const double &temperature) const {
      const double T394 = pow(temperature,0.394);
      const double rho17437 = pow(density,17.437);
      const double p3 = 3.934 + 0.099*T394 - 0.097*rho17437 + 0.075*T394 * rho17437;
      return p3;
    }
	
	/** computes and returns the coefficient q1 according to Zhou-paper
	 *	@param density
	 *	@return q1
	 */
    double getQ1(const double& density) const {
      const double rho2 = density*density;
      const double rho3 = rho2*density;
      const double rho4 = rho2*rho2;
      return -30.471 + 113.879*density - 207.205*rho2 + 184.242*rho3 - 62.879*rho4;
    }
	
	/** computes and returns the coefficient q2 according to Zhou-paper
	 *	@param density
	 *	@return q2
	 */
    double getQ2(const double& density) const {
      const double rho2 = density*density;
      const double rho3 = rho2*density;
      const double rho4 = rho2*rho2;
      return 6.938 - 25.788*density + 46.773*rho2 - 41.768*rho3 + 14.394*rho4;
    }
	
	/** computes and returns the coefficient q3 according to Zhou-paper
	 *	@param density
	 *	@return q3
	 */
    double getQ3(const double& density) const {
      const double rho2 = density*density;
      const double rho3 = rho2*density;
      const double rho4 = rho2*rho2;
      return 39.634 - 147.821*density + 269.519*rho2 - 239.066*rho3 + 81.439*rho4;
    }

     /** an entry is true, if this is an open boundary. For enumeration of boundaries, see BoundaryForceConfiguration */
	const tarch::la::Vector<2*dim,bool> _boundary;   
    /** lower left corner of MD domain */
    const tarch::la::Vector<dim,double> _domainLowerLeft;
	/** upper right corner of MD domain */
    const tarch::la::Vector<dim,double> _domainUpperRight;

    coupling::interface::MDSolverInterface<LinkedCell,dim> * const _mdSolverInterface;
    // helper variables to speed up evaluation, see Zhou paper
    const double _p1;
    const double _p2;
    const double _p3;
    const double _q1;
    const double _q2;
    const double _q3;

    /** factor to scale the force according to MD units */
	const double _forceFactor;
	/** factor to scale length according to MD units */
    const double _sigma;

    const unsigned int _energyResolution;
    double* _energyTable;
};
#endif // _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_ZHOUBOUNDARYFORCE_H_
