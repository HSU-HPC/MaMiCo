#ifndef _MOLECULARDYNAMICS_COUPLING_SOLVERS_PRECICESOLVER_H_
#define _MOLECULARDYNAMICS_COUPLING_SOLVERS_PRECICESOLVER_H_

#include "coupling/CouplingMDDefinitions.h"
#include "tarch/la/Vector.h"
#include <cmath>
#include <sstream>
#include <fstream>
#include <iomanip>
#include "coupling/solvers/CouetteSolver.h"
#include "coupling/datastructures/MacroscopicCell.h"
#include "coupling/IndexConversion.h"
#include "precice/SolverInterface.hpp"
#include "coupling/interface/MacroscopicSolverInterface.h"
#include "coupling/indexing/CellIndex.h"

namespace coupling {
namespace solvers{

using namespace indexing;

template <unsigned int dim>
class PreciceAdapter: 
	public AbstractCouetteSolver<3>, 
	public interface::MacroscopicSolverInterface<dim> 
{
public:
	PreciceAdapter(
	  const double channelHeight,
	  const double dx,
	  const double dt,
	  const unsigned int plotEveryTimestep,
	  const std::string filestem,
	  const unsigned int overlap) :
	coupling::solvers::AbstractCouetteSolver<3>(),
	coupling::interface::MacroscopicSolverInterface<dim>(),
	_channelHeight(channelHeight),
	_dx(dx), 
	_dt(dt), 
	_plotEveryTimestep(plotEveryTimestep), 
	_filestem(filestem)
	{
	}

    virtual ~PreciceAdapter()
    {
		if (_interface!=NULL){_interface->finalize();delete _interface;_interface=NULL;}
    }
    
    void setInterface(const unsigned int overlap)
    {
    	_overlap=overlap;
    }
    
    void setCouplingMesh( tarch::la::Vector<3, double> mdDomainOffset,
						const tarch::la::Vector<3, double> macroscopicCellSize)
    {
	    _globalNumberMacroscopicCells = CellIndex<dim,IndexTrait::noGhost>().numberCellsInDomain;
    	_interface = new precice::SolverInterface("mamico","precice-config.xml",0,1);
		_numberOfCells = 0;
		for (CellIndex<dim, IndexTrait::vector, IndexTrait::noGhost> macroscopicCellIndex : coupling::indexing::CellIndex<dim,IndexTrait::vector, IndexTrait::noGhost>()) {
			// @TODO not sure the static_cast is casting what it should cast, it might be a giraffe.
			if (sendMacroscopicQuantityToMDSolver(static_cast<tarch::la::Vector<3, unsigned int>>(macroscopicCellIndex.get())))
				_numberOfCells++;
		}
		double *coords = new double[dim*_numberOfCells];
		unsigned int interfaceCellIndex=0;
		for (CellIndex<dim, IndexTrait::vector, IndexTrait::noGhost> macroscopicCellIndex : coupling::indexing::CellIndex<dim,IndexTrait::vector, IndexTrait::noGhost>()) {
			tarch::la::Vector<3, unsigned int> macroscopicCellVectorIndex = static_cast<tarch::la::Vector<3, unsigned int>>(macroscopicCellIndex.get());
			if (sendMacroscopicQuantityToMDSolver(macroscopicCellVectorIndex))
			{
				for (unsigned int currentDim = 0; currentDim < dim; currentDim++) {
					coords[dim*interfaceCellIndex+currentDim]=mdDomainOffset[currentDim] + macroscopicCellVectorIndex[currentDim]*_dx + 0.5*macroscopicCellSize[currentDim];
				}
				interfaceCellIndex++;
			}
		}

		/*
		std::cout << "[";
		for (unsigned int interfaceCellIndex = 0; interfaceCellIndex < _numberOfCells; interfaceCellIndex++) {
			std::cout << "[" << coords[dim*interfaceCellIndex] << "," << coords[dim*interfaceCellIndex+1] << "," << coords[dim*interfaceCellIndex+2] << "]";
		}
		std::cout << "]" << std::endl;
		*/
		
		int meshID = _interface->getMeshID("mamico-mesh");
		_vertexIDs = new int[_numberOfCells];
		_interface->setMeshVertices(meshID, _numberOfCells, coords, _vertexIDs); 
		delete[] coords;
		_precice_dt = _interface->initialize();
		std::cout << "PreciceSolver constructor called" << std::endl;
    }
    
	void setWallVelocity(const tarch::la::Vector<3,double> wallVelocity) override
	{
		_wallVelocity = wallVelocity;
	}
	
	tarch::la::Vector<3,double> getVelocity(tarch::la::Vector<3,double> pos) const override
	{
		tarch::la::Vector<3,double> vel(0.0);
		return vel;
	}
	
	double getDensity(tarch::la::Vector<3,double> pos) const
	{
		return 0.0;
	}
	
	void advance(double dt) override 
	{
		std::cout << "PreciceSolver::advance called with dt = " << dt << std::endl;
		if (_interface->isCouplingOngoing()) {
			int meshID = _interface->getMeshID("mamico-mesh");
			double* velocities = new double[_numberOfCells*dim];
			if (_interface->isReadDataAvailable())
			{
				// velocity from the conitnuum solver
				int dataID = _interface->getDataID("Velocity", meshID);
				_interface->readBlockVectorData(dataID, _numberOfCells, _vertexIDs, velocities);
			}
			
			// Solving the time step
			// Normally does nothing, everything is done on the MD side
			for (unsigned int vertexIndex = 0; vertexIndex < _numberOfCells; vertexIndex++)
			{
				velocities[dim*vertexIndex] = 1.0;
				velocities[dim*vertexIndex+1] = 1.0;
				velocities[dim*vertexIndex+2] = 1.0;
			}
			
			if (_interface->isWriteDataRequired(dt))
			{
				// Velocit from the md solver
				int dataID = _interface->getDataID("MDVelocity", meshID);
				_interface->writeBlockVectorData(dataID, _numberOfCells, _vertexIDs, velocities);
			}
			
			double computed_dt = std::min(_precice_dt, dt);
			_precice_dt = _interface->advance(computed_dt);
    	}
	}
	
	bool receiveMacroscopicQuantityFromMDSolver(tarch::la::Vector<dim, unsigned int> globalCellIndex) override
	{
		bool recv = true;
		for (unsigned int d = 0; d < dim; d++) {
		  recv = recv && (globalCellIndex[d] > _overlap) &&
		         (globalCellIndex[d] <
		          _globalNumberMacroscopicCells[d] + 1 - _overlap);
		}
		return recv;
	}
	
	bool sendMacroscopicQuantityToMDSolver(tarch::la::Vector<3, unsigned int> globalCellIndex) override
	{
		bool outer = false;
		for (unsigned int d = 0; d < 3; d++) {
		  outer = outer || (globalCellIndex[d] < 1) ||
		          (globalCellIndex[d] > _globalNumberMacroscopicCells[d]);
		}
		return (!outer) &&
		       (!receiveMacroscopicQuantityFromMDSolver(globalCellIndex));
	}
  
	std::vector<unsigned int> getRanks(tarch::la::Vector<dim, unsigned int> globalCellIndex) override 
	{
		std::vector<unsigned int> ranks;
		ranks.push_back(0);
		return ranks;
	  }

private:    
    const double _channelHeight;
    const double _dx;
    const double _dt;
    const unsigned int _plotEveryTimestep;
    const std::string _filestem;
    unsigned int _numberOfCells;
    unsigned int _counter{0};
	tarch::la::Vector<3,double> _wallVelocity;
	tarch::la::Vector<3, unsigned int> _globalNumberMacroscopicCells;
	unsigned int _overlap;
	
    precice::SolverInterface *_interface=NULL;
	double _precice_dt;
	int *_vertexIDs;
};
}
}

#endif
