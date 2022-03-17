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
#include <filesystem>

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
		if (_coords!=NULL){delete[] _coords;}
		if (_velocities!=NULL){delete[] _velocities;}
    }
    
    void setInterface(const unsigned int overlap)
    {
    	_overlap=overlap;
    }
    
    void setCouplingMesh(const tarch::la::Vector<3, double> mdDomainOffset,
						const tarch::la::Vector<3, double> macroscopicCellSize)
    {
    	_mdDomainOffset = mdDomainOffset;
    	_macroscopicCellSize = macroscopicCellSize;
    	_interface = new precice::SolverInterface("mamico","../precice-config.xml",0,1);
		_numberOfCells = 0;
		tarch::la::Vector<3, int> lowerBoundary = CellIndex<3>::lowerBoundary.get();
		tarch::la::Vector<3, int> upperBoundary = CellIndex<3>::upperBoundary.get();
		for (CellIndex<dim, IndexTrait::vector> cellIndex : CellIndex<dim, IndexTrait::vector>()) 
		{
			tarch::la::Vector<3, int> cellVectorIndex = cellIndex.get();
			bool isInnerCell = true;
			bool isGhostCell = false;
			for (unsigned int currentDim = 0; currentDim < dim; currentDim++) 
			{
				isInnerCell &= cellVectorIndex[currentDim] >= lowerBoundary[currentDim] + (int)_overlap + 1;
				isInnerCell &= cellVectorIndex[currentDim] <= upperBoundary[currentDim] - (int)_overlap - 1 ;
				isGhostCell |= cellVectorIndex[currentDim] < 1;
				isGhostCell |= cellVectorIndex[currentDim] >= upperBoundary[currentDim] ;				
			}
			if (!isInnerCell && !isGhostCell)
				_numberOfCells++;
		}
		std::cout << "Number of coupling volumes: " << _numberOfCells << std::endl;
		_coords = new double[dim*_numberOfCells];
		unsigned int interfaceCellIndex=0;
		for (CellIndex<dim, IndexTrait::vector> cellIndex : CellIndex<dim, IndexTrait::vector>()) 
		{
			tarch::la::Vector<3, int> cellVectorIndex = cellIndex.get();
			bool isInnerCell = true;
			bool isGhostCell = false;
			for (unsigned int currentDim = 0; currentDim < dim; currentDim++) 
			{
				isInnerCell &= cellVectorIndex[currentDim] >= lowerBoundary[currentDim] + (int)_overlap + 1;
				isInnerCell &= cellVectorIndex[currentDim] <= upperBoundary[currentDim] - (int)_overlap - 1 ;
				isGhostCell |= cellVectorIndex[currentDim] < 1;
				isGhostCell |= cellVectorIndex[currentDim] >= upperBoundary[currentDim] ;				
			}
			if (!isInnerCell && !isGhostCell)
			{
				for (unsigned int currentDim = 0; currentDim < dim; currentDim++) {
					_coords[dim*interfaceCellIndex+currentDim]=mdDomainOffset[currentDim] + cellVectorIndex[currentDim]*_dx - _dx + 0.5*macroscopicCellSize[currentDim];
				}
				interfaceCellIndex++;
			}
		}
	    /*
		std::cout << "[";
		for (unsigned int interfaceCellIndex = 0; interfaceCellIndex < _numberOfCells; interfaceCellIndex++) {
			std::cout << "[" << _coords[dim*interfaceCellIndex] << "," << _coords[dim*interfaceCellIndex+1] << "," << _coords[dim*interfaceCellIndex+2] << "]";
		}
		std::cout << "]" << std::endl;
		*/
		
		int meshID = _interface->getMeshID("mamico-mesh");
		_vertexIDs = new int[_numberOfCells];
		_interface->setMeshVertices(meshID, _numberOfCells, _coords, _vertexIDs); 
		_precice_dt = _interface->initialize();
		
		_velocities = new double[_numberOfCells*dim];
		std::cout << "PreciceSolver constructor called" << std::endl;
    }
    
	void setWallVelocity(const tarch::la::Vector<3,double> wallVelocity) override
	{
		_wallVelocity = wallVelocity;
	}
	
	tarch::la::Vector<3,double> getVelocity(tarch::la::Vector<3,double> pos) const override
	{
		std::cout << "getVelocity:" << pos << std::endl;
		tarch::la::Vector<3,double> vel(0.0);
		unsigned int cellIndex=0;
		while (cellIndex < _numberOfCells && _coords[dim*cellIndex] != pos[0] && _coords[dim*cellIndex+1] != pos[1] && _coords[dim*cellIndex+2] != pos[2])
			cellIndex++;
		if (cellIndex == _numberOfCells)
		{
			std::cout << "ERROR PreciceAdapter::getVelocity(): position " << pos << " not found in coupling mesh" << std::endl;
      		//exit(EXIT_FAILURE);
		} /*else 
		{
			for (unsigned int currentDim = 0; currentDim < dim; currentDim++) {
				vel[dim] = _velocities[dim*cellIndex+currentDim];
			}
			//std::cout << "precice Velocity:" << vel << " at point:" << pos << std::endl;
		}*/
		std::cout << "getVelocity:" << pos << std::endl;	
		return vel;
	}
	
	double getDensity(tarch::la::Vector<3,double> pos) const
	{
		return 1.0;
	}
	
	void advance(double dt) override 
	{
		std::cout << "PreciceSolver::advance called with dt = " << dt << std::endl;
		if (_interface->isCouplingOngoing()) {
			int meshID = _interface->getMeshID("mamico-mesh");
			if (_interface->isReadDataAvailable())
			{
				std::cout << "Reading Velocity from MaMiCo !" << std::endl;
				// velocity from the conitnuum solver
				int dataID = _interface->getDataID("Velocity", meshID);
				_interface->readBlockVectorData(dataID, _numberOfCells, _vertexIDs, _velocities);
				
				/*
				std::cout << "[";
				for (unsigned int interfaceCellIndex = 0; interfaceCellIndex < _numberOfCells; interfaceCellIndex++) {
					std::cout << "[" << _velocities[dim*interfaceCellIndex] << "," << _velocities[dim*interfaceCellIndex+1] << "," << _velocities[dim*interfaceCellIndex+2] << "]";
				}
				std::cout << "]" << std::endl;
				*/
				
			}
			
			// Solving the time step
			// Normally does nothing, everything is done on the MD side
			double* mdVelocities = new double[_numberOfCells*dim];
			for (unsigned int vertexIndex = 0; vertexIndex < _numberOfCells; vertexIndex++)
			{
				mdVelocities[dim*vertexIndex] = 0.;
				mdVelocities[dim*vertexIndex+1] = 0.;
				mdVelocities[dim*vertexIndex+2] = 0.;
			}
			
			if (_interface->isWriteDataRequired(dt))
			{
				std::cout << "Writing MD Velocity from MaMiCo !" << std::endl;
				// Velocit from the md solver
				int dataID = _interface->getDataID("MDVelocity", meshID);
				_interface->writeBlockVectorData(dataID, _numberOfCells, _vertexIDs, mdVelocities);
			}
			delete[] mdVelocities;
			double computed_dt = std::min(_precice_dt, dt);
			_precice_dt = _interface->advance(computed_dt);
    	}
	}
	
	bool receiveMacroscopicQuantityFromMDSolver(tarch::la::Vector<dim, unsigned int> globalCellIndex) override
	{
		const tarch::la::Vector<3, unsigned int> globalNumberMacroscopicCells = CellIndex<dim>().numberCellsInDomain;
		bool recv = true;
		for (unsigned int d = 0; d < dim; d++) {		
		  recv = recv && (globalCellIndex[d] > _overlap) &&
		         (globalCellIndex[d] < globalNumberMacroscopicCells[d] + 1 - _overlap);
		}
		return recv;
	}
	
	bool sendMacroscopicQuantityToMDSolver(tarch::la::Vector<3, unsigned int> globalCellIndex) 
	{
		const tarch::la::Vector<3, unsigned int> globalNumberMacroscopicCells = CellIndex<dim>().numberCellsInDomain;
		bool ghost = false;
		for (unsigned int d = 0; d < 3; d++) {
		  ghost = ghost || (globalCellIndex[d] < 1) ||
		          (globalCellIndex[d] > globalNumberMacroscopicCells[d]);
		}
		return (!ghost) &&
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
	unsigned int _overlap;
	tarch::la::Vector<3, double> _mdDomainOffset;
	tarch::la::Vector<3, double> _macroscopicCellSize;
	
    precice::SolverInterface *_interface=NULL;
	double _precice_dt;
	int *_vertexIDs;
	double *_coords;
	double *_velocities;
};
}
}

#endif
