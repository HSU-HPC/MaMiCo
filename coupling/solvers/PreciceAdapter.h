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
    	int rank = 0;
    	int size = 1;
    	#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
			MPI_Comm_rank(MPI_COMM_WORLD, &rank);
			MPI_Comm_size(MPI_COMM_WORLD, &size);
		#endif
    	_interface = new precice::SolverInterface("mamico","../precice-config.xml",rank,size);
		_numberOfCells = 0;
		for (CellIndex<dim, IndexTrait::vector> cellIndex : CellIndex<dim, IndexTrait::vector>()) 
		{
			if (sendMacroscopicQuantityToMDSolver(static_cast<tarch::la::Vector<dim, unsigned int>>(cellIndex.get())))
				_numberOfCells++;
		}
		std::cout << "Number of coupling volumes: " << _numberOfCells << std::endl;
		_coords = new double[dim*_numberOfCells];
		unsigned int interfaceCellIndex=0;
		for (CellIndex<dim, IndexTrait::vector> cellIndex : CellIndex<dim, IndexTrait::vector>()) 
		{
			tarch::la::Vector<3, unsigned int> cellVectorIndex = static_cast<tarch::la::Vector<dim, unsigned int>>(cellIndex.get());
			if (sendMacroscopicQuantityToMDSolver(cellVectorIndex))
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
    }
    
	void setWallVelocity(const tarch::la::Vector<3,double> wallVelocity) override
	{
		_wallVelocity = wallVelocity;
	}
	
	tarch::la::Vector<3,double> getVelocity(tarch::la::Vector<3,double> pos) const override
	{
		tarch::la::Vector<3,double> vel(0.0);
		unsigned int cellIndex=0;
		while (cellIndex < _numberOfCells && !(_coords[dim*cellIndex] == pos[0] && _coords[dim*cellIndex+1] == pos[1] && _coords[dim*cellIndex+2] == pos[2]))
			cellIndex++;
		if (cellIndex < _numberOfCells)
		{
			for (unsigned int currentDim = 0; currentDim < dim; currentDim++) {
				vel[currentDim] = _velocities[dim*cellIndex+currentDim];
			}
		}
		std::cout << "pos:" << pos << ", vel:" << vel << std::endl;
		return vel;
	}
	
	double getDensity(tarch::la::Vector<3,double> pos) const
	{
		return 1.0;
	}
	
	void advance(double dt) override 
	{
		if (_interface->isCouplingOngoing()) {
			int meshID = _interface->getMeshID("mamico-mesh");
			if (_interface->isReadDataAvailable())
			{
				std::cout << "Reading external velocities in MaMiCo !" << std::endl;
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
				/*std::cout << "Writing MD Velocity from MaMiCo !" << std::endl;
				// Velocit from the md solver
				int dataID = _interface->getDataID("MDVelocity", meshID);
				_interface->writeBlockVectorData(dataID, _numberOfCells, _vertexIDs, mdVelocities);
				*/
			}
			delete[] mdVelocities;
			double computed_dt = std::min(_precice_dt, dt);
			_precice_dt = _interface->advance(computed_dt);
    	}
	}
	
	bool receiveMacroscopicQuantityFromMDSolver(tarch::la::Vector<dim, unsigned int> globalCellIndex) override
	{
		tarch::la::Vector<3, int> lowerBoundary = CellIndex<3>::lowerBoundary.get();
		tarch::la::Vector<3, int> upperBoundary = CellIndex<3>::upperBoundary.get();
		bool isInnerCell = true;
		for (unsigned int currentDim = 0; currentDim < dim; currentDim++) 
		{
			isInnerCell &= (int)globalCellIndex[currentDim] >= lowerBoundary[currentDim] + (int)_overlap + 1;
			isInnerCell &= (int)globalCellIndex[currentDim] <= upperBoundary[currentDim] - (int)_overlap - 1 ;
		}
		return isInnerCell;
	}
	
	bool sendMacroscopicQuantityToMDSolver(tarch::la::Vector<3, unsigned int> globalCellIndex) 
	{
		tarch::la::Vector<3, int> upperBoundary = CellIndex<3>::upperBoundary.get();
		bool isGhostCell = false;
		for (unsigned int currentDim = 0; currentDim < dim; currentDim++) 
		{
			isGhostCell |= globalCellIndex[currentDim] < 1;
			isGhostCell |= (int)globalCellIndex[currentDim] >= upperBoundary[currentDim] ;
		}
		return (!isGhostCell) &&
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
