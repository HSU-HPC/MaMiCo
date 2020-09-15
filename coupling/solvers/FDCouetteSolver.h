#ifndef _MOLECULARDYNAMICS_COUPLING_SOLVERS_FINITEDIFFERENCE_H_
#define _MOLECULARDYNAMICS_COUPLING_SOLVERS_FINITEDIFFERENCE_H_

#include "coupling/solvers/NumericalSolver.h"
#if (COUPLING_MD_PARALLEL==COUPLING_MD_YES)
#include <mpi.h>
#endif

namespace coupling{
	namespace solvers{
		class FiniteDifferenceSolver;
	}
}

class coupling::solvers::FiniteDifferenceSolver: public coupling::solvers::NumericalSolver {
public:
	FiniteDifferenceSolver(
		const double channelheight,
		tarch::la::Vector<3,double> wallVelocity,
		const double kinVisc,
		const double dx,
		const double dt,
		const int plotEveryTimestep,
		const std::string filestem,
		const tarch::la::Vector<3,unsigned int> processes,
		const unsigned int numThreads=1):
		coupling::solvers::NumericalSolver(channelheight, dx, dt, kinVisc,
		plotEveryTimestep, filestem, processes), _omega(dt*kinVisc/(dx*dx)),
		_wallVelocity(wallVelocity)
		{
			// return if required
      if (skipRank()){return;}

      _velold = new double[3*(_domainSizeX+2)*(_domainSizeY+2)*(_domainSizeZ+2)];

      #if defined(_OPENMP)
      omp_set_num_threads(numThreads);
      #endif
      #if (COUPLING_MD_DEBUG==COUPLING_MD_YES)
      std::cout << "Domain size=" << _domainSizeX << "," << _domainSizeY << "," << _domainSizeZ << std::endl;
      std::cout << "tau=" << 1.0/_omega << std::endl;
      std::cout << "wallVelocity=" << _wallVelocity << std::endl;
      for (int z = 0; z < _domainSizeZ+2; z++){ for (int y = 0; y < _domainSizeY+2; y++){ for (int x = 0; x < _domainSizeX+2; x++){
        std::cout << x << "," << y << "," << z << "FLAG=" << _flag[get(x,y,z)] << std::endl;
      }}}
      #endif
      // check pointers
			if ( (!_vel) || (!_density) || (!_flag) ){
        std::cout << "ERROR FiniteDifferenceSolver: nullptr!" << std::endl; exit(EXIT_FAILURE);
      }
      #if (COUPLING_MD_PARALLEL==COUPLING_MD_YES)
      if ( (_sendBufferX==nullptr) || (_recvBufferX==nullptr) || (_sendBufferY==nullptr) || (_recvBufferY==nullptr) || (_sendBufferZ==nullptr) || (_recvBufferZ==nullptr) ){
        std::cout << "ERROR FiniteDifferenceSolver: nullptr in send/recv!" << std::endl; exit(EXIT_FAILURE);
      }
      #endif
    }

	~FiniteDifferenceSolver(){
		if (_vel){delete [] _vel; _vel=nullptr;}
		if (_density){delete [] _density; _density=nullptr;}
		if (_flag){delete [] _flag; _flag=nullptr;}
		if (_velold){delete [] _flag; _flag=nullptr;}
		#if (COUPLING_MD_PARALLEL==COUPLING_MD_YES)
		if (_sendBufferX){delete [] _sendBufferX; _sendBufferX=nullptr;}
		if (_sendBufferY){delete [] _sendBufferY; _sendBufferY=nullptr;}
		if (_sendBufferZ){delete [] _sendBufferZ; _sendBufferZ=nullptr;}
		if (_recvBufferX){delete [] _recvBufferX; _recvBufferX=nullptr;}
		if (_recvBufferY){delete [] _recvBufferY; _recvBufferY=nullptr;}
		if (_recvBufferZ){delete [] _recvBufferZ; _recvBufferZ=nullptr;}
		#endif
	}

	void advance (double dt) override {
		const int timesteps=floor( dt/_dt+0.5 );
		if ( fabs(timesteps*_dt-dt)/_dt > 1.0e-8 ){std::cout << "ERROR FiniteDifferenceSolver::advance(): time steps and dt do not match!" << std::endl; exit(EXIT_FAILURE);}
		for (int i = 0; i < timesteps; i++){
			setBeyondWall();
			stream();
			plot();
			//plottxt(); // Use this function as output to analyze the result
			_time += _dt;
			_counter++;}
	}

  // get current velocity for a certain position (vector pos)
  tarch::la::Vector<3,double> getVelocity (tarch::la::Vector<3,double> pos)const override{
		const tarch::la::Vector<3,double> domainOffset(_coords[0]*_dx*_avgDomainSizeX, _coords[1]*_dx*_avgDomainSizeY, _coords[2]*_dx*_avgDomainSizeZ);
		#if (COUPLING_MD_DEBUG==COUPLING_MD_YES)
		if (   (pos[0]<domainOffset[0]) || (pos[0]>domainOffset[0]+_domainSizeX*_dx)
				|| (pos[1]<domainOffset[1]) || (pos[1]>domainOffset[1]+_domainSizeY*_dx)
				|| (pos[2]<domainOffset[2]) || (pos[2]>domainOffset[2]+_domainSizeZ*_dx) ){
			std::cout << "ERROR FiniteDifferenceSolver::getVelocity(): Position " << pos << " out of range!" << std::endl; exit(EXIT_FAILURE);
		}
		#endif
		// compute index for respective cell (_dx+... for ghost cells); use coords to store local cell coordinates
		tarch::la::Vector<3,unsigned int> coords;
		for (unsigned int d = 0; d < 3; d++){ coords[d] = (unsigned int) ((_dx+pos[d]-domainOffset[d])/_dx);}
		const int index = get(coords[0],coords[1],coords[2]);
		tarch::la::Vector<3,double> vel(0.0);
		// extract velocity
		for (int d = 0; d < 3; d++){ vel[d] = _vel[3*index+d]; }
		return vel;
	}

	virtual void setWallVelocity(const tarch::la::Vector<3,double> wallVelocity) override{
			_wallVelocity = wallVelocity;
		}

  // get current density for a certain position
  double getDensity(tarch::la::Vector<3,double> pos) const override{
		tarch::la::Vector<3,unsigned int> coords;
		const tarch::la::Vector<3,double> domainOffset(_coords[0]*_dx*_avgDomainSizeX,_coords[1]*_dx*_avgDomainSizeY,_coords[2]*_dx*_avgDomainSizeZ);
		#if (COUPLING_MD_DEBUG==COUPLING_MD_YES)
		if (   (pos[0]<domainOffset[0]) || (pos[0]>domainOffset[0]+_domainSizeX*_dx)
			 || (pos[1]<domainOffset[1]) || (pos[1]>domainOffset[1]+_domainSizeY*_dx)
			 || (pos[2]<domainOffset[2]) || (pos[2]>domainOffset[2]+_domainSizeZ*_dx) ){
		 std::cout << "ERROR FiniteDifferenceSolver::getDensity(): Position " << pos << " out of range!" << std::endl; exit(EXIT_FAILURE);
		}
		#endif
		// compute index for respective cell (_dx+... for ghost cells); use coords to store local cell coordinates
		for (unsigned int d = 0; d < 3; d++){ coords[d] = (unsigned int) ((_dx+pos[d]-domainOffset[d])/_dx);}
		const int index = get(coords[0],coords[1],coords[2]);
		return _density[index];
	}

  // set values for MD within the solver as boundary conditions for next timestep
  void setMDBoundaryValues(
		std::vector<coupling::datastructures::MacroscopicCell<3>* >& recvBuffer,const unsigned int * const recvIndices,
		const coupling::IndexConversion<3>& indexConversion)override{
		if (skipRank()){return ;}
		for (unsigned int i = 0; i < recvBuffer.size(); i++){
			// determine cell index of this cell in FD domain
			tarch::la::Vector<3,unsigned int> globalCellCoords = indexConversion.getGlobalVectorCellIndex(recvIndices[i]);
			globalCellCoords[0] = (globalCellCoords[0]+_offset[0]) - _coords[0]*_avgDomainSizeX;
			globalCellCoords[1] = (globalCellCoords[1]+_offset[1]) - _coords[1]*_avgDomainSizeY;
			globalCellCoords[2] = (globalCellCoords[2]+_offset[2]) - _coords[2]*_avgDomainSizeZ;
			const int index = get(globalCellCoords[0],globalCellCoords[1],globalCellCoords[2]);
			#if (COUPLING_MD_DEBUG==COUPLING_MD_YES)
			if (_flag[index]!=MD_BOUNDARY){std::cout << "ERROR FiniteDifferenceSolver::setMDBoundaryValues(): Cell " << index << " is no MD boundary cell!" << std::endl; exit(EXIT_FAILURE);}
			#endif
			tarch::la::Vector<3,double> localVel( (1.0/recvBuffer[i]->getMacroscopicMass())*recvBuffer[i]->getMacroscopicMomentum());
				for (unsigned int d = 0; d<3; d++){ _vel[3*index+d] = localVel[d]; }
		}
	}

private:
  // set correct values at the walls, this is necessary since values will be used at cell centers but the wall is
  // at the cell boundary
	void setBeyondWall(){
		for(int i=0; i<_yO; i++){
			_vel[3*i] = 2*_wallVelocity[0] - _vel[3*(i+_yO)];
		}
	}

	// calcuation of velocity field for next timestap
	void stream(){
		double* swap = _velold;
		_velold = _vel;
		_vel = swap;
		for(int i=_yO; i<(_domainSizeX+2)*(_domainSizeY+2)*(_domainSizeZ+2)-_yO; i++){
			if (_flag[i] == FLUID){
				_vel[3*i] = _omega*(_velold[3*(i-_yO)]-2*_velold[3*i]+_velold[3*(i+_yO)])+_velold[3*i];
			}
		}
	}

	const double _omega; // dt*kinVisc/dx/dx - factor for streaming
	tarch::la::Vector<3,double> _wallVelocity; // velocity of moving wall of Couette flow
	double *_velold{nullptr}; // field to store velocities from old timesteps
};
#endif
