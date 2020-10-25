#ifndef _MOLECULARDYNAMICS_COUPLING_SOLVERS_FINITEDIFFERENCE_H_
#define _MOLECULARDYNAMICS_COUPLING_SOLVERS_FINITEDIFFERENCE_H_

#include "tarch/la/Vector.h"
#include <cmath>
#if (COUPLING_MD_PARALLEL==COUPLING_MD_YES)
#include <mpi.h>
#endif
#include "coupling/solvers/CouetteSolver.h"

namespace coupling{
	solvers{
		class FiniteDifferenceSolver;
	}}

class coupling::solvers::FiniteDifferenceSolver: public coupling::solvers::AbstractCouetteSolver<3> {
private:
	enum Flag{FLUID=0,NO_SLIP=1,MOVING_WALL=2,PERIODIC_X_LEFT=3,PERIODIC_X_RIGHT=3,PERIODIC_Y_FRONT=3,PERIODIC_Y_BACK=4,PARALLEL_BOUNDARY=5};
	enum NbFlag{LEFT=0,RIGHT=1,BACK=2,FRONT=3,BOTTOM=4,TOP=5};

public:
	FiniteDifferenceSolver(const double channelheight, const tarch::la::Vector<3,double> wallVelocity, const double kinVisc,
		const double dx, const double dt, const int plotEveryTimestep, const std::string filestem, const tarch::la::Vector<3,unsigned int> processes,
		const unsigned int numThreads=1):
		coupling::solvers::AbstractCouetteSolver(), _wallVelocity(wallVelocity), _time(0.0), _counter(0), _omega(dt*kinVisc/(dx*dx)),
		_dx(dx), _dt(dt), _plotEveryTimestep(plotEveryTimestep), _processes(processes), _coords(getProcessCoordinates()), _totalcellnumber(NULL),
		_vel(NULL), _density(NULL), _flag(NULL), _parallelNeighbours(-1), _offset(-1),  _cells{NULL},
		_domain({getLowerCorner(channelheight, 0), getLowerCorner(channelheight, 0), getLowerCorner(channelheight, 0), getUpperCorner(channelheight, 0),
		getUpperCorner(channelheight, 1), getUpperCorner(channelheight, 2)})
		/*#if (COUPLING_MD_PARALLEL==COUPLING_MD_YES)
    _sendBufferX(NULL), _recvBufferX(NULL), _sendBufferY(NULL), _recvBufferY(NULL), _sendBufferZ(NULL), _recvBufferZ(NULL)
    #endif*/
		{
			// determine parallel neighbours
			determineParallelNeighbours(); // cells are changed depending on the amount of ghost cells
			_totalcellnumber = _cells[0]*_cells[1]*_cells[2];// needs to be executed after determineParallelNeighbours

			_vel = new double[ 3*_totalcellnumber]; //add ghost cells for calculation
			_density = new double[_totalcellnumber];
			_flag = new Flag[_totalcellnumber];
			/*#if (COUPLING_MD_PARALLEL==COUPLING_MD_YES)
			_sendBufferX = new double [5*(_domainSizeY+2)*(_domainSizeZ+2)];
			_recvBufferX = new double [5*(_domainSizeY+2)*(_domainSizeZ+2)];
			_sendBufferY = new double [5*(_domainSizeX+2)*(_domainSizeZ+2)];
			_recvBufferY = new double [5*(_domainSizeX+2)*(_domainSizeZ+2)];
			_sendBufferZ = new double [5*(_domainSizeX+2)*(_domainSizeY+2)];
			_recvBufferZ = new double [5*(_domainSizeX+2)*(_domainSizeY+2)];
			#endif */

			// check pointers
			if ( (_vel==NULL) || (_density==NULL) || (_flag==NULL) ){
				std::cout << "ERROR FiniteDifferenceSolver: NULL ptr!" << std::endl; exit(EXIT_FAILURE);
			}

			/*#if (COUPLING_MD_PARALLEL==COUPLING_MD_YES)
			if ( (_sendBufferX==NULL) || (_recvBufferX==NULL) || (_sendBufferY==NULL) || (_recvBufferY==NULL) || (_sendBufferZ==NULL) || (_recvBufferZ==NULL) ){
				std::cout << "ERROR FiniteDifferenceSolver: NULL ptr in send/recv!" << std::endl; exit(EXIT_FAILURE);
			}
			#endif*/

			// return if required
			if (skipRank()){return;}

			// init everything with zero velocity, unit density; flags are set to FLUID
			for (int i = 0; i < _totalcellnumber; i++){
				const int index1 = get(i); // added, before get(i), one function + saving should be faster
				for (int d = 0; d <  3; d++){ _vel[index1*3+d]  = 0.0;}
				_density[index1] = 1.0;
				_flag[index1] = FLUID;
			} // one time appling get(i) for checking should be enough within the constructor

			// correct boundary flags based on physical description (Couette scenario)
			// bottom - moving wall
			for (int i = 0; i < _cells[0]*_cells[1]; i++){ _flag[get(i)] = MOVING_WALL; } //TODO remove all get(i) within the constructor following this
			//TODO set wall velocities here??
			// top - noslip
			for (int i = (_cells[0]*_cells[1]*(_cells[2]-1); i < _totalcellnumber; i++){ _flag[get(i)] = NO_SLIP;} // before initialised with dsx+2*dsy+2*dsz+1
			// front - periodic
			for (int z = 1; z < _cells[2]-1; z++){for (int x = 0; x < _cells[0]; x++){ _flag[get(x,0,z)] = PERIODIC_Y_FRONT; } }
			// back - periodic
			for (int z = 1; z < _cells[2]-1; z++){for (int x = 0; x < _cells[0]; x++){ _flag[get(x,_cells[1]-1,z)] = PERIODIC_Y_BACK; } }
			// left - periodic
			for (int z = 1; z < _cells[2]-1; z++){for (int y = 1; y < _cells[1]-1; y++){ _flag[get(0,y,z)] = PERIODIC_X_LEFT; }}
			for (int z = 1; z < _cells[2]-1; z++){for (int y = 1; y < _cells[1]-1; y++){ _flag[get(_cells[0]-1,y,z)] = PERIODIC_X_RIGHT; }}

			// correct boundary flags in case of MPI-parallel simulations (Couette scenario)
			//setParallelBoundaryFlags();

			//check if all border cells non Fluid
			for(i=0; i<xO; i++){
				if(_flag[i]==FLUID){std::cout << "ERROR FiniteDifferenceSolver::initialisation: setup of domain is wrong, x-border cell has flag FLUID" << std::endl; exit(EXIT_FAILURE);}}
			//TODO: same for all other planes
		}

	~FiniteDifferenceSolver(){
		if (_vel !=NULL){delete [] _vel; _vel=NULL;}
		if (_density!=NULL){delete [] _density; _density=NULL;}
		if (_flag!=NULL){delete [] _flag; _flag=NULL;}
		/*#if (COUPLING_MD_PARALLEL==COUPLING_MD_YES)
		if (_sendBufferX!=NULL){delete [] _sendBufferX; _sendBufferX=NULL;}
		if (_sendBufferY!=NULL){delete [] _sendBufferY; _sendBufferY=NULL;}
		if (_sendBufferZ!=NULL){delete [] _sendBufferZ; _sendBufferZ=NULL;}
		if (_recvBufferX!=NULL){delete [] _recvBufferX; _recvBufferX=NULL;}
		if (_recvBufferY!=NULL){delete [] _recvBufferY; _recvBufferY=NULL;}
		if (_recvBufferZ!=NULL){delete [] _recvBufferZ; _recvBufferZ=NULL;}
		#endif*/
	}

	void advance override (double dt){ // changed: rather override than virtual
		const int timesteps=floor( dt/_dt+0.5 );
		if ( fabs(timesteps*_dt-dt)/_dt > 1.0e-8 ){std::cout << "ERROR FiniteDifferenceSolver::advance(): time steps and dt do not match!" << std::endl; exit(EXIT_FAILURE);}
		for (int i = 0; i < timesteps; i++){
			plot();
			stream();
			communicate(); // exchange between neighbouring MPI subdomains
			_time += _dt;
			_counter++;}
	}

	tarch::la::Vector<3,double> getVelocity override (tarch::la::Vector<3,double> pos){ // same as with the function above
		const tarch::la::Vector<3,double> domainOffset(_coords[0]*_dx*_avgDomainSizeX, _coords[1]*_dx*_avgDomainSizeY, _coords[2]*_dx*_avgDomainSizeZ);
		// check pos-data for process locality (todo: put this in debug mode in future releases)
		if (   (pos[0]<domainOffset[0]) || (pos[0]>domainOffset[0]+_domainSizeX*_dx)
				|| (pos[1]<domainOffset[1]) || (pos[1]>domainOffset[1]+_domainSizeY*_dx)
				|| (pos[2]<domainOffset[2]) || (pos[2]>domainOffset[2]+_domainSizeZ*_dx) ){
			std::cout << "ERROR FiniteDifferenceSolver::getVelocity(): Position " << pos << " out of range!" << std::endl; exit(EXIT_FAILURE);
		}
		// compute index for respective cell (_dx+... for ghost cells); use coords to store local cell coordinates
		for (unsigned int d = 0; d < 3; d++){ localcoords[d] = (unsigned int) ((_dx+pos[d]-domainOffset[d])/_dx);}
		const int index = get(localcoords[0],localcoords[1],localcoords[2]);
		tarch::la::Vector<3,double> vel(0.0);
		// extract and scale velocity to "real"=MD units
		for (int d = 0; d < 3; d++){ vel[d] = _dx/_dt*_vel[3*index+d]; }
		return vel;
	}

private:
	/** returns i and performs checks in debug mode */
	int get(int i) const {return i;}

	/** returns linearized index and performs checks in debug mode */
	int get(int x,int y,int z) const{return x + (_domainSizeX+2)*(y+(_domainSizeY+2)*z);}

	// Get upper, right, back corner
	int getUpperCorner(double channelheight, int d) const {
		#if (COUPLING_MD_PARALLEL==COUPLING_MD_YES)
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD,&rank);
		// if this rank is outside the range given by _processes: return 0 // -> cannot use method skipRank() here, since __processes may not be initialized yet!
		if (skipRank(){return 0;}
		#endif
		// if this is not the last process along this direction: just return avg. number of cells
		if (_coords[d]<_processes[d]-1){
			_cells[d] = floor( 1+(channelheight+0.5)/dx )/_processes[d];
			return _domain[d-3]+floor(channelheight/_processes[d]);}
		// otherwise: add the cells that have not been distributed yet
		else {
			int globalDomainSize = floor( 1+(channelheight+0.5)/dx );
			_cells[d] = globalDomainSize/_processes[d] + globalDomainSize%_processes[d];
			return _domain[d-3]+floor(channelheight/_processes[d]) + channelheight%_processes[d];}
	}

	// get lower, left, front corner
	int getLowerCorner(double channelheight, int d){
		if(_processes[d]==1){return 0}
		#if (COUPLING_MD_PARALLEL==COUPLING_MD_YES)
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD,&rank);
		// if this rank is outside the range given by _processes: return 0
		if (skipRank()){return 0;}
		else if(_processes[d]){
			if(rank == 0){ return 0}
			else{
				return (rank-1)*floor(channelheight/_processes[d])
		}}
		#endif
		// TODO make else producing error and telling something is wrong
	}

	/** determines the "avg" domain size which is the domain size on each MPI process, except for potentially the last one (the last one may include additional cells) */
	int getAvgDomainSize(double channelheight,double dx, tarch::la::Vector<3,unsigned int> processes, int d) const {
		int globalDomainSize = floor( (channelheight+0.5)/dx );
		return globalDomainSize/processes[d];
	}

	tarch::la::Vector<3,unsigned int> getProcessCoordinates() const{
		tarch::la::Vector<3,unsigned int> coords(0);
		/*#if (COUPLING_MD_PARALLEL==COUPLING_MD_YES)
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD,&rank);
		// determine rank coordinates
		coords[2] = ((unsigned int)rank)/(_processes[0]*_processes[1]);
		coords[1] =(((unsigned int)rank)-coords[2]*_processes[0]*_processes[1])/_processes[0];
		coords[0] =((unsigned int)rank)-coords[2]*_processes[0]*_processes[1]-coords[1]*_processes[0];
		#endif*/
		return coords;

	/** create vtk plot if required */ // Copied from LBCouetteSolver & changed in parts
	void plot()  {
		// only plot output if this is the correct timestep
		if (_plotEveryTimestep==-1){ return;}
		if (_counter%_plotEveryTimestep!=0){return;}
		int rank = 0; // rank in MPI-parallel simulations
		#if (COUPLING_MD_PARALLEL==COUPLING_MD_YES)
		MPI_Comm_rank(MPI_COMM_WORLD,&rank);
		#endif

		std::stringstream ss; ss << _filestem << "_" << rank << "_" << _counter << ".vtk";
		std::ofstream file(ss.str().c_str());
		if (!file.is_open()){std::cout << "ERROR LBCouetteSolver::plot(): Could not open file " << ss.str() << "!" << std::endl; exit(EXIT_FAILURE);}
		std::stringstream flag, density, velocity;

		file << "# vtk DataFile Version 2.0" << std::endl;
		file << "MaMiCo FiniteDifferenceSolver" << std::endl;
		file << "ASCII" << std::endl << std::endl;
		file << "DATASET STRUCTURED_GRID" << std::endl;
		file << "DIMENSIONS " << _domainSizeX+1 << " " << _domainSizeY+1 << " " << _domainSizeZ+1 << std::endl;
		file << "POINTS " << (_domainSizeX+1)*(_domainSizeY+1)*(_domainSizeZ+1) << " float" << std::endl;

		flag << "CELL_DATA " << _domainSizeX*_domainSizeY*_domainSizeZ << std::endl;
		flag << "SCALARS flag float 1" << std::endl;
		flag << "LOOKUP_TABLE default" << std::endl;

		density  << std::setprecision(12);
		density  << "SCALARS density float 1 " << std::endl;
		density  << "LOOKUP_TABLE default" << std::endl;

		velocity << std::setprecision(12);
		velocity << "VECTORS velocity float" << std::endl;
		// loop over domain (incl. boundary) and write point coordinates
		for (int z = 0; z < _domainSizeZ+1; z++){ for (int y = 0; y < _domainSizeY+1; y++){ for (int x = 0; x < _domainSizeX+1; x++){
			file << (_coords[0]*_avgDomainSizeX+x)*_dx << " " <<  (_coords[1]*_avgDomainSizeY+y)*_dx  << " " << (_coords[2]*_avgDomainSizeZ+z)*_dx << std::endl;
		}}}
		// loop over domain (incl. boundary)
		for (int z = 1; z < _domainSizeZ+1; z++){ for (int y = 1; y < _domainSizeY+1; y++){ for (int x = 1; x < _domainSizeX+1; x++){
			const int index=get(x,y,z);
			// write information to streams
			flag << _flag[index] << std::endl;
			density << _density[index] << std::endl;
			velocity << _vel[3*index] << " " << _vel[3*index+1] << " " << _vel[3*index+2] << std::endl;
		}}}

		file << std::endl;
		file << flag.str() << std::endl << std::endl;
		flag.str("");
		file << density.str() << std::endl;
		density.str("");
		file << velocity.str() << std::endl;
		velocity.str("");
		file.close();
	}

	// calcuation of velocity field for next timestap
	void stream(){
		double velold[] = _vel;
		for(i=0; i<_totalcellnumber; i=i+3){ // TODO: change the implementation, there shouldn't be any if in here, rather a switch
			if (_flag[i] == FLUID){
				_vel[i] += _omega*(velold[i-3]-2*velold[i]+velold[i+3]);
				}
			else if(_flag[i] == PERIODIC_X_LEFT){
				_vel[i] += _omega*(velold[i+3*(_cells[0]-1)]-2*velold[i]+velold[i+3]);
				//_vel[i+1] += _omega*(velold[])
				//_vel[i+2] +=
			}
			else if(_flag[i] == PERIODIC_X_RIGHT){

			}
			else if(_flag[i] == PERIODIC_Y_FRONT){

			}
			else if(_flag[i] == PERIODIC_Y_BACK){

			}
			else if(_flag[i] == NO_SLIP){
				_vel[i] = 0;
				_vel[i+1] = 0;
				_vel[i+2] = 0;
			}
			else if(_flag[i] == MOVING_WALL){
				_vel[i] = _wallVelocity[0];
				_vel[i+1] = _wallVelocity[1];
				_vel[i+2] = _wallVelocity[2];
			} //Maybe set NULL to recognize when something goes wrong
			else{}
		}
	}

	//comunication with other domain (get values for ghost cells)
	void communicate(){

	}

	bool skipRank() const {
		int rank = 0;
		/*#if (COUPLING_MD_PARALLEL==COUPLING_MD_YES)
		MPI_Comm_rank(MPI_COMM_WORLD,&rank);
		#endif*/
		return ((unsigned int)rank>_processes[0]*_processes[1]*_processes[2]-1);
	}

	void determineParallelNeighbours(){
		#if (COUPLING_MD_PARALLEL==COUPLING_MD_YES)
		int rank;
		int size;
		MPI_Comm_rank(MPI_COMM_WORLD,&rank);
		MPI_Comm_size(MPI_COMM_WORLD,&size);
		// check if enough ranks are available
		if (_processes[0]*_processes[1]*_processes[2]>(unsigned int)size){
			std::cout << "ERROR FiniteDifferenceSolver::determineParallelNeighbours(): Not enough ranks available!" << std::endl; exit(EXIT_FAILURE);
		}
		// neighbour dependencies based on Couette problem
		// left,right: periodic
		if(processes[0]!=1){
		_parallelNeighbours[LEFT] = ((_coords[0]+_processes[0]-1) % _processes[0]) + _processes[0]*(_coords[1]+_processes[1]*_coords[2]);
		_parallelNeighbours[RIGHT] = ((_coords[0]+_processes[0]+1) % _processes[0]) + _processes[0]*(_coords[1]+_processes[1]*_coords[2]);
		cells[0] += 2;}
		// back,front: periodic
		if(processes[1]!=1){
		_parallelNeighbours[FRONT]= _coords[0] + _processes[0]*( ((_coords[1]+_processes[1]-1) % _processes[1]) + _processes[1]*_coords[2]);
		_parallelNeighbours[BACK] = _coords[0] + _processes[0]*( ((_coords[1]+_processes[1]+1) % _processes[1]) + _processes[1]*_coords[2]);
		cells[1] += 2;}
		// top: either neighbour or MPI_PROC_NULL
		if(processes[2]!=1){}
		if (_coords[2]==_processes[2]-1){ _parallelNeighbours[TOP] = MPI_PROC_NULL; } else { _parallelNeighbours[TOP] = _coords[0]+_processes[0]*(_coords[1]+_processes[1]*(_coords[2]+1)); }
		// bottom either neighbour or MPI_PROC_NULL
		if (_coords[2]==0)              { _parallelNeighbours[BOTTOM]=MPI_PROC_NULL;} else { _parallelNeighbours[BOTTOM]=_coords[0]+_processes[0]*(_coords[1]+_processes[1]*(_coords[2]-1));}
		cells[2] += 2;}
		std::cout << "Parallel neighbours for rank " << rank << ": " << _parallelNeighbours << std::endl;
		#endif
	}

	/** sets parallel boundary flags according to Couette scenario */
	void setParallelBoundaryFlags(){
		#if (COUPLING_MD_PARALLEL==COUPLING_MD_YES)
		const tarch::la::Vector<3,unsigned int> coords(getProcessCoordinates());
		// bottom - moving wall
		if (coords[2]!=0){
			for (int i = 0; i < (_domainSizeX+2)*(_domainSizeY+2); i++){ _flag[get(i)] = PARALLEL_BOUNDARY; }
		}
		// top - noslip
		if (coords[2]!=_processes[2]-1){
			for (int i = (_domainSizeX+2)*(_domainSizeY+2)*(_domainSizeZ+1); i < (_domainSizeX+2)*(_domainSizeY+2)*(_domainSizeZ+2); i++){ _flag[get(i)] = PARALLEL_BOUNDARY;}
		}
		// for all other boundaries front,back,left,right, we use parallel boundaries. If we send from one processes to itself, this is still the same as periodic conditions
		// front - periodic
		for (int z = 1; z < _domainSizeZ+1; z++){for (int x = 0; x < _domainSizeX+2; x++){ _flag[get(x,0,z)] = PARALLEL_BOUNDARY; } }
		// back - periodic
		for (int z = 1; z < _domainSizeZ+1; z++){for (int x = 0; x < _domainSizeX+2; x++){ _flag[get(x,_domainSizeY+1,z)] = PARALLEL_BOUNDARY; } }
		// left - periodic
		for (int z = 1; z < _domainSizeZ+1; z++){for (int y = 1; y < _domainSizeY+1; y++){ _flag[get(0,y,z)] = PARALLEL_BOUNDARY; }}
		for (int z = 1; z < _domainSizeZ+1; z++){for (int y = 1; y < _domainSizeY+1; y++){ _flag[get(_domainSizeX+1,y,z)] = PARALLEL_BOUNDARY; }}
		#endif // COUPLING_MD_PARALLEL
	}


	const double _omega; // relaxation frequency
	const tarch::la::Vector<3,double> _wallVelocity; // velocity of moving wall of Couette flow
	const int _domainSizeX; // domain size in x-direction
	const int _domainSizeY; // domain size in y-direction
	const int _domainSizeZ; // domain size in z-direction
	const int _avgDomainSizeX; // avg. domain size in MPI-parallel simulation in x-direction
	const int _avgDomainSizeY; // "" in y-direction
	const int _avgDomainSizeZ; // "" in z-direction
	double _time; // simulation time
	int _counter; // time step counter
	const double _dx; // actual mesh size for x, y, z
	const double _dt; // actual time step size
	const int _plotEveryTimestep; // number of time steps between vtk plots
	const std::string _filestem; // file stem for vtk plot
	double *_vel;  // velocity field
	double *_density; // density
	Flag *_flag; // flag field // TODO could be const?!
	#if (COUPLING_MD_PARALLEL==COUPLING_MD_YES)
	double *_sendBufferX; // buffer to send data from left/right to right/left neighbour
	double *_recvBufferX;
	double *_sendBufferY; // buffer to receive data from from left/right neighbour
	double *_recvBufferY;
	double *_sendBufferZ;
	double *_recvBufferZ;
	#endif
	tarch::la::Vector<3,unsigned int> _processes; // domain decomposition on MPI rank basis // TODO could be const?
	tarch::la::Vector<3,unsigned int> _coords; // TODO could be const
	tarch::la::Vector<6,int> _parallelNeighbours; // neighbour ranks // TODO could be const
	const int _totalcellnumber; //absolut number of cells with ghost cells
	const tarch::la::Vector<3,double> _domainOffset; //position of domain for calculation (without ghost cells) // TODO do I need a ghost cell information?
	// for coupling with MD
	tarch::la::Vector<3,int> _offset;
	tarch::la::Vector<3,int> _globalNumberMacroscopicCells;
}

#endif
