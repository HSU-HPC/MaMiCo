// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_SOLVERS_NUMERICALSOLVER_H_
#define _MOLECULARDYNAMICS_COUPLING_SOLVERS_NUMERICALSOLVER_H_

#include "coupling/CouplingMDDefinitions.h"
#include "tarch/la/Vector.h"
#include <cmath>
#include <sstream>
#include <fstream>
#include <iomanip>
#if (COUPLING_MD_PARALLEL==COUPLING_MD_YES)
#include <mpi.h>
#endif
#include "coupling/solvers/CouetteSolver.h"
#include "coupling/datastructures/MacroscopicCell.h"
#include "coupling/IndexConversion.h"

namespace coupling {
  namespace solvers{
    class NumericalSolver;
  }
}

/** It implements the basic functions for a numerical Couette flow solver.
  * For impelmenting a Couette flow solver, you can use it in the hierachy.
  *  @author Philipp Neumann & Helene Wittenberg
 */
class coupling::solvers::NumericalSolver: public coupling::solvers::AbstractCouetteSolver<3> {
  public:
    NumericalSolver(
      const double channelheight,
      const double dx,
      const double dt,
      const double kinVisc,
      const int plotEveryTimestep,
      const std::string filestem,
      const tarch::la::Vector<3,unsigned int> processes) :
    coupling::solvers::AbstractCouetteSolver<3>(), _channelheight(channelheight),
    _dx(dx), _dt(dt), _kinVisc(kinVisc), _processes(processes),
    _plotEveryTimestep(plotEveryTimestep), _filestem(filestem)
    {
      _vel = new double[ 3*(_domainSizeX+2)*(_domainSizeY+2)*(_domainSizeZ+2)];
      _density = new double[(_domainSizeX+2)*(_domainSizeY+2)*(_domainSizeZ+2)];
      _flag = new Flag[(_domainSizeX+2)*(_domainSizeY+2)*(_domainSizeZ+2)];
      #if (COUPLING_MD_PARALLEL==COUPLING_MD_YES)
      _sendBufferX = new double [5*(_domainSizeY+2)*(_domainSizeZ+2)];
      _recvBufferX = new double [5*(_domainSizeY+2)*(_domainSizeZ+2)];
      _sendBufferY = new double [5*(_domainSizeX+2)*(_domainSizeZ+2)];
      _recvBufferY = new double [5*(_domainSizeX+2)*(_domainSizeZ+2)];
      _sendBufferZ = new double [5*(_domainSizeX+2)*(_domainSizeY+2)];
      _recvBufferZ = new double [5*(_domainSizeX+2)*(_domainSizeY+2)];
      #endif
      // zero velocity, unit density; flags are set to FLUID
      for (int i = 0; i < (_domainSizeX+2)*(_domainSizeY+2)*(_domainSizeZ+2); i++){
        for (int d = 0; d <  3; d++){ _vel[i*3+d]  = (double)0.0;}
        _density[i] = 1.0;
        _flag[i] = FLUID;
      }
      // determine parallel neighbours
      determineParallelNeighbours();

      // correct boundary flags based on physical description (Couette scenario)
      // bottom - moving wall
      for (int i = 0; i < (_domainSizeX+2)*(_domainSizeY+2); i++){ _flag[i] = MOVING_WALL; }
      // top - noslip
      for (int i = (_domainSizeX+2)*(_domainSizeY+2)*(_domainSizeZ+1); i < (_domainSizeX+2)*(_domainSizeY+2)*(_domainSizeZ+2); i++){ _flag[i] = NO_SLIP;}
      // front - periodic
      for (int z = 1; z < _domainSizeZ+1; z++){for (int x = 0; x < _domainSizeX+2; x++){ _flag[get(x,0,z)] = PERIODIC; } }
      // back - periodic
      for (int z = 1; z < _domainSizeZ+1; z++){for (int x = 0; x < _domainSizeX+2; x++){ _flag[get(x,_domainSizeY+1,z)] = PERIODIC; } }
      // left - periodic
      for (int z = 1; z < _domainSizeZ+1; z++){for (int y = 1; y < _domainSizeY+1; y++){ _flag[get(0,y,z)] = PERIODIC; }}
      for (int z = 1; z < _domainSizeZ+1; z++){for (int y = 1; y < _domainSizeY+1; y++){ _flag[get(_domainSizeX+1,y,z)] = PERIODIC; }}

      // correct boundary flags in case of MPI-parallel simulations (Couette scenario)
      setParallelBoundaryFlags();
    }

    virtual ~NumericalSolver(){
      if (_vel !=NULL){delete [] _vel; _vel=NULL;}
      if (_density!=NULL){delete [] _density; _density=NULL;}
      if (_flag!=NULL){delete [] _flag; _flag=NULL;}
      #if (COUPLING_MD_PARALLEL==COUPLING_MD_YES)
      if (_sendBufferX!=NULL){delete [] _sendBufferX; _sendBufferX=NULL;}
      if (_sendBufferY!=NULL){delete [] _sendBufferY; _sendBufferY=NULL;}
      if (_sendBufferZ!=NULL){delete [] _sendBufferZ; _sendBufferZ=NULL;}
      if (_recvBufferX!=NULL){delete [] _recvBufferX; _recvBufferX=NULL;}
      if (_recvBufferY!=NULL){delete [] _recvBufferY; _recvBufferY=NULL;}
      if (_recvBufferZ!=NULL){delete [] _recvBufferZ; _recvBufferZ=NULL;}
      #endif
    }

    /** flags the domain boundary cells. mdDomainOffset and mdDomainSize correspond to lower/left/front corner of the MD domain and the size of the domain. overlapStrip is the number of cells
     *  that the LB and MD domain shall overlap, i.e. the number of "MD cells" that lie within the LB domain but which should not be flagged (and thus be handled by both solvers).*/
    void setMDBoundary(tarch::la::Vector<3,double> mdDomainOffset,tarch::la::Vector<3,double> mdDomainSize,unsigned int overlapStrip,
    const coupling::IndexConversion<3>& indexConversion, const unsigned int* const recvIndice, unsigned int size){
      if (skipRank()){return;}

      for (int d = 0; d < 3; d++){
        _offset[d] = (floor(mdDomainOffset[d]/_dx + 0.5));
        if (fabs(_offset[d]*_dx - mdDomainOffset[d])/_dx>1.0e-8 ){std::cout << "ERROR NumericalSolver::setMDBoundary(): offset does not match!" << std::endl; exit(EXIT_FAILURE);}
        _globalNumberMacroscopicCells[d] = (floor(mdDomainSize[d]/_dx + 0.5));
        if (fabs(_globalNumberMacroscopicCells[d]*_dx - mdDomainSize[d])/_dx>1.0e-8 ){std::cout << "ERROR NumericalSolver::setMDBoundary(): globalNumber does not match!" << std::endl; exit(EXIT_FAILURE);}
      }

      // flag local domain
      for (int z = 0; z < _domainSizeZ+2; z++){
        for (int y = 0; y < _domainSizeY+2; y++){
          for (int x = 0; x < _domainSizeX+2; x++){
            // determine global cell coordinates of the process-local (sub)domain
            tarch::la::Vector<3,int> globalCoords(x-1+_coords[0]*_avgDomainSizeX,y-1+_coords[1]*_avgDomainSizeY,z-1+_coords[2]*_avgDomainSizeZ);
            bool isMDCell=true;
            for (int d = 0; d < 3; d++){
              isMDCell = isMDCell && (globalCoords[d]>_offset[d]+(int)overlapStrip-1) && (globalCoords[d]<_offset[d]+_globalNumberMacroscopicCells[d]-(int)overlapStrip);
            }
            if (isMDCell){
              _flag[get(x,y,z)] = MD_BOUNDARY;
            }
          }
        }
      }
    }

    /** Function applies the values received from the MD-solver within the conntinuum solver */
    virtual void setMDBoundaryValues(std::vector<coupling::datastructures::MacroscopicCell<3>* >& recvBuffer,const unsigned int * const recvIndices,
    const coupling::IndexConversion<3>& indexConversion)=0;

    /** getters required by LB Couette Solver Interface */
    tarch::la::Vector<3,unsigned int> getNumberProcesses() const { return _processes;}

    tarch::la::Vector<3,unsigned int> getAvgNumberLBCells() const { tarch::la::Vector<3,unsigned int> avgCells(_avgDomainSizeX,_avgDomainSizeY,_avgDomainSizeZ); return avgCells;}

    virtual double getDensity(tarch::la::Vector<3,double> pos) const = 0;

    virtual void setWallVelocity(const tarch::la::Vector<3,double> wallVelocity) = 0;

  private:
    tarch::la::Vector<3,unsigned int> getProcessCoordinates() const{
      tarch::la::Vector<3,unsigned int> coords(0);
      #if (COUPLING_MD_PARALLEL==COUPLING_MD_YES)
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD,&rank);
      // determine rank coordinates
      coords[2] = ((unsigned int)rank)/(_processes[0]*_processes[1]);
      coords[1] =(((unsigned int)rank)-coords[2]*_processes[0]*_processes[1])/_processes[0];
      coords[0] =((unsigned int)rank)-coords[2]*_processes[0]*_processes[1]-coords[1]*_processes[0];
      #endif
      return coords;
    }

    void determineParallelNeighbours(){
      #if (COUPLING_MD_PARALLEL==COUPLING_MD_YES)
      int rank;
      int size;
      MPI_Comm_rank(MPI_COMM_WORLD,&rank);
      MPI_Comm_size(MPI_COMM_WORLD,&size);
      // check if enough ranks are available
      if (_processes[0]*_processes[1]*_processes[2]>(unsigned int)size){
        std::cout << "ERROR NumericalSolver::determineParallelNeighbours(): Not enough ranks available!" << std::endl; exit(EXIT_FAILURE);
      }
      // neighbour dependencies based on Couette problem
      // left,right: periodic
      _parallelNeighbours[LEFT] = ((_coords[0]+_processes[0]-1) % _processes[0]) + _processes[0]*(_coords[1]+_processes[1]*_coords[2]);
      _parallelNeighbours[RIGHT]= ((_coords[0]+_processes[0]+1) % _processes[0]) + _processes[0]*(_coords[1]+_processes[1]*_coords[2]);
      // back,front: periodic
      _parallelNeighbours[FRONT]= _coords[0] + _processes[0]*( ((_coords[1]+_processes[1]-1) % _processes[1]) + _processes[1]*_coords[2]);
      _parallelNeighbours[BACK] = _coords[0] + _processes[0]*( ((_coords[1]+_processes[1]+1) % _processes[1]) + _processes[1]*_coords[2]);
      // top: either neighbour or MPI_PROC_NULL
      if (_coords[2]==_processes[2]-1){ _parallelNeighbours[TOP] = MPI_PROC_NULL; } else { _parallelNeighbours[TOP] = _coords[0]+_processes[0]*(_coords[1]+_processes[1]*(_coords[2]+1)); }
      // bottom either neighbour or MPI_PROC_NULL
      if (_coords[2]==0)              { _parallelNeighbours[BOTTOM]=MPI_PROC_NULL;} else { _parallelNeighbours[BOTTOM]=_coords[0]+_processes[0]*(_coords[1]+_processes[1]*(_coords[2]-1));}
      // std::cout << "Parallel neighbours for rank " << rank << ": " << _parallelNeighbours << std::endl;
      #endif
    }

    /** sets parallel boundary flags according to Couette scenario */
    void setParallelBoundaryFlags(){
      #if (COUPLING_MD_PARALLEL==COUPLING_MD_YES)
      // bottom - moving wall
      if (_coords[2]!=0){
        for (int i = 0; i < (_domainSizeX+2)*(_domainSizeY+2); i++){ _flag[i] = PARALLEL_BOUNDARY; }
      }
      // top - noslip
      if (_coords[2]!=_processes[2]-1){
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

    /** determines the local domain size on this rank where channelheight is the domain length in direction d. */
    int getDomainSize(double channelheight, double dx, tarch::la::Vector<3,unsigned int> processes, int d) const {
      int globalDomainSize = floor( (channelheight+0.5)/dx );
      tarch::la::Vector<3,unsigned int> coords(0);
      #if (COUPLING_MD_PARALLEL==COUPLING_MD_YES)
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD,&rank);
      // if this rank is outside the range given by processes: return 0
      // -> cannot use method skipRank() here, since _processes may not be initialized yet!
      if ((unsigned int)rank > processes[0]*processes[1]*processes[2]-1){return 0;}
      // determine rank coordinates
      coords[2] = ((unsigned int)rank)/(processes[0]*processes[1]);
      coords[1] =(((unsigned int)rank)-coords[2]*processes[0]*processes[1])/processes[0];
      coords[0] =((unsigned int)rank)-coords[2]*processes[0]*processes[1]-coords[1]*processes[0];
      #endif
      // if this is not the last process along this direction: just return avg. number of cells
      if (coords[d]<processes[d]-1){
        return globalDomainSize/processes[d];
      // otherwise: add the cells that have not been distributed yet
      } else {return globalDomainSize/processes[d] + globalDomainSize%processes[d];}
    }

    /** determines the "avg" domain size which is the domain size on each MPI process, except for potentially the last one (the last one may include additional cells) */
    int getAvgDomainSize(double channelheight,double dx, tarch::la::Vector<3,unsigned int> processes, int d) const {
      int globalDomainSize = floor( (channelheight+0.5)/dx );
      return globalDomainSize/processes[d];
    }

  protected:
    /** returns i and performs checks in debug mode */
    int get(int i) const {
      #if (COUPLING_MD_DEBUG==COUPLING_MD_YES)
      if (i<0){std::cout << "ERROR NumericalSolver::get(i): i<0!" << std::endl; exit(EXIT_FAILURE);}
      if (i>(_domainSizeX+2)*(_domainSizeY+2)*(_domainSizeZ+2)){std::cout << "ERROR NumericalSolver::get(i): i>max. Value!" << std::endl; exit(EXIT_FAILURE);}
      #endif
      return i;
    }

    /** returns linearized index and performs checks in debug mode */
    int get(int x,int y,int z) const{
      #if (COUPLING_MD_DEBUG==COUPLING_MD_YES)
      if (x<0){std::cout << "ERROR NumericalSolver::get(x,y,z): x<0!" << std::endl; exit(EXIT_FAILURE);}
      if (x>_domainSizeX+1){std::cout << "ERROR NumericalSolver::get(x,y,z): x>max. Value!" << std::endl; exit(EXIT_FAILURE);}
      if (y<0){std::cout << "ERROR NumericalSolver::get(x,y,z): y<0!" << std::endl; exit(EXIT_FAILURE);}
      if (y>_domainSizeY+1){std::cout << "ERROR NumericalSolver::get(x,y,z): y>max. Value!" << std::endl; exit(EXIT_FAILURE);}
      if (z<0){std::cout << "ERROR NumericalSolver::get(x,y,z): z<0!" << std::endl; exit(EXIT_FAILURE);}
      if (z>_domainSizeZ+1){std::cout << "ERROR NumericalSolver::get(x,y,z): z>max. Value!" << std::endl; exit(EXIT_FAILURE);}
      #endif
      return x + (_domainSizeX+2)*(y+(_domainSizeY+2)*z);
    }

    /** returns index in 2D parallel buffer with buffer dimensions lengthx+2,lengthy+2. Performs checks in debug mode */
    int getParBuf(int x,int y,int lengthx,int lengthy) const {
      #if (COUPLING_MD_DEBUG==COUPLING_MD_YES)
      if (x<0 || x>lengthx+1){std::cout << "ERROR NumericalSolver::getParBuf(...): x out of range!" << std::endl; exit(EXIT_FAILURE);}
      if (y<0 || y>lengthy+1){std::cout << "ERROR NumericalSolver::getParBuf(...): y out of range!" << std::endl; exit(EXIT_FAILURE);}
      #endif
      return x+(lengthx+2)*y;
    }

    /** create vtk plot if required */
    void plot() const {
      // only plot output if this is the correct timestep
      if (_plotEveryTimestep==-1){ return;}
      if (_counter%_plotEveryTimestep!=0){return;}

      //const tarch::la::Vector<3,unsigned int> coords(getProcessCoordinates()); // offset of domain for MPI-parallel simulations
      int rank = 0; // rank in MPI-parallel simulations
      #if (COUPLING_MD_PARALLEL==COUPLING_MD_YES)
      MPI_Comm_rank(MPI_COMM_WORLD,&rank);
      #endif

      std::stringstream ss; ss << _filestem << "_" << rank << "_" << _counter << ".vtk";
      std::ofstream file(ss.str().c_str());
      if (!file.is_open()){std::cout << "ERROR NumericalSolver::plot(): Could not open file " << ss.str() << "!" << std::endl; exit(EXIT_FAILURE);}
      std::stringstream flag, density, velocity;

      file << "# vtk DataFile Version 2.0" << std::endl;
      file << "MaMiCo NumericalSolver" << std::endl;
      file << "ASCII" << std::endl << std::endl;
      file << "DATASET STRUCTURED_GRID" << std::endl;
      file << "DIMENSIONS " << _domainSizeX+3 << " " << _domainSizeY+3 << " " << _domainSizeZ+3 << std::endl; // everything +1 cause of change in index
      file << "POINTS " << (_domainSizeX+3)*(_domainSizeY+3)*(_domainSizeZ+3) << " float" << std::endl;

      flag << "CELL_DATA " << (_domainSizeX+2)*(_domainSizeY+2)*(_domainSizeZ+2) << std::endl;
      flag << "SCALARS flag float 1" << std::endl;
      flag << "LOOKUP_TABLE default" << std::endl;

      density  << std::setprecision(12);
      density  << "SCALARS density float 1 " << std::endl;
      density  << "LOOKUP_TABLE default" << std::endl;

      velocity << std::setprecision(12);
      velocity << "VECTORS velocity float" << std::endl;

      // loop over domain (incl. boundary) and write point coordinates

      for (int z = -1; z < _domainSizeZ+2; z++){ for (int y = -1; y < _domainSizeY+2; y++){ for (int x = -1; x < _domainSizeX+2; x++){
        file << ((int)(_coords[0]*_avgDomainSizeX)+x)*_dx << " " <<  ((int)(_coords[1]*_avgDomainSizeY)+y)*_dx  << " " << ((int)(_coords[2]*_avgDomainSizeZ)+z)*_dx << std::endl;
      }}}

      // loop over domain (incl. boundary)
      for (int z = 1; z < _domainSizeZ+1; z++){ for (int y = 1; y < _domainSizeY+1; y++){ for (int x = 1; x < _domainSizeX+1; x++){ //CHANGE: start index used to be one
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

    /** returns true, if this rank is not of relevance for the LB simulation */
    bool skipRank() const {
      int rank = 0;
      #if (COUPLING_MD_PARALLEL==COUPLING_MD_YES)
      MPI_Comm_rank(MPI_COMM_WORLD,&rank);
      #endif
      return ((unsigned int)rank>_processes[0]*_processes[1]*_processes[2]-1);
    }

    enum Flag{FLUID=0,NO_SLIP=1,MOVING_WALL=2,PERIODIC=3,MD_BOUNDARY=4,PARALLEL_BOUNDARY=5};
    enum NbFlag{LEFT=0,RIGHT=1,BACK=2,FRONT=3,BOTTOM=4,TOP=5};

    const double _channelheight; //
    const double _dx; // actual mesh size
    const double _dt; // actual time step size
    const double _kinVisc;
    tarch::la::Vector<3,unsigned int> _processes; // domain decomposition on MPI rank basis
    const int _plotEveryTimestep; // number of time steps between vtk plots
    const std::string _filestem; // file stem for vtk plot
    const int _domainSizeX { getDomainSize(_channelheight,_dx,_processes,0) }; // domain size in x-direction
    const int _domainSizeY { getDomainSize(_channelheight,_dx,_processes,1) }; // domain size in y-direction
    const int _domainSizeZ { getDomainSize(_channelheight,_dx,_processes,2) }; // domain size in z-direction
    const int _avgDomainSizeX { getAvgDomainSize(_channelheight,_dx,_processes,0) }; // avg. domain size in MPI-parallel simulation in x-direction
    const int _avgDomainSizeY { getAvgDomainSize(_channelheight,_dx,_processes,1) }; // "" in y-direction
    const int _avgDomainSizeZ { getAvgDomainSize(_channelheight,_dx,_processes,2) }; // "" in z-direction
    const tarch::la::Vector<3,unsigned int> _coords{getProcessCoordinates()};
    int _counter{0}; // time step counter
    double *_vel{NULL};  // velocity field
    double *_density{NULL}; // density
    Flag *_flag{NULL}; // flag field
    #if (COUPLING_MD_PARALLEL==COUPLING_MD_YES)
    double *_sendBufferX{NULL}; // buffer to send data from left/right to right/left neighbour
    double *_recvBufferX{NULL};
    double *_sendBufferY{NULL}; // buffer to receive data from from left/right neighbour
    double *_recvBufferY{NULL};
    double *_sendBufferZ{NULL};
    double *_recvBufferZ{NULL};
    #endif
    const int _xO{_domainSizeX+2}; // offset for y-direction (lexicographic grid ordering)
    const int _yO{(_domainSizeX+2)*(_domainSizeY+2)}; // offset for z-direction
    tarch::la::Vector<6,int> _parallelNeighbours{(-1)}; // neighbour ranks
    // for coupling with MD
    tarch::la::Vector<3,int> _offset{(-1)};
    tarch::la::Vector<3,int> _globalNumberMacroscopicCells{(-1)}; // OH
};

#endif// _MOLECULARDYNAMICS_COUPLING_SOLVERS_NUMERICALSOLVER_H_
