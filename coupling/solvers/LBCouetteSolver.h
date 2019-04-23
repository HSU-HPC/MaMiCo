// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_SOLVERS_LBCOUETTESOLVER_H_
#define _MOLECULARDYNAMICS_COUPLING_SOLVERS_LBCOUETTESOLVER_H_

#include "coupling/CouplingMDDefinitions.h"
#include "tarch/la/Vector.h"
#include <cmath>
#include <sstream>
#include <fstream>
#include <iomanip>
#if defined(_OPENMP)
#include <omp.h>
#endif
#if (COUPLING_MD_PARALLEL==COUPLING_MD_YES)
#include <mpi.h>
#endif
#include "coupling/solvers/CouetteSolver.h"
#include "coupling/datastructures/MacroscopicCell.h"
#include "coupling/IndexConversion.h"

namespace coupling {
  namespace solvers{
    class LBCouetteSolver;
  }
}


/** implements an analytic Couette flow solver.
 *  In our scenario, the lower wall is accelerated and the upper wall stands still.
 *  The lower wall is located at zero height.
 *  @author Philipp Neumann
 */
class coupling::solvers::LBCouetteSolver: public coupling::solvers::AbstractCouetteSolver<3> {
  private:
    enum Flag{FLUID=0,NO_SLIP=1,MOVING_WALL=2,PERIODIC=3,MD_BOUNDARY=4,PARALLEL_BOUNDARY=5};
    enum NbFlag{LEFT=0,RIGHT=1,BACK=2,FRONT=3,BOTTOM=4,TOP=5};

  public:
    LBCouetteSolver(
      const double channelheight,
      const tarch::la::Vector<3,double> wallVelocity,
      const double kinVisc,
      const double dx,
      const double dt,
      const int plotEveryTimestep,
      const std::string filestem,
      const tarch::la::Vector<3,unsigned int> processes,
      const unsigned int numThreads=1
    ):
    coupling::solvers::AbstractCouetteSolver<3>(),
    _omega(1.0/(3.0*(kinVisc*dt/(dx*dx))+0.5)),
    _wallVelocity((dt/dx)*wallVelocity),
    _domainSizeX(getDomainSize(channelheight,dx,processes,0)),
    _domainSizeY(getDomainSize(channelheight,dx,processes,1)),
    _domainSizeZ(getDomainSize(channelheight,dx,processes,2)),
    _avgDomainSizeX(getAvgDomainSize(channelheight,dx,processes,0)),
    _avgDomainSizeY(getAvgDomainSize(channelheight,dx,processes,1)),
    _avgDomainSizeZ(getAvgDomainSize(channelheight,dx,processes,2)),
    _time(0.0), _counter(0),
    _dx(dx),_dt(dt), _plotEveryTimestep(plotEveryTimestep), _filestem(filestem),
    _pdf1(NULL), _pdf2(NULL), _vel(NULL), _density(NULL), _flag(NULL),
    #if (COUPLING_MD_PARALLEL==COUPLING_MD_YES)
    _sendBufferX(NULL), _recvBufferX(NULL), _sendBufferY(NULL), _recvBufferY(NULL), _sendBufferZ(NULL), _recvBufferZ(NULL),
    #endif
    _C{{ 0,-1,-1}, {-1, 0,-1}, { 0, 0,-1}, { 1, 0,-1}, { 0, 1,-1},
       {-1,-1, 0}, { 0,-1, 0}, { 1,-1, 0}, {-1, 0, 0}, { 0, 0, 0}, { 1, 0, 0}, {-1, 1, 0}, { 0, 1, 0}, { 1, 1, 0},
       { 0,-1, 1}, {-1, 0, 1}, { 0, 0, 1}, { 1, 0, 1}, { 0, 1, 1}},
    _W{1.0/36.0, 1.0/36.0, 1.0/18.0, 1.0/36.0, 1.0/36.0,
       1.0/36.0, 1.0/18.0, 1.0/36.0, 1.0/18.0, 1.0/ 3.0, 1.0/18.0, 1.0/36.0, 1.0/18.0, 1.0/36.0,
       1.0/36.0, 1.0/36.0, 1.0/18.0, 1.0/36.0, 1.0/36.0},
    _xO(_domainSizeX+2), _yO((_domainSizeX+2)*(_domainSizeY+2)),
    _processes(processes),
    _parallelNeighbours(-1),
    _offset(-1),
    _globalNumberMacroscopicCells(-1)
    {
      _pdf1 = new double[19*(_domainSizeX+2)*(_domainSizeY+2)*(_domainSizeZ+2)];
      _pdf2 = new double[19*(_domainSizeX+2)*(_domainSizeY+2)*(_domainSizeZ+2)];
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

      #if defined(_OPENMP)
      omp_set_num_threads(numThreads);
      #endif
      #if (COUPLING_MD_DEBUG==COUPLING_MD_YES)
      std::cout << "Domain size=" << _domainSizeX << "," << _domainSizeY << "," << _domainSizeZ << std::endl;
      std::cout << "tau=" << 1.0/_omega << std::endl;
      std::cout << "wallVelocity=" << _wallVelocity << std::endl;
      #endif
      // check pointers
      if ( (_pdf1==NULL) || (_pdf2==NULL) || (_vel==NULL) || (_density==NULL) || (_flag==NULL) ){
        std::cout << "ERROR LBCouetteSolver: NULL ptr!" << std::endl; exit(EXIT_FAILURE);
      }
      #if (COUPLING_MD_PARALLEL==COUPLING_MD_YES)
      if ( (_sendBufferX==NULL) || (_recvBufferX==NULL) || (_sendBufferY==NULL) || (_recvBufferY==NULL) || (_sendBufferZ==NULL) || (_recvBufferZ==NULL) ){
        std::cout << "ERROR LBCouetteSolver: NULL ptr in send/recv!" << std::endl; exit(EXIT_FAILURE);
      }
      #endif

      // return if required
      if (skipRank()){return;}

      // init everything with lattice weights, zero velocity, unit density; flags are set to FLUID
      #if defined(_OPENMP)
      #pragma omp parallel for
      #endif
      for (int i = 0; i < (_domainSizeX+2)*(_domainSizeY+2)*(_domainSizeZ+2); i++){
        for (int q = 0; q < 19; q++){ _pdf1[get(i)*19+q] = _W[q]; _pdf2[get(i)*19+q] = _W[q]; }
        for (int d = 0; d <  3; d++){ _vel[get(i)*3+d]  = 0.0;}
        _density[get(i)] = 1.0;
        _flag[get(i)] = FLUID;
      }
      // determine parallel neighbours
      determineParallelNeighbours();

      // correct boundary flags based on physical description (Couette scenario)
      // bottom - moving wall
      #if defined(_OPENMP)
      #pragma omp parallel for
      #endif
      for (int i = 0; i < (_domainSizeX+2)*(_domainSizeY+2); i++){ _flag[get(i)] = MOVING_WALL; }
      // top - noslip
      #if defined(_OPENMP)
      #pragma omp parallel for
      #endif
      for (int i = (_domainSizeX+2)*(_domainSizeY+2)*(_domainSizeZ+1); i < (_domainSizeX+2)*(_domainSizeY+2)*(_domainSizeZ+2); i++){ _flag[get(i)] = NO_SLIP;}
      // front - periodic
      #if defined(_OPENMP)
      #pragma omp parallel for
      #endif
      for (int z = 1; z < _domainSizeZ+1; z++){for (int x = 0; x < _domainSizeX+2; x++){ _flag[get(x,0,z)] = PERIODIC; } }
      // back - periodic
      #if defined(_OPENMP)
      #pragma omp parallel for
      #endif
      for (int z = 1; z < _domainSizeZ+1; z++){for (int x = 0; x < _domainSizeX+2; x++){ _flag[get(x,_domainSizeY+1,z)] = PERIODIC; } }
      // left - periodic
      #if defined(_OPENMP)
      #pragma omp parallel for
      #endif
      for (int z = 1; z < _domainSizeZ+1; z++){for (int y = 1; y < _domainSizeY+1; y++){ _flag[get(0,y,z)] = PERIODIC; }}
      #if defined(_OPENMP)
      #pragma omp parallel for
      #endif
      for (int z = 1; z < _domainSizeZ+1; z++){for (int y = 1; y < _domainSizeY+1; y++){ _flag[get(_domainSizeX+1,y,z)] = PERIODIC; }}

      // correct boundary flags in case of MPI-parallel simulations (Couette scenario)
      setParallelBoundaryFlags();

      #if (COUPLING_MD_DEBUG==COUPLING_MD_YES)
      for (int z = 0; z < _domainSizeZ+2; z++){ for (int y = 0; y < _domainSizeY+2; y++){ for (int x = 0; x < _domainSizeX+2; x++){
        std::cout << x << "," << y << "," << z << "FLAG=" << _flag[get(x,y,z)] << std::endl;
      }}}
      #endif
    }
    virtual ~LBCouetteSolver(){
      if (_pdf1!=NULL){delete [] _pdf1; _pdf1=NULL;}
      if (_pdf2!=NULL){delete [] _pdf2; _pdf2=NULL;}
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

    /** advances one time step dt in time and triggers vtk plot if required */
    virtual void advance(double dt){
      if (skipRank()){return;}
      const int timesteps=floor( dt/_dt+0.5 );
      if ( fabs(timesteps*_dt-dt)/_dt > 1.0e-8 ){std::cout << "ERROR LBCouetteSolver::advance(): time steps and dt do not match!" << std::endl; exit(EXIT_FAILURE);}
      for (int i = 0; i < timesteps; i++){
        plot();
        collidestream();
        communicate(); // exchange between neighbouring MPI subdomains
        _time += _dt;
        _counter++;
      }
    }

    virtual void setWallVelocity(const tarch::la::Vector<3,double> wallVelocity){
      _wallVelocity = (_dt/_dx)*wallVelocity;
    }

    /** flags the domain boundary cells. mdDomainOffset and mdDomainSize correspond to lower/left/front corner of the MD domain and the size of the domain. overlapStrip is the number of cells
     *  that the LB and MD domain shall overlap, i.e. the number of "MD cells" that lie within the LB domain but which should not be flagged (and thus be handled by both solvers).
     */
    void setMDBoundary(tarch::la::Vector<3,double> mdDomainOffset,tarch::la::Vector<3,double> mdDomainSize,unsigned int overlapStrip){
      if (skipRank()){return;}
      const tarch::la::Vector<3,unsigned int> coords(getProcessCoordinates());

      for (int d = 0; d < 3; d++){
        _offset[d] = (floor(mdDomainOffset[d]/_dx + 0.5));
        if (fabs(_offset[d]*_dx - mdDomainOffset[d])/_dx>1.0e-8 ){std::cout << "ERROR LBCouetteSolver::setMDBoundary(): offset does not match!" << std::endl; exit(EXIT_FAILURE);}
        _globalNumberMacroscopicCells[d] = (floor(mdDomainSize[d]/_dx + 0.5));
        if (fabs(_globalNumberMacroscopicCells[d]*_dx - mdDomainSize[d])/_dx>1.0e-8 ){std::cout << "ERROR LBCouetteSolver::setMDBoundary(): globalNumber does not match!" << std::endl; exit(EXIT_FAILURE);}
      }

      // flag local domain
      for (int z = 0; z < _domainSizeZ+2; z++){
        for (int y = 0; y < _domainSizeY+2; y++){
          for (int x = 0; x < _domainSizeX+2; x++){
            // determine global cell coordinates of the process-local (sub)domain
            tarch::la::Vector<3,int> globalCoords(x-1+coords[0]*_avgDomainSizeX,y-1+coords[1]*_avgDomainSizeY,z-1+coords[2]*_avgDomainSizeZ);
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

    void setMDBoundaryValues(
      std::vector<coupling::datastructures::MacroscopicCell<3>* >& recvBuffer,const unsigned int * const recvIndices,
      const coupling::IndexConversion<3>& indexConversion
    ){
      if (skipRank()){return ;}
      // determine process coordinates
      const tarch::la::Vector<3,unsigned int> coords(getProcessCoordinates());

      // loop over all received cells
      const unsigned int size = (unsigned int) recvBuffer.size();
      for (unsigned int i = 0; i < size; i++){
        // determine cell index of this cell in LB domain
        tarch::la::Vector<3,unsigned int> globalCellCoords = indexConversion.getGlobalVectorCellIndex(recvIndices[i]);
        globalCellCoords[0] = (globalCellCoords[0]+_offset[0]) - coords[0]*_avgDomainSizeX;
        globalCellCoords[1] = (globalCellCoords[1]+_offset[1]) - coords[1]*_avgDomainSizeY;
        globalCellCoords[2] = (globalCellCoords[2]+_offset[2]) - coords[2]*_avgDomainSizeZ;
//std::cout << "Process coords: " << coords << ":  GlobalCellCoords for index " << indexConversion.getGlobalVectorCellIndex(recvIndices[i]) << ": " << globalCellCoords << std::endl;
        const int index = get(globalCellCoords[0],globalCellCoords[1],globalCellCoords[2]);
        #if (COUPLING_MD_DEBUG==COUPLING_MD_YES)
        if (_flag[index]!=MD_BOUNDARY){std::cout << "ERROR LBCouetteSolver::setMDBoundaryValues(): Cell " << index << " is no MD boundary cell!" << std::endl; exit(EXIT_FAILURE);}
        #endif
        // set velocity value and pdfs in MD boundary cell (before streaming); the boundary velocities are interpolated between the neighbouring and this cell. This interpolation is valid for
        // FLUID-MD_BOUNDARY neighbouring relations only.
        // determine local velocity received from MaMiCo and convert it to LB units; store the velocity in _vel
        tarch::la::Vector<3,double> localVel( (1.0/recvBuffer[i]->getMacroscopicMass())*(_dt/_dx)*recvBuffer[i]->getMacroscopicMomentum());
        for (unsigned int d = 0; d<3; d++){ _vel[3*index+d] = localVel[d]; }
        // loop over all pdfs and set them according to interpolated moving-wall conditions
        for (unsigned int q = 0; q < 19; q++){
          // index of neighbour cell; only if cell is located inside local domain
          if (   ((int)globalCellCoords[0]+_C[q][0] > 0) && ((int)globalCellCoords[0]+_C[q][0]<_domainSizeX+1)
              && ((int)globalCellCoords[1]+_C[q][1] > 0) && ((int)globalCellCoords[1]+_C[q][1]<_domainSizeY+1)
              && ((int)globalCellCoords[2]+_C[q][2] > 0) && ((int)globalCellCoords[2]+_C[q][2]<_domainSizeZ+1)){
            const int nbIndex = get((_C[q][0]+globalCellCoords[0]),(_C[q][1]+globalCellCoords[1]),(_C[q][2]+globalCellCoords[2]));
            const tarch::la::Vector<3,double> interpolVel(0.5*(_vel[3*index]+_vel[3*nbIndex]),0.5*(_vel[3*index+1]+_vel[3*nbIndex+1]),0.5*(_vel[3*index+2]+_vel[3*nbIndex+2]));
            _pdf1[19*index+q] = _pdf1[19*nbIndex+18-q] - 6.0*_W[q]*_density[nbIndex]*(_C[18-q][0]*interpolVel[0]+_C[18-q][1]*interpolVel[1]+_C[18-q][2]*interpolVel[2]);
          }
        }
      }
    }

    /**
     * returns velocity at a certain position.
     */
    virtual tarch::la::Vector<3,double> getVelocity(tarch::la::Vector<3,double> pos) const {
      tarch::la::Vector<3,unsigned int> coords(getProcessCoordinates());
      const tarch::la::Vector<3,double> domainOffset(coords[0]*_dx*_avgDomainSizeX,coords[1]*_dx*_avgDomainSizeY,coords[2]*_dx*_avgDomainSizeZ);

      // check pos-data for process locality (todo: put this in debug mode in future releases)
      if (   (pos[0]<domainOffset[0]) || (pos[0]>domainOffset[0]+_domainSizeX*_dx)
          || (pos[1]<domainOffset[1]) || (pos[1]>domainOffset[1]+_domainSizeY*_dx)
          || (pos[2]<domainOffset[2]) || (pos[2]>domainOffset[2]+_domainSizeZ*_dx) ){
        std::cout << "ERROR LBCouetteSolver::getVelocity(): Position " << pos << " out of range!" << std::endl; exit(EXIT_FAILURE);
      }
      // compute index for respective cell (_dx+... for ghost cells); use coords to store local cell coordinates
      for (unsigned int d = 0; d < 3; d++){ coords[d] = (unsigned int) ((_dx+pos[d]-domainOffset[d])/_dx);}
      const int index = get(coords[0],coords[1],coords[2]);
      tarch::la::Vector<3,double> vel(0.0);
      // extract and scale velocity to "real"=MD units
      for (int d = 0; d < 3; d++){ vel[d] = _dx/_dt*_vel[3*index+d]; }
      #if (COUPLING_MD_DEBUG==COUPLING_MD_YES)
      std::cout << "Position " << pos << " corresponds to cell: " << coords << "; vel=" << vel << std::endl;
      #endif
      return vel;
    }

    /**
     * returns density IN LB UNITS at a certain position.
     */
     double getDensity(tarch::la::Vector<3,double> pos) const {
      tarch::la::Vector<3,unsigned int> coords(getProcessCoordinates());
      const tarch::la::Vector<3,double> domainOffset(coords[0]*_dx*_avgDomainSizeX,coords[1]*_dx*_avgDomainSizeY,coords[2]*_dx*_avgDomainSizeZ);

      // check pos-data for process locality (todo: put this in debug mode in future releases)
      if (   (pos[0]<domainOffset[0]) || (pos[0]>domainOffset[0]+_domainSizeX*_dx)
          || (pos[1]<domainOffset[1]) || (pos[1]>domainOffset[1]+_domainSizeY*_dx)
          || (pos[2]<domainOffset[2]) || (pos[2]>domainOffset[2]+_domainSizeZ*_dx) ){
        std::cout << "ERROR LBCouetteSolver::getDensity(): Position " << pos << " out of range!" << std::endl; exit(EXIT_FAILURE);
      }
      // compute index for respective cell (_dx+... for ghost cells); use coords to store local cell coordinates
      for (unsigned int d = 0; d < 3; d++){ coords[d] = (unsigned int) ((_dx+pos[d]-domainOffset[d])/_dx);}
      const int index = get(coords[0],coords[1],coords[2]);
      return _density[index];
    }

    /** getters required by LB Couette Solver Interface */
    tarch::la::Vector<3,unsigned int> getNumberProcesses() const { return _processes;}
    tarch::la::Vector<3,unsigned int> getAvgNumberLBCells() const { tarch::la::Vector<3,unsigned int> avgCells(_avgDomainSizeX,_avgDomainSizeY,_avgDomainSizeZ); return avgCells;}

  private:
    /** returns i and performs checks in debug mode */
    int get(int i) const {
      #if (COUPLING_MD_DEBUG==COUPLING_MD_YES)
      if (i<0){std::cout << "ERROR LBCouetteSolver::get(i): i<0!" << std::endl; exit(EXIT_FAILURE);}
      if (i>(_domainSizeX+2)*(_domainSizeY+2)*(_domainSizeZ+2)){std::cout << "ERROR LBCouetteSolver::get(i): i>max. Value!" << std::endl; exit(EXIT_FAILURE);}
      #endif
      return i;
    }

    /** returns linearized index and performs checks in debug mode */
    int get(int x,int y,int z) const{
      #if (COUPLING_MD_DEBUG==COUPLING_MD_YES)
      if (x<0){std::cout << "ERROR LBCouetteSolver::get(x,y,z): x<0!" << std::endl; exit(EXIT_FAILURE);}
      if (x>_domainSizeX+1){std::cout << "ERROR LBCouetteSolver::get(x,y,z): x>max. Value!" << std::endl; exit(EXIT_FAILURE);}
      if (y<0){std::cout << "ERROR LBCouetteSolver::get(x,y,z): y<0!" << std::endl; exit(EXIT_FAILURE);}
      if (y>_domainSizeY+1){std::cout << "ERROR LBCouetteSolver::get(x,y,z): y>max. Value!" << std::endl; exit(EXIT_FAILURE);}
      if (z<0){std::cout << "ERROR LBCouetteSolver::get(x,y,z): z<0!" << std::endl; exit(EXIT_FAILURE);}
      if (z>_domainSizeZ+1){std::cout << "ERROR LBCouetteSolver::get(x,y,z): z>max. Value!" << std::endl; exit(EXIT_FAILURE);}
      #endif
      return x + (_domainSizeX+2)*(y+(_domainSizeY+2)*z);
    }

    /** returns index in 2D parallel buffer with buffer dimensions lengthx+2,lengthy+2. Performs checks in debug mode */
    int getParBuf(int x,int y,int lengthx,int lengthy) const {
      #if (COUPLING_MD_DEBUG==COUPLING_MD_YES)
      if (x<0 || x>lengthx+1){std::cout << "ERROR LBCouetteSolver::getParBuf(...): x out of range!" << std::endl; exit(EXIT_FAILURE);}
      if (y<0 || y>lengthy+1){std::cout << "ERROR LBCouetteSolver::getParBuf(...): y out of range!" << std::endl; exit(EXIT_FAILURE);}
      #endif
      return x+(lengthx+2)*y;
    }

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

    /** collide-stream algorithm */
    void collidestream(){
      #if defined(_OPENMP)
      #pragma omp parallel for
      #endif
      for (int z = 1; z < _domainSizeZ+1; z++){
        for (int y = 1; y < _domainSizeY+1; y++){
          for (int x = 1; x < _domainSizeX+1; x++){
            const int index = get(x,y,z);
            if (_flag[index]==FLUID){
              stream(index);
              collide(index,x,y,z);
            }
          }
        }
      }
      // swap fields
      double *swap=_pdf1;
      _pdf1 = _pdf2;
      _pdf2 = swap;
    }
    
    /** stream from pdf1 to pdf2 */
    void stream(int index){
      const int pI = 19*index;
      for (int q = 0; q < 9; q++){
         const int nb= 19*(_C[q][0]+_C[q][1]*_xO+_C[q][2]*_yO);
         _pdf2[pI+q]    = _pdf1[pI+q-nb];
         _pdf2[pI+18-q] = _pdf1[pI+18-q+nb];
      }
      _pdf2[pI+9] = _pdf1[pI+9];
    }

    /** collide within pdf2 */
    void collide(int index,int x, int y, int z){
      // index of start of cell-local pdfs in AoS
      const int pI = 19*index;

      // compute and store density, velocity
      double *vel = &_vel[3*index];
      computeDensityAndVelocity(vel,_density[index],&_pdf2[pI]);

      // collide (BGK); always handle pdfs no. q and inv(q)=18-q in one step
      const double u2 = 1.0 - 1.5*(vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2]);
      // pdf 0,18
      double cu = -vel[1]-vel[2];
      int nb = -_xO-_yO;
      double feq = _W[0]*_density[index]*(u2 + 3.0*cu + 4.5*cu*cu);
      _pdf2[pI] -= _omega*(_pdf2[pI] - feq);
      boundary(_pdf2,pI,x,y,z,0,_flag[index+nb],pI+19*nb);
      feq -= 6.0*_W[0]*_density[index]*cu;
      _pdf2[pI+18] -= _omega*(_pdf2[pI+18] - feq);
      boundary(_pdf2,pI,x,y,z,18,_flag[index-nb],pI-19*nb);
      // pdf 1,17
      cu = -vel[0]-vel[2];
      nb = -1-_yO;
      feq = _W[1]*_density[index]*(u2 + 3.0*cu + 4.5*cu*cu);
      _pdf2[pI+1] -= _omega*(_pdf2[pI+1] - feq);
      boundary(_pdf2,pI,x,y,z,1,_flag[index+nb],pI+19*nb);
      feq -= 6.0*_W[1]*_density[index]*cu;
      _pdf2[pI+17] -= _omega*(_pdf2[pI+17] - feq);
      boundary(_pdf2,pI,x,y,z,17,_flag[index-nb],pI-19*nb);
      // pdf 2,16
      cu = -vel[2];
      nb = -_yO;
      feq = _W[2]*_density[index]*(u2 + 3.0*cu + 4.5*cu*cu);
      _pdf2[pI+2] -= _omega*(_pdf2[pI+2] - feq);
      boundary(_pdf2,pI,x,y,z,2,_flag[index+nb],pI+19*nb);
      feq -= 6.0*_W[2]*_density[index]*cu;
      _pdf2[pI+16] -= _omega*(_pdf2[pI+16] - feq);
      boundary(_pdf2,pI,x,y,z,16,_flag[index-nb],pI-19*nb);
      // pdf 3,15
      cu = vel[0]-vel[2];
      nb = 1-_yO;
      feq = _W[3]*_density[index]*(u2 + 3.0*cu + 4.5*cu*cu);
      _pdf2[pI+3] -= _omega*(_pdf2[pI+3] - feq);
      boundary(_pdf2,pI,x,y,z,3,_flag[index+nb],pI+19*nb);
      feq -= 6.0*_W[3]*_density[index]*cu;
      _pdf2[pI+15] -= _omega*(_pdf2[pI+15] - feq);
      boundary(_pdf2,pI,x,y,z,15,_flag[index-nb],pI-19*nb);
      // pdf 4,14
      cu = vel[1]-vel[2];
      nb = _xO-_yO;
      feq = _W[4]*_density[index]*(u2 + 3.0*cu + 4.5*cu*cu);
      _pdf2[pI+4] -= _omega*(_pdf2[pI+4] - feq);
      boundary(_pdf2,pI,x,y,z,4,_flag[index+nb],pI+19*nb);
      feq -= 6.0*_W[4]*_density[index]*cu;
      _pdf2[pI+14] -= _omega*(_pdf2[pI+14] - feq);
      boundary(_pdf2,pI,x,y,z,14,_flag[index-nb],pI-19*nb);
      // pdf 5,13
      cu = -vel[0]-vel[1];
      nb = -1-_xO;
      feq = _W[5]*_density[index]*(u2 + 3.0*cu + 4.5*cu*cu);
      _pdf2[pI+5] -= _omega*(_pdf2[pI+5] - feq);
      boundary(_pdf2,pI,x,y,z,5,_flag[index+nb],pI+19*nb);
      feq -= 6.0*_W[5]*_density[index]*cu;
      _pdf2[pI+13] -= _omega*(_pdf2[pI+13] - feq);
      boundary(_pdf2,pI,x,y,z,13,_flag[index-nb],pI-19*nb);
      // pdf 6,12
      cu = -vel[1];
      nb = -_xO;
      feq = _W[6]*_density[index]*(u2 + 3.0*cu + 4.5*cu*cu);
      _pdf2[pI+6] -= _omega*(_pdf2[pI+6] - feq);
      boundary(_pdf2,pI,x,y,z,6,_flag[index+nb],pI+19*nb);
      feq -= 6.0*_W[6]*_density[index]*cu;
      _pdf2[pI+12] -= _omega*(_pdf2[pI+12] - feq);
      boundary(_pdf2,pI,x,y,z,12,_flag[index-nb],pI-19*nb);
      // pdf 7,11
      cu = vel[0]-vel[1];
      nb = 1-_xO;
      feq = _W[7]*_density[index]*(u2 + 3.0*cu + 4.5*cu*cu);
      _pdf2[pI+7] -= _omega*(_pdf2[pI+7] - feq);
      boundary(_pdf2,pI,x,y,z,7,_flag[index+nb],pI+19*nb);
      feq -= 6.0*_W[7]*_density[index]*cu;
      _pdf2[pI+11] -= _omega*(_pdf2[pI+11] - feq);
      boundary(_pdf2,pI,x,y,z,11,_flag[index-nb],pI-19*nb);
      // pdf 8,10
      cu = -vel[0];
      nb = -1;
      feq = _W[8]*_density[index]*(u2 + 3.0*cu + 4.5*cu*cu);
      _pdf2[pI+8] -= _omega*(_pdf2[pI+8] - feq);
      boundary(_pdf2,pI,x,y,z,8,_flag[index+nb],pI+19*nb);
      feq -= 6.0*_W[8]*_density[index]*cu;
      _pdf2[pI+10] -= _omega*(_pdf2[pI+10] - feq);
      boundary(_pdf2,pI,x,y,z,10,_flag[index-nb],pI-19*nb);
      // pdf 9
      _pdf2[pI+9] -= _omega*(_pdf2[pI+9] - _W[9]*_density[index]*u2);
    }

    // boundary treatment; pdf - distributions (typically _pdf2), index - start index for current cell in pdf-array, x/y/z - cell coordinates of current cell,
    // q - distribution function no., flag - boundary flag of neighbouring cell, nbIndex - index of neighbouring cell
    void boundary(double * const pdf, int index, int x, int y, int z, int q, const Flag &flag,int nbIndex){
      if (flag!=FLUID){
      if (flag==NO_SLIP){
        // half-way bounce back
        pdf[nbIndex+18-q] = pdf[index+q];
      } else if (flag==MOVING_WALL){
        // half-way bounce back + moving wall acceleration (only x-direction for wall supported at the moment)
        pdf[nbIndex+18-q] = pdf[index+q] - 6.0*_W[q]*_density[index/19]*(_C[q][0]*_wallVelocity[0]+_C[q][1]*_wallVelocity[1]+_C[q][2]*_wallVelocity[2]);
      } else if (flag==PERIODIC){
        // periodic treatment
        int target[3] = {x,y,z};
        if (target[0]+_C[q][0]==0){target[0] = _domainSizeX+1;} else if (target[0]+_C[q][0]==_domainSizeX+1){target[0] = 0;}
        if (target[1]+_C[q][1]==0){target[1] = _domainSizeY+1;} else if (target[1]+_C[q][1]==_domainSizeY+1){target[1] = 0;}
        if (target[2]+_C[q][2]==0){target[2] = _domainSizeZ+1;} else if (target[2]+_C[q][2]==_domainSizeZ+1){target[2] = 0;}
        const int periodicNb = target[0] + (_domainSizeX+2)*(target[1] + (_domainSizeY+2)*target[2]);
        pdf[19*periodicNb+q] = pdf[index+q];
      }
      }
    }

    /** create vtk plot if required */
    void plot()  {
      // only plot output if this is the correct timestep
      if (_plotEveryTimestep==-1){ return;}
      if (_counter%_plotEveryTimestep!=0){return;}

      const tarch::la::Vector<3,unsigned int> coords(getProcessCoordinates()); // offset of domain for MPI-parallel simulations
      int rank = 0; // rank in MPI-parallel simulations
      #if (COUPLING_MD_PARALLEL==COUPLING_MD_YES)
      MPI_Comm_rank(MPI_COMM_WORLD,&rank);
      #endif

      std::stringstream ss; ss << _filestem << "_" << rank << "_" << _counter << ".vtk";
      std::ofstream file(ss.str().c_str());
      if (!file.is_open()){std::cout << "ERROR LBCouetteSolver::plot(): Could not open file " << ss.str() << "!" << std::endl; exit(EXIT_FAILURE);}
      std::stringstream flag, density, velocity;

      file << "# vtk DataFile Version 2.0" << std::endl;
      file << "MaMiCo LBCouetteSolver" << std::endl;
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
        file << (coords[0]*_avgDomainSizeX+x)*_dx << " " <<  (coords[1]*_avgDomainSizeY+y)*_dx  << " " << (coords[2]*_avgDomainSizeZ+z)*_dx << std::endl;
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

    /** compute density and velocity on pdf */
    void computeDensityAndVelocity(double * const vel, double &density, const double * const pdf){
      vel[0] = -(pdf[1]+pdf[5]+pdf[8]+pdf[11]+pdf[15]);
      density=   pdf[3]+pdf[7]+pdf[10]+pdf[13]+pdf[17];
      vel[1] = (pdf[4]+pdf[11]+pdf[12]+pdf[13]+pdf[18]) - (pdf[0]+pdf[5]+pdf[6]+pdf[7]+pdf[14]);
      vel[0] = density + vel[0];
      density= density + pdf[0]+pdf[1]+pdf[2] + pdf[4]+pdf[5]+pdf[6] + pdf[8]+pdf[9] + pdf[11]+pdf[12] + pdf[14]+pdf[15]+pdf[16] + pdf[18];
      vel[2] = (pdf[14]+pdf[15]+pdf[16]+pdf[17]+pdf[18]) - (pdf[0]+pdf[1]+pdf[2]+pdf[3]+pdf[4]);
      vel[0] = vel[0]/density;
      vel[1] = vel[1]/density;
      vel[2] = vel[2]/density;
    }

    void determineParallelNeighbours(){
      #if (COUPLING_MD_PARALLEL==COUPLING_MD_YES)
      int rank;
      int size;
      MPI_Comm_rank(MPI_COMM_WORLD,&rank);
      MPI_Comm_size(MPI_COMM_WORLD,&size);
      // check if enough ranks are available
      if (_processes[0]*_processes[1]*_processes[2]>(unsigned int)size){
        std::cout << "ERROR LBCouetteSolver::determineParallelNeighbours(): Not enough ranks available!" << std::endl; exit(EXIT_FAILURE);
      }
      // determine rank coordinates
      const tarch::la::Vector<3,unsigned int> coords(getProcessCoordinates());

      // neighbour dependencies based on Couette problem
      // left,right: periodic
      _parallelNeighbours[LEFT] = ((coords[0]+_processes[0]-1) % _processes[0]) + _processes[0]*(coords[1]+_processes[1]*coords[2]);
      _parallelNeighbours[RIGHT]= ((coords[0]+_processes[0]+1) % _processes[0]) + _processes[0]*(coords[1]+_processes[1]*coords[2]);
      // back,front: periodic
      _parallelNeighbours[FRONT]= coords[0] + _processes[0]*( ((coords[1]+_processes[1]-1) % _processes[1]) + _processes[1]*coords[2]);
      _parallelNeighbours[BACK] = coords[0] + _processes[0]*( ((coords[1]+_processes[1]+1) % _processes[1]) + _processes[1]*coords[2]);
      // top: either neighbour or MPI_PROC_NULL
      if (coords[2]==_processes[2]-1){ _parallelNeighbours[TOP] = MPI_PROC_NULL; } else { _parallelNeighbours[TOP] = coords[0]+_processes[0]*(coords[1]+_processes[1]*(coords[2]+1)); }
      // bottom either neighbour or MPI_PROC_NULL
      if (coords[2]==0)              { _parallelNeighbours[BOTTOM]=MPI_PROC_NULL;} else { _parallelNeighbours[BOTTOM]=coords[0]+_processes[0]*(coords[1]+_processes[1]*(coords[2]-1));}
      std::cout << "Parallel neighbours for rank " << rank << ": " << _parallelNeighbours << std::endl;
      #endif
    }

    /** sets parallel boundary flags according to Couette scenario */
    void setParallelBoundaryFlags(){
      #if (COUPLING_MD_PARALLEL==COUPLING_MD_YES)
      const tarch::la::Vector<3,unsigned int> coords(getProcessCoordinates());

      // bottom - moving wall
      if (coords[2]!=0){
        #if defined(_OPENMP)
        #pragma omp parallel for
        #endif
        for (int i = 0; i < (_domainSizeX+2)*(_domainSizeY+2); i++){ _flag[get(i)] = PARALLEL_BOUNDARY; }
      }
      // top - noslip
      if (coords[2]!=_processes[2]-1){
        #if defined(_OPENMP)
        #pragma omp parallel for
        #endif
        for (int i = (_domainSizeX+2)*(_domainSizeY+2)*(_domainSizeZ+1); i < (_domainSizeX+2)*(_domainSizeY+2)*(_domainSizeZ+2); i++){ _flag[get(i)] = PARALLEL_BOUNDARY;}
      }
      // for all other boundaries front,back,left,right, we use parallel boundaries. If we send from one processes to itself, this is still the same as periodic conditions
      // front - periodic
      #if defined(_OPENMP)
      #pragma omp parallel for
      #endif
      for (int z = 1; z < _domainSizeZ+1; z++){for (int x = 0; x < _domainSizeX+2; x++){ _flag[get(x,0,z)] = PARALLEL_BOUNDARY; } }
      // back - periodic
      #if defined(_OPENMP)
      #pragma omp parallel for
      #endif
      for (int z = 1; z < _domainSizeZ+1; z++){for (int x = 0; x < _domainSizeX+2; x++){ _flag[get(x,_domainSizeY+1,z)] = PARALLEL_BOUNDARY; } }
      // left - periodic
      #if defined(_OPENMP)
      #pragma omp parallel for
      #endif
      for (int z = 1; z < _domainSizeZ+1; z++){for (int y = 1; y < _domainSizeY+1; y++){ _flag[get(0,y,z)] = PARALLEL_BOUNDARY; }}
      #if defined(_OPENMP)
      #pragma omp parallel for
      #endif
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

    /** takes care of communication across one face in one direction.
     *  pdf - pdf field
     *  sendBuffer, recvBuffer - send and recv buffer
     *  nbFlagTo - direction into which message is sent
     *  nbFlagFrom - direction from which message is received
     *  startSend,endSend - 3-d coordinates that define a range to be sent to neighbouring process
     *  startRecv,endRecv - 3-d coordinates that define a range to be received from neighbouring process
     */
    void communicatePart(double *pdf, double *sendBuffer, double *recvBuffer, NbFlag nbFlagTo, NbFlag nbFlagFrom,
      tarch::la::Vector<3,int> startSend, tarch::la::Vector<3,int> endSend,
      tarch::la::Vector<3,int> startRecv, tarch::la::Vector<3,int> endRecv
    ){
      #if (COUPLING_MD_PARALLEL==COUPLING_MD_YES)
      // directions that point to LEFT/RIGHT,... -> same ordering as enums!
      const int directions[6][5] = { { 1, 5, 8,11,15}, { 3, 7,10,13,17},
                                     { 4,11,12,13,18}, { 0, 5, 6, 7,14},
                                     { 0, 1, 2, 3, 4}, {14,15,16,17,18} };
      MPI_Request requests[2];
      MPI_Status  status[2];
      tarch::la::Vector<2,int> plane;
      tarch::la::Vector<2,int> domainSize;
      // find out plane coordinates
      if        (nbFlagTo==LEFT || nbFlagTo==RIGHT){plane[0]=1; plane[1]=2; domainSize[0] = _domainSizeY; domainSize[1] = _domainSizeZ;
      } else if (nbFlagTo==FRONT|| nbFlagTo==BACK) {plane[0]=0; plane[1]=2; domainSize[0] = _domainSizeX; domainSize[1] = _domainSizeZ;
      } else if (nbFlagTo==TOP || nbFlagTo==BOTTOM){plane[0]=0; plane[1]=1; domainSize[0] = _domainSizeX; domainSize[1] = _domainSizeY; }
      else { std::cout << "ERROR LBCouetteSolver::communicatePart: d >2 or d < 0!" << std::endl; exit(EXIT_FAILURE);}

      // extract data and write to send buffer
      tarch::la::Vector<3,int> coords(0);
      for (coords[2] = startSend[2]; coords[2] < endSend[2]; coords[2]++){
      for (coords[1] = startSend[1]; coords[1] < endSend[1]; coords[1]++){
      for (coords[0] = startSend[0]; coords[0] < endSend[0]; coords[0]++){
        for (int q = 0; q < 5; q++){
          sendBuffer[q+5*getParBuf(coords[plane[0]],coords[plane[1]],domainSize[0],domainSize[1])] = pdf[directions[nbFlagTo][q]+19*get(coords[0],coords[1],coords[2])];
        }
      }}}

      // send and receive data
      MPI_Irecv(recvBuffer,(domainSize[0]+2)*(domainSize[1]+2)*5,MPI_DOUBLE,_parallelNeighbours[nbFlagFrom],1000,MPI_COMM_WORLD,&requests[0]);
      MPI_Isend(sendBuffer,(domainSize[0]+2)*(domainSize[1]+2)*5,MPI_DOUBLE,_parallelNeighbours[nbFlagTo],  1000,MPI_COMM_WORLD,&requests[1]);
      MPI_Waitall(2,requests,status);

      // write data back to pdf field
      if (_parallelNeighbours[nbFlagFrom]!=MPI_PROC_NULL){
        for (coords[2] = startRecv[2]; coords[2] < endRecv[2]; coords[2]++){
        for (coords[1] = startRecv[1]; coords[1] < endRecv[1]; coords[1]++){
        for (coords[0] = startRecv[0]; coords[0] < endRecv[0]; coords[0]++){
          for (int q = 0; q < 5; q++){
            if (_flag[get(coords[0],coords[1],coords[2])] == PARALLEL_BOUNDARY){
              pdf[directions[nbFlagTo][q]+19*get(coords[0],coords[1],coords[2])] = recvBuffer[q+5*getParBuf(coords[plane[0]],coords[plane[1]],domainSize[0],domainSize[1])];
            }
          }
        }}}
      }
      #endif
    }

    void communicate(){
      #if (COUPLING_MD_PARALLEL==COUPLING_MD_YES)
      // send from right to left
      communicatePart(_pdf1,_sendBufferX,_recvBufferX,LEFT,RIGHT,
        tarch::la::Vector<3,int>(1,1,1),             tarch::la::Vector<3,int>(2,_domainSizeY+1,_domainSizeZ+1),
        tarch::la::Vector<3,int>(_domainSizeX+1,1,1),tarch::la::Vector<3,int>(_domainSizeX+2,_domainSizeY+1,_domainSizeZ+1));
      // send from left to right
      communicatePart(_pdf1,_sendBufferX,_recvBufferX,RIGHT,LEFT,
        tarch::la::Vector<3,int>(_domainSizeX,1,1),tarch::la::Vector<3,int>(_domainSizeX+1,_domainSizeY+1,_domainSizeZ+1),
        tarch::la::Vector<3,int>(0,1,1),           tarch::la::Vector<3,int>(1,_domainSizeY+1,_domainSizeZ+1));

      // send from back to front
      communicatePart(_pdf1,_sendBufferY,_recvBufferY,FRONT,BACK,
        tarch::la::Vector<3,int>(0,1,1),              tarch::la::Vector<3,int>(_domainSizeX+2,2,_domainSizeZ+1),
        tarch::la::Vector<3,int>(0,_domainSizeY+1,1), tarch::la::Vector<3,int>(_domainSizeX+2,_domainSizeY+2,_domainSizeZ+1));
      // send from front to back
      communicatePart(_pdf1,_sendBufferY,_recvBufferY,BACK,FRONT,
        tarch::la::Vector<3,int>(0,_domainSizeY,1), tarch::la::Vector<3,int>(_domainSizeX+2,_domainSizeY+1,_domainSizeZ+1),
        tarch::la::Vector<3,int>(0,0,1),            tarch::la::Vector<3,int>(_domainSizeX+2,1,_domainSizeZ+1));

      // send from top to bottom
      communicatePart(_pdf1,_sendBufferZ,_recvBufferZ,BOTTOM,TOP,
        tarch::la::Vector<3,int>(0,0,1),              tarch::la::Vector<3,int>(_domainSizeX+2,_domainSizeY+2,2),
        tarch::la::Vector<3,int>(0,0,_domainSizeZ+1), tarch::la::Vector<3,int>(_domainSizeX+2,_domainSizeY+2,_domainSizeZ+2));
      // send from bottom to top
      communicatePart(_pdf1,_sendBufferZ,_recvBufferZ,TOP,BOTTOM,
        tarch::la::Vector<3,int>(0,0,_domainSizeZ), tarch::la::Vector<3,int>(_domainSizeX+2,_domainSizeY+2,_domainSizeZ+1),
        tarch::la::Vector<3,int>(0,0,0), tarch::la::Vector<3,int>(_domainSizeX+2,_domainSizeY+2,1));
      #endif
    }

    /** returns true, if this rank is not of relevance for the LB simulation */
    bool skipRank() const {
      int rank = 0;
      #if (COUPLING_MD_PARALLEL==COUPLING_MD_YES)
      MPI_Comm_rank(MPI_COMM_WORLD,&rank);
      #endif
      return ((unsigned int)rank>_processes[0]*_processes[1]*_processes[2]-1);
    }

    const double _omega; // relaxation frequency
    tarch::la::Vector<3,double> _wallVelocity; // velocity of moving wall of Couette flow
    const int _domainSizeX; // domain size in x-direction
    const int _domainSizeY; // domain size in y-direction
    const int _domainSizeZ; // domain size in z-direction
    const int _avgDomainSizeX; // avg. domain size in MPI-parallel simulation in x-direction
    const int _avgDomainSizeY; // "" in y-direction
    const int _avgDomainSizeZ; // "" in z-direction
    double _time; // simulation time
    int _counter; // time step counter
    const double _dx; // actual mesh size
    const double _dt; // actual time step size
    const int _plotEveryTimestep; // number of time steps between vtk plots
    const std::string _filestem; // file stem for vtk plot
    double *_pdf1; // field 1
    double *_pdf2; // field 2
    double *_vel;  // velocity field
    double *_density; // density
    Flag *_flag; // flag field
    #if (COUPLING_MD_PARALLEL==COUPLING_MD_YES)
    double *_sendBufferX; // buffer to send data from left/right to right/left neighbour
    double *_recvBufferX;
    double *_sendBufferY; // buffer to receive data from from left/right neighbour
    double *_recvBufferY;
    double *_sendBufferZ;
    double *_recvBufferZ;
    #endif
    const int _C[19][3]; // lattice velocities
    const double _W[19]; // lattice weights
    const double _xO; // offset for y-direction (lexicographic grid ordering)
    const double _yO; // offset for z-direction
    tarch::la::Vector<3,unsigned int> _processes; // domain decomposition on MPI rank basis
    tarch::la::Vector<6,int> _parallelNeighbours; // neighbour ranks
    // for coupling with MD
    tarch::la::Vector<3,int> _offset;
    tarch::la::Vector<3,int> _globalNumberMacroscopicCells;
};

#endif // _MOLECULARDYNAMICS_COUPLING_SOLVERS_LBCOUETTESOLVER_H_
