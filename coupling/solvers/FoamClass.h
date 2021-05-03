// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _COUPLING_SOLVERS_IcoFoam_H_
#define _COUPLING_SOLVERS_FoamTest_H_

#include "coupling/solvers/CouetteSolver.h"
// Includes from OpenFOAM fdCFD.H, unnecessary includes are removed
#include "fvc.H"
#include "fvm.H"
#include "constrainHbyA.H"
#include "constrainPressure.H"
#include "adjustPhi.H"
#include "findRefCell.H"
#include "pisoControl.H"

namespace coupling{
  namespace solvers{
    class IcoFoam;
  }
}

/**
 * The code is an evolution of the solver IcoFoam in OpenFOAM(R) 7,
 * where additional functionality for MaMiCO coupling is added.
 * It is an incompressible CFD solver for the Couette scenario.
 * The implementation is for a equidistant mesh.
 * Due to the Couette szenario 12 boundaries for the continuum are assumed,
 * 6 of them to be boundaries with the MD.
 * @author Helene Wittenberg
 */
class coupling::solvers::IcoFoam: public coupling::solvers::AbstractCouetteSolver<3> {
public:
  IcoFoam(int rank, int plotEveryTimestep, double channelheight, std::string dict, std::string folder, tarch::la::Vector<12,
    unsigned int> boundariesWithMD, tarch::la::Vector<3,double> uWall):
  AbstractCouetteSolver<3>(),
  runTime(Foam::Time::controlDictName, dict,folder),
  mesh(Foam::IOobject(Foam::fvMesh::defaultRegion,runTime.timeName(),runTime,Foam::IOobject::MUST_READ)),
  transportProperties(Foam::IOobject("transportProperties",runTime.constant(),mesh, Foam::IOobject::MUST_READ_IF_MODIFIED,Foam::IOobject::NO_WRITE)),
  nu("nu", Foam::dimViscosity, transportProperties.lookup("nu")),
  p(Foam::IOobject("p", runTime.timeName(), mesh, Foam::IOobject::MUST_READ, Foam::IOobject::AUTO_WRITE), mesh),
  U(Foam::IOobject("U", runTime.timeName(), mesh, Foam::IOobject::MUST_READ, Foam::IOobject::AUTO_WRITE), mesh),
  phi(Foam::IOobject("phi", runTime.timeName(), mesh, Foam::IOobject::READ_IF_PRESENT, Foam::IOobject::AUTO_WRITE), Foam::fvc::flux(U)),
  piso(mesh),
  _boundariesWithMD(boundariesWithMD),
  _dx(std::cbrt(Foam::max(mesh.cellVolumes()))),
  _channelheight(channelheight),
  _rank(rank),
  _plotEveryTimestep(plotEveryTimestep),
  _uWall(uWall){
    if(skipRank()){return;}
    Foam::setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);
    mesh.setFluxRequired(p.name());
  }
  virtual ~IcoFoam(){
    if(_boundary2RecvBufferIndicesOuter){ delete [] _boundary2RecvBufferIndicesOuter; _boundary2RecvBufferIndicesOuter=NULL;}
    if(_boundary2RecvBufferIndicesInner){ delete [] _boundary2RecvBufferIndicesInner; _boundary2RecvBufferIndicesInner=NULL;}
    if(_boundaryIndices){ delete [] _boundaryIndices; _boundaryIndices=NULL;}
  }

  // advance the solver dt in time, this may include several timesteps,
  // the sequence is based on the OpenFOAM IcoFoam solver
  void advance(double dt)override{
    if(skipRank()){return;}
    unsigned int number = floor(dt/runTime.deltaTValue()+0.5);
    for(unsigned int i=0; i<number; i++){
      using namespace Foam;
      ++runTime;
      Info<< "Time = " << runTime.timeName() << nl << endl;

      // Momentum predictor
      fvVectorMatrix UEqn(fvm::ddt(U) + fvm::div(phi, U)-  fvm::laplacian(nu, U));
      if (piso.momentumPredictor()){
          solve(UEqn == -fvc::grad(p));
      }
      // --- PISO loop
      while (piso.correct()){
        volScalarField rAU(1.0/UEqn.A());
        volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
        surfaceScalarField phiHbyA("phiHbyA",fvc::flux(HbyA)+fvc::interpolate(rAU)*fvc::ddtCorr(U, phi));
        adjustPhi(phiHbyA, U, p);

        // Update the pressure BCs to ensure flux consistency
        constrainPressure(p, U, phiHbyA, rAU);

        // Non-orthogonal pressure corrector loop
        while (piso.correctNonOrthogonal()){
          // Pressure corrector
          fvScalarMatrix pEqn(fvm::laplacian(rAU, p) == fvc::div(phiHbyA));
          pEqn.setReference(pRefCell, pRefValue);
          pEqn.solve();
          if (piso.finalNonOrthogonalIter()){
            phi = phiHbyA - pEqn.flux();
          }
        }
        U = HbyA - rAU*fvc::grad(p);
        U.correctBoundaryConditions();
      }
      // runTime.write(); // writes the original OpenFOAM output
      plottxt();
      _timestepCounter++;
    }
  }

  const double getError()const override{
    double actualError;
    double error = 0.0;
    double errorRMS = 0.0;
    double maxError = 0.0;
    const double pi = 3.141592653589793238;
    const unsigned int s = static_cast<unsigned int>(U.mesh().nCells());
    for(unsigned int i = 0; i< s; i++){
      double u_analytical = 0.0;
      const double pos = U.mesh().C()[i][2];
      for (int k = 1; k < 30; k++){
        u_analytical += 1.0/k * std::sin(k*pi*pos/_channelheight) * std::exp(-k*k * pi*pi/(_channelheight*_channelheight) * nu.value() * runTime.timeOutputValue());
      }
      u_analytical = _uWall[0]*(1.0-pos/_channelheight - 2.0/pi * u_analytical);
      actualError = std::abs(u_analytical-static_cast<double>(U[i][0]));
      error += actualError;
      errorRMS += actualError*actualError;
      maxError = actualError>maxError? actualError: maxError;
    }
    error /= s;
    errorRMS = std::sqrt(errorRMS/s);
    std::string filename = "couetteError_CFD.txt";
    std::ofstream file(filename, std::ios_base::app );
    if (!file.is_open()){std::cout << "ERROR CouetteTest::write2CSV(): Could not open file " << filename << "!" << std::endl; exit(EXIT_FAILURE);}
    file << _timestepCounter << " " << error << " "<< errorRMS << " " << maxError<< std::endl;
    file.close();
    return error;
  }

  tarch::la::Vector<3,double> getVelocity(tarch::la::Vector<3,double> pos)const override{
    const Foam::vector foamPosition(pos[0],pos[1],pos[2]);
    int foamIndice = U.mesh().findCell(foamPosition);
    if (foamIndice > 0){
      return tarch::la::Vector<3,double>(U[foamIndice][0], U[foamIndice][1], U[foamIndice][2]);
    }
    else{return tarch::la::Vector<3,double>(0, 0, 0);}
  };

  // Changes the velocity on the moving wall (refers to Couette szenario)
  void setWallVelocity(tarch::la::Vector<3,double> wallVelocity)override{
    const unsigned int pointsInBoundary = U.boundaryFieldRef()[0].size();
    for(unsigned int i=0; i<pointsInBoundary; i++){
      U.boundaryFieldRef()[0][i].x()=wallVelocity[0];
      U.boundaryFieldRef()[0][i].y()=wallVelocity[1];
      U.boundaryFieldRef()[0][i].z()=wallVelocity[2];
    }
  };

  // Gets the next cell center beside the boundary, necessary to set the boundary condition from MD data
  const tarch::la::Vector<3,double> getOuterPointFromBoundary(const int layer, const int index){
     const Foam::vectorField FoamCoord = U.boundaryFieldRef()[layer].patch().Cf()[index]+(U.boundaryFieldRef()[layer].patch().nf()*_dx*0.5);
     const tarch::la::Vector<3,double> FoamCoordVector(FoamCoord[0][0],FoamCoord[0][1],FoamCoord[0][2]);
     return FoamCoordVector;
  }

  // Gets the next cell center beside the boundary, necessary to set the boundary condition from MD data
  const tarch::la::Vector<3,double> getInnerPointFromBoundary(const int layer, const int index){
     const Foam::vectorField FoamCoord = U.boundaryFieldRef()[layer].patch().Cf()[index]-(U.boundaryFieldRef()[layer].patch().nf()*_dx*0.5);
     const tarch::la::Vector<3,double> FoamCoordVector(FoamCoord[0][0],FoamCoord[0][1],FoamCoord[0][2]);
     return FoamCoordVector;
  }

  // Applies the MD data (just velocities) as boundary condition, the mapping between the conntinuum and the MD is provided by the setMDBoundary()
  void setMDBoundaryValues(std::vector<coupling::datastructures::MacroscopicCell<3>* >& recvBuffer,
  const unsigned int * const recvIndices, const coupling::IndexConversion<3>& indexConversion){
    if(skipRank()){return;}
    for(unsigned int i=0; i < _numberBoundaryPoints; i++){
      unsigned int outer = _boundary2RecvBufferIndicesOuter[i];
      unsigned int inner = _boundary2RecvBufferIndicesInner[i];
      tarch::la::Vector<3,double> localOuterVel( (1.0/recvBuffer[outer]->getMacroscopicMass())*recvBuffer[outer]->getMacroscopicMomentum() );
      tarch::la::Vector<3,double> localInnerVel( (1.0/recvBuffer[inner]->getMacroscopicMass())*recvBuffer[inner]->getMacroscopicMomentum() );
      _boundaryIndices[i]->x() = (localOuterVel[0]+localInnerVel[0])*0.5;
      _boundaryIndices[i]->y() = (localOuterVel[1]+localInnerVel[1])*0.5;
      _boundaryIndices[i]->z() = (localOuterVel[2]+localInnerVel[2])*0.5;
    }
  }

  // Settup for the mapping from MD to continuum data. Velocity values are necessary directly on the boundary. Therefore the MD values from the two cells beside
  // the boundary are interpolated. The function looks for the index of every cell that data is necessary from. This indices will be stored in two arrays.
  void setMDBoundary(tarch::la::Vector<3,double> mdDomainOffset,tarch::la::Vector<3,double> mdDomainSize,unsigned int overlapStrip,
  const coupling::IndexConversion<3>& indexConversion, const unsigned int* const recvIndice, unsigned int size){
    if(skipRank()){return;}
    unsigned int innerMDBoundaryIndex=0;
    while (_boundariesWithMD[innerMDBoundaryIndex] == 0){innerMDBoundaryIndex++;}
    _numberBoundaryPoints = 6*U.boundaryFieldRef()[innerMDBoundaryIndex].size();
    _boundary2RecvBufferIndicesOuter = new unsigned int [_numberBoundaryPoints];
    _boundary2RecvBufferIndicesInner = new unsigned int [_numberBoundaryPoints];
    _boundaryIndices = new Foam::vector* [_numberBoundaryPoints];
    unsigned int counter = 0;
    for (unsigned int boundary=0; boundary < 12; boundary++){
      if(_boundariesWithMD[boundary]==1){
        unsigned int MDPointsPerBoundary=_numberBoundaryPoints/6;
        for (unsigned int j = 0; j < MDPointsPerBoundary; j++){
          _boundaryIndices[counter] = &(U.boundaryFieldRef()[boundary][j]);
          const unsigned int globalIndexOuter = indexConversion.getGlobalCellIndex(indexConversion.getGlobalVectorCellIndex(getOuterPointFromBoundary(boundary, j)));
          const unsigned int globalIndexInner = indexConversion.getGlobalCellIndex(indexConversion.getGlobalVectorCellIndex(getInnerPointFromBoundary(boundary, j)));
          for(unsigned int k = 0; k < size; k++){
            if(globalIndexOuter==recvIndice[k]){
              _boundary2RecvBufferIndicesOuter[counter] = k;
              goto endloop;
            }
          }
          std::cout << "IcoFoam: Within the mapping of the FoamBoundary and the SimpleMD cells there was an error" << std::endl;
          endloop:
          for(unsigned int k = 0; k < size; k++){
            if(globalIndexInner==recvIndice[k]){
              _boundary2RecvBufferIndicesInner[counter] = k;
              goto endloop2;
            }
          }
          std::cout << "IcoFoam: Within the mapping of the FoamBoundary and the SimpleMD cells there was an error" << std::endl;
          endloop2:
          counter++;
        }
      }
    }
  }

private:
  /** create txt plot if required */
  void plottxt() {
    if(_plotEveryTimestep < 1 || _timestepCounter % _plotEveryTimestep > 0) return;
    std::stringstream ss; ss << "velocity_" << _timestepCounter << ".txt";
    std::ofstream file(ss.str().c_str());
    if (!file.is_open()){std::cout << "ERROR NumericalSolver::plottxt(): Could not open file " << ss.str() << "!" << std::endl; exit(EXIT_FAILURE);}
    std::stringstream velocity;

    // loop over domain (incl. boundary)
    double y=_channelheight/2;
    double x=_channelheight/2;
    for (double z = _dx/2; z < _channelheight; z=z+_dx){
      const int foamIndice = U.mesh().findCell(Foam::vector(x,y,z));
      // write information to streams
      if(foamIndice>0){
      velocity << U[foamIndice][0] << ", " << U[foamIndice][1] << ", " << U[foamIndice][2] << std::endl;
    }}
    file << velocity.str() << std::endl;
    file.close();
  }

  // /** create vtk plot if required */
  // void plot() const {
  //   // only plot output if this is the correct timestep
  //   if (_plotEveryTimestep==-1){ return;}
  //   if (_timestepCounter%_plotEveryTimestep!=0){return;}
  //
  //   std::stringstream ss; ss << "Continuum_Velocity_IcoFoam_" << _rank << "_" << _timestepCounter << ".vtk";
  //   std::ofstream file(ss.str().c_str());
  //   if (!file.is_open()){std::cout << "ERROR NumericalSolver::plot(): Could not open file " << ss.str() << "!" << std::endl; exit(EXIT_FAILURE);}
  //   std::stringstream velocity;
  //
  //   file << "# vtk DataFile Version 2.0" << std::endl;
  //   file << "MaMiCo FoamSolver" << std::endl;
  //   file << "ASCII" << std::endl << std::endl;
  //   file << "DATASET STRUCTURED_GRID" << std::endl; It is not a structured grid, the md domain in the middle is missing. ToDo check which other type is fitting
  //   int pointsPerDimension = _channelheight/_dx+2;
  //   file << "DIMENSIONS " << pointsPerDimension << " " << pointsPerDimension << " " << pointsPerDimension << std::endl; // everything +1 cause of change in index
  //   file << "POINTS " << (pointsPerDimension)*(pointsPerDimension)*(pointsPerDimension) << " float" << std::endl;
  //
  //   velocity << std::setprecision(12);
  //   velocity << "VECTORS velocity float" << std::endl;
  //
  //   // // loop over domain (incl. boundary)
  //   // int size = U.size();
  //   // for (int i = 0; i < size; i++){
  //   //   // write information to streams;
  //   //   file << U.mesh()[i][0] << ", " << U.mesh()[i][1] << ", " << U.mesh()[i][2] << std::endl; // boundary is missing, how to include? Just not do?
  //   //   velocity << U[i][0] << ", " << U[i][1] << ", " << U[i][2] << std::endl;
  //   // }
  //
  //   file << std::endl;
  //   file << velocity.str() << std::endl;
  //   velocity.str("");
  //   file.close();
  // }

  // The solver runs sequentially on rank 0. Therefore the function checks if the acutal rank is zero.
  bool skipRank(){
    return !(_rank==0);
  }

  // the following are original OpenFOAM variables, their names shall not be changed
  Foam::Time runTime;
  Foam::fvMesh mesh;
  Foam::IOdictionary transportProperties;
  Foam::dimensionedScalar nu;
  Foam::volScalarField p;
  Foam::volVectorField U;
  Foam::surfaceScalarField phi;
  Foam::pisoControl piso;
  // this are additional variables, they can be changed
  // the entries define which boundaries are for coupling with MD
  // 0 means no MD boundary and 1 means MD boundary
  tarch::la::Vector<12, unsigned int> _boundariesWithMD;
  float _dx; // mesh size
  double _channelheight; // overall height of the Couette channel
  unsigned int *_boundary2RecvBufferIndicesOuter; // pointer to an array with data for communication
  unsigned int *_boundary2RecvBufferIndicesInner; // pointer to an array with data for communication
  Foam::vector **_boundaryIndices; // pointer to OpenFOAM data for communication
  int _rank; // rank of the actual process
  int _plotEveryTimestep; // every n-th time step should be plotted
  int _timestepCounter{0}; // actual time step number
  // the following are original OpenFOAM variables, their names shall not be changed
  Foam::label pRefCell{0};
  Foam::scalar pRefValue{0.0};
  Foam::scalar cumulativeContErr{0};
  unsigned int _numberBoundaryPoints; // the number of CFD boundary points which need data from the MD
  tarch::la::Vector<3,double> _uWall;
};
#endif // _COUPLING_SOLVERS_ICOFOAM_H
