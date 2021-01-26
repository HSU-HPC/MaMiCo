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
 * TODO: sendBuffer request too expensive? Change hard coded variables.
 * The code is an evolution of the solver icoFoam in OpenFOAM(R) 7,
 * where additional functionality for MaMiCO coupling is added.
 * It is an incompressible CFD solver for the Couette scenario.
 * @author Helene Wittenberg
 */
class coupling::solvers::IcoFoam: public coupling::solvers::AbstractCouetteSolver<3> {
public:
  IcoFoam(int rank, int plotEveryTimestep, double channelheight, std::string dict, std::string folder, tarch::la::Vector<12, unsigned int> boundariesWithMD):
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
  _boundary2RecvBufferIndicesOuter(new unsigned int [_numberBoundaryPoints]),
  _boundary2RecvBufferIndicesInner(new unsigned int [_numberBoundaryPoints]),
  _boundaryIndices(new Foam::vector* [_numberBoundaryPoints]),
  _rank(rank),
    _plotEveryTimestep(plotEveryTimestep){
    if(skipRank()){return;}
    unsigned int innerMDBoundaryIndex=0;
    while (_boundariesWithMD[innerMDBoundaryIndex] == 0){innerMDBoundaryIndex++;}
    _numberBoundaryPoints = 6*U.boundaryFieldRef()[innerMDBoundaryIndex].size();
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
      // runTime.write();
      plottxt();
      _timestepCounter++;
    }
  }

  tarch::la::Vector<3,double> getVelocity(tarch::la::Vector<3,double> pos)const override{
    const Foam::vector foamPosition(pos[0],pos[1],pos[2]);
    int foamIndice = U.mesh().findCell(foamPosition);
    if (foamIndice > 0){
      return tarch::la::Vector<3,double>(U[foamIndice][0], U[foamIndice][1], U[foamIndice][2]);
    }
    else{return tarch::la::Vector<3,double>(0, 0, 0);}
  };

  void setWallVelocity(tarch::la::Vector<3,double> wallVelocity)override{
    const unsigned int pointsInBoundary = U.boundaryFieldRef()[0].size();
    for(unsigned int i=0; i<pointsInBoundary; i++){
      U.boundaryFieldRef()[0][i].x()=wallVelocity[0];
      U.boundaryFieldRef()[0][i].y()=wallVelocity[1];
      U.boundaryFieldRef()[0][i].z()=wallVelocity[2];
    }
  };

  const tarch::la::Vector<3,double> getOuterPointFromBoundary(const int layer, const int index){
     const Foam::vectorField FoamCoord = U.boundaryFieldRef()[layer].patch().Cf()[index]+(U.boundaryFieldRef()[layer].patch().nf()*_dx*0.5);
     const tarch::la::Vector<3,double> FoamCoordVector(FoamCoord[0][0],FoamCoord[0][1],FoamCoord[0][2]);
     return FoamCoordVector;
  }

  const tarch::la::Vector<3,double> getInnerPointFromBoundary(const int layer, const int index){
     const Foam::vectorField FoamCoord = U.boundaryFieldRef()[layer].patch().Cf()[index]-(U.boundaryFieldRef()[layer].patch().nf()*_dx*0.5);
     const tarch::la::Vector<3,double> FoamCoordVector(FoamCoord[0][0],FoamCoord[0][1],FoamCoord[0][2]);
     return FoamCoordVector;
  }

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

  void setMDBoundary(tarch::la::Vector<3,double> mdDomainOffset,tarch::la::Vector<3,double> mdDomainSize,unsigned int overlapStrip,
  const coupling::IndexConversion<3>& indexConversion, const unsigned int* const recvIndice, unsigned int size){
    if(skipRank()){return;}
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
  tarch::la::Vector<12, unsigned int> _boundariesWithMD;
  float _dx;
  double _channelheight;
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
};
#endif // _COUPLING_SOLVERS_ICOFOAM_H
