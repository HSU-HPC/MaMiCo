// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _COUPLING_SOLVERS_IcoFoam_H_
#define _COUPLING_SOLVERS_FoamTest_H_

// Includes from OpenFOAM fdCFD.H
// #include "parRun.H"
#include "Time.H"
#include "fvMesh.H"
#include "fvc.H"
// #include "fvMatrices.H"
#include "fvm.H"
// #include "linear.H"
#include "uniformDimensionedFields.H"
#include "calculatedFvPatchFields.H"
#include "extrapolatedCalculatedFvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "zeroGradientFvPatchFields.H"
#include "fixedFluxPressureFvPatchScalarField.H"
#include "constrainHbyA.H"
#include "constrainPressure.H"
#include "adjustPhi.H"
#include "findRefCell.H"
// #include "IOMRFZoneList.H"
// #include "constants.H"
#include "OSspecific.H"
#include "argList.H"
#include "timeSelector.H"
#include "pisoControl.H"

namespace coupling{
  namespace solvers{
    class IcoFoam;
  }
}

/**
 *
 *  @author Helene Wittenberg
 */
class coupling::Solvers::IcoFoam: {
public:
  IcoFoam(): {
    Foam::setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);
    mesh.setFluxRequired(p.name());
  }
  virtual ~IcoFoam(){}

  void advance(){
    using namespace Foam;
    ++runTime;
    Info<< "Time = " << runTime.timeName() << nl << endl;

    // scalar CoNum = 0.0;
    // scalar meanCoNum = 0.0;
    //
    // {
    //     scalarField sumPhi
    //     (
    //         fvc::surfaceSum(mag(phi))().primitiveField()
    //     );
    //
    //     CoNum = 0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();
    //
    //     meanCoNum =
    //         0.5*(gSum(sumPhi)/gSum(mesh.V().field()))*runTime.deltaTValue();
    // }
    //
    // Info<< "Courant Number mean: " << meanCoNum
    //     << " max: " << CoNum << endl;

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
      // {
      //     volScalarField contErr(fvc::div(phi));
      //
      //     scalar sumLocalContErr = runTime.deltaTValue()*mag(contErr)().weightedAverage(mesh.V()).value();
      //
      //     scalar globalContErr = runTime.deltaTValue()*contErr.weightedAverage(mesh.V()).value();
      //     cumulativeContErr += globalContErr;
      //
      //     // Info<< "time step continuity errors : sum local = " << sumLocalContErr
      //     //     << ", global = " << globalContErr
      //     //     << ", cumulative = " << cumulativeContErr
      //     //     << endl;
      // }
      U = HbyA - rAU*fvc::grad(p);
      U.correctBoundaryConditions();
    }
    //runTime.write();
    plottxt();
  }

  const tarch::la::Vector<,double> getPointFromBoundary(const int layer, const int index)const{
     const Foam::vectorField FoamCoord = U.boundaryFieldRef()[layer].patch().Cf()[index]+(U.boundaryFieldRef()[layer].patch().nf()*1.25);
     const tarch::la::Vector<3,double> FoamCoordVector(FoamCoord[0][0],FoamCoord[0][1],FoamCoord[0][2]))
     return FoamCoordVector;
  }

  const Foam::vector** getReference4BoundaryField(const int layer, const int index)const{
    return &(U.boundaryFieldRef()[layer][index]);
  }

  unsigned int getIndex4Vector(const Foam::vector Vector)const{
    return U.mesh().findCell(foamPosition) > 0 ? mesh.findCell(foamPosition) : 0;
  }

  tarch::la::Vector<3,double> getVelocityByIndex(const unsigned int index){
    return tarch::la::Vector<3,double>(U[index][0], U[index][1], U[index][2]);
  }

private:
  /** create vtk plot if required */
  void plottxt() {
    std::stringstream ss; ss << "velocity_" << runTime.timeOutputValue()/0.25 << ".txt";
    std::ofstream file(ss.str().c_str());
    if (!file.is_open()){std::cout << "ERROR NumericalSolver::plottxt(): Could not open file " << ss.str() << "!" << std::endl; exit(EXIT_FAILURE);}
    std::stringstream velocity;

    // loop over domain (incl. boundary)
    double y=25;
    double x=25;
    for (double z = 1.25; z < 50.0; z=z+2.5){
      const int foamIndice = U.mesh().findCell(Foam::vector(x,y,z));
      // write information to streams
      if(foamIndice>0){
      velocity << U[foamIndice][0] << ", " << U[foamIndice][1] << ", " << U[foamIndice][2] << std::endl;
    }}
    file << velocity.str() << std::endl;
    file.close();
  }

  Foam::Time runTime{Foam::Time(Foam::Time::controlDictName, "/home/helene/Dokumente/mamico-dev/coupling/tests", "build_foam")};
  Foam::fvMesh mesh{Foam::fvMesh(Foam::IOobject(Foam::fvMesh::defaultRegion,runTime.timeName(),runTime,Foam::IOobject::MUST_READ))};
  Foam::IOdictionary transportProperties{Foam::IOdictionary(Foam::IOobject("transportProperties",runTime.constant(),mesh, Foam::IOobject::MUST_READ_IF_MODIFIED,Foam::IOobject::NO_WRITE))};
  Foam::dimensionedScalar nu{Foam::dimensionedScalar("nu", Foam::dimViscosity, transportProperties.lookup("nu"))};
  Foam::volScalarField p{Foam::volScalarField(Foam::IOobject("p", runTime.timeName(), mesh, Foam::IOobject::MUST_READ, Foam::IOobject::AUTO_WRITE), mesh)};
  Foam::volVectorField U{Foam::volVectorField(Foam::IOobject("U", runTime.timeName(), mesh, Foam::IOobject::MUST_READ, Foam::IOobject::AUTO_WRITE), mesh)};
  Foam::surfaceScalarField phi{Foam::surfaceScalarField(Foam::IOobject("phi", runTime.timeName(), mesh, Foam::IOobject::READ_IF_PRESENT, Foam::IOobject::AUTO_WRITE), Foam::fvc::flux(U))};
  Foam::label pRefCell{0};
  Foam::scalar pRefValue{0.0};
  Foam::pisoControl piso{Foam::pisoControl(mesh)};
  Foam::scalar cumulativeContErr{0};
};
#endif // _COUPLING_SOLVERS_ICOFOAM_H
