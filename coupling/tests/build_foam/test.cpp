#include "fvCFD.H"
#include "pisoControl.H"
#include <vector>

void runFoam(Foam::surfaceScalarField& phi, Foam::fvMesh& mesh,
  Foam::volVectorField& U, Foam::dimensionedScalar& nu, Foam::pisoControl& piso,
  Foam::volScalarField& p, Foam::label& pRefCell, Foam::scalar& pRefValue,
  Foam::scalar& cumulativeContErr, Foam::Time& runTime){
  Info<< "Time = " << runTime.timeName() << nl << endl;

  #include "CourantNo.H"

  // Momentum predictor

  fvVectorMatrix UEqn
  (
      fvm::ddt(U)
    + fvm::div(phi, U)
    - fvm::laplacian(nu, U)
  );

  if (piso.momentumPredictor())
  {
      solve(UEqn == -fvc::grad(p));
  }

  // --- PISO loop
  while (piso.correct())
  {
      volScalarField rAU(1.0/UEqn.A());
      volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
      surfaceScalarField phiHbyA
      (
          "phiHbyA",
          fvc::flux(HbyA)
        + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
      );

      adjustPhi(phiHbyA, U, p);

      // Update the pressure BCs to ensure flux consistency
      constrainPressure(p, U, phiHbyA, rAU);

      // Non-orthogonal pressure corrector loop
      while (piso.correctNonOrthogonal())
      {
          // Pressure corrector

          fvScalarMatrix pEqn
          (
              fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
          );

          pEqn.setReference(pRefCell, pRefValue);

          pEqn.solve();

          if (piso.finalNonOrthogonalIter())
          {
              phi = phiHbyA - pEqn.flux();
          }
      }

      #include "continuityErrs.H"

      U = HbyA - rAU*fvc::grad(p);
      U.correctBoundaryConditions();
  }

  runTime.write();

  Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
      << "  ClockTime = " << runTime.elapsedClockTime() << " s"
      << nl << endl;
}

Foam::vector stdVector2FoamVector(std::vector<double>& stdVector){
  Foam::vector FoamVector(stdVector[0], stdVector[1], stdVector[2]);
  return FoamVector;
}

int main(){
  int argc=1;
  char a[]={"test"};
  char* b=&a[0];
  char** argv = &b;

  //#include "listOptions.H"

  Foam::argList args(argc, argv);
  if (!args.checkRootCase())
  {
      Foam::FatalError.exit();
  }

  //#include "listOutput.H"

  Foam::Time runTime(Foam::Time::controlDictName, args);
  Info << "Path" << endl;
  Info << runTime.path() << endl;
  Info << "ControlDict" << endl;
  Info << runTime.controlDict() << endl;
  Info << "dbDir" << endl;
  Info << runTime.dbDir() << endl;

  Foam::fvMesh mesh(Foam::IOobject(
          Foam::fvMesh::defaultRegion,
          runTime.timeName(),
          runTime,
          Foam::IOobject::MUST_READ
  ));

  pisoControl piso(mesh);

  IOdictionary transportProperties(IOobject(
          "transportProperties",
          runTime.constant(),
          mesh,
          IOobject::MUST_READ_IF_MODIFIED,
          IOobject::NO_WRITE
  ));

  dimensionedScalar nu(
      "nu",
      dimViscosity,
      transportProperties.lookup("nu")
  );
  Info << transportProperties.lookup("nu") << endl;
  Info << runTime.timeName() << endl;

  volScalarField p(IOobject(
          "p",
          runTime.timeName(),
          mesh,
          IOobject::MUST_READ,
          IOobject::AUTO_WRITE
      ), mesh
  );

  volVectorField U(IOobject(
          "U",
          runTime.timeName(),
          mesh,
          IOobject::MUST_READ,
          IOobject::AUTO_WRITE
      ), mesh
  );

  surfaceScalarField phi(IOobject(
          "phi",
          runTime.timeName(),
          mesh,
          IOobject::READ_IF_PRESENT,
          IOobject::AUTO_WRITE
      ), fvc::flux(U)
  );

  label pRefCell = 0;
  scalar pRefValue = 0.0;
  setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);
  mesh.setFluxRequired(p.name());

  #ifndef initContinuityErrs_H
  #define initContinuityErrs_H
  scalar cumulativeContErr = 0;
  #endif

  // U.boundaryFieldRef()[7][0].x() = 1.0;
  // U.boundaryFieldRef()[8][1].x() = 1.0;
  // Info << mesh.boundary()["fronti"].Cf() << endl;
  // Info << U.mesh().C()[0] << endl;
  // Info << U[0] << endl;
  // Info << mesh.C()[0] << endl;
  //
  // std::vector<double> origPoint = {8.5, 9, 9.1};//mesh.C()[10];
  // Foam::vector myPoint(stdVector2FoamVector(origPoint));
  // Info << U.mesh().findCell(myPoint) << endl;

  Info<< "\nStarting time loop\n" << endl;

  while (runTime.loop())
  {
    runFoam(phi, mesh, U, nu, piso, p, pRefCell, pRefValue, cumulativeContErr, runTime);
  }

  Info<< "End\n" << endl;
}
