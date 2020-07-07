#include "fvCFD.H"
#include "pisoControl.H"

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

int main(){
  int argc=1;
  char a[]={"test"};
  char* b=&a[0];
  char** argv = &b;

  #include "setRootCaseLists.H"
  #include "createTime.H"
  #include "createMesh.H"

  pisoControl piso(mesh);

  #include "createFields.H"
  #include "initContinuityErrs.H"

  U.boundaryFieldRef()[7][0].x() = 1.0;
  U.boundaryFieldRef()[7][2].x() = 2.0;
  U.boundaryFieldRef()[8][1].x() = 1.0;
  U.boundaryFieldRef()[8][3].x() = 2.0;
  U.boundaryFieldRef()[9][2].x() = 1.0;
  U.boundaryFieldRef()[9][4].x() = 2.0;
  U.boundaryFieldRef()[10][3].x() = 1.0;
  U.boundaryFieldRef()[10][5].x() = 2.0;

  Info << mesh.boundary()["fronti"].Cf() << endl;
  Info << mesh.boundary()["backi"].Cf() << endl;
  Info << mesh.boundary()["lefti"].Cf() << endl;
  Info << mesh.boundary()["righti"].Cf() << endl;
  Info << mesh.boundary()["topi"].Cf() << endl;
  Info << mesh.boundary()["buttomi"].Cf() << endl;
  Info << mesh.boundaryMesh()["buttomi"].localPoints() << endl;

  Info<< "\nStarting time loop\n" << endl;

  // while (runTime.loop())
  // {
  //   runFoam(phi, mesh, U, nu, piso, p, pRefCell, pRefValue, cumulativeContErr, runTime);
  // }

  Info<< "End\n" << endl;
}
