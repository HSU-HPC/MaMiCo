// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _COUPLING_SOLVERS_RHOCENTRALINTERFACE4EVAPORATION_H_
#define _COUPLING_SOLVERS_RHOCENTRALINTERFACE4EVAPORATION_H_

// Includes necessary for OpenFOAM
#ifndef namespaceFoam
#define namespaceFoam
#endif
#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "psiThermo.H"
#include "turbulentFluidThermoModel.H"
#include "fixedRhoFvPatchScalarField.H"
#include "directionInterpolate.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"

namespace coupling {
namespace solvers {
class RhoCentralInterface4Evaporation;
}
}

/** interface to run rhoCentralFoam from the OpenFOAM library as macro solver for the Evaporation
 *  scenario
 *  @author Helene Wittenberg
 */
class coupling::solvers::RhoCentralInterface4Evaporation {
public:
  /* Constructs an instance of the class and initialises all the necessary variables for OpenFOAM
   * The initialisation is mostly copied from the following files, some things are sligtly adapted
   * "createTime.H", "createDynamicFvMesh.H", "createFields.H", "createFieldRefs.H", "createTimeControls.H"
   * and "readFluxScheme.H"
   * @param directory path to the openFoam setup folder
   * @param folder name of the openFoam setup folder
   * @param rank rank of the current proccess
   */
  RhoCentralInterface4Evaporation(std::string directory, std::string folder, const int rank)
  :
  runTime(directory, folder),
  meshPtr(Foam::dynamicFvMesh::New(Foam::IOobject(Foam::fvMesh::defaultRegion, runTime.timeName(), runTime, Foam::IOobject::MUST_READ))),
  mesh(meshPtr()),
  adjustTimeStep(runTime.controlDict().getOrDefault("adjustTimeStep", false)),
  maxCo(runTime.controlDict().getOrDefault<Foam::scalar>("maxCo", 1)),
  maxDeltaT(runTime.controlDict().getOrDefault<Foam::scalar>("maxDeltaT", Foam::GREAT)),
  pThermo(Foam::psiThermo::New(mesh)),
  thermo(pThermo()),
  e(thermo.he()),
  U(Foam::IOobject("U", runTime.timeName(), mesh, Foam::IOobject::MUST_READ, Foam::IOobject::AUTO_WRITE), mesh),
  rho(Foam::IOobject("rho", runTime.timeName(), mesh, Foam::IOobject::NO_READ, Foam::IOobject::AUTO_WRITE), thermo.rho()),
  rhoU(Foam::IOobject("rhoU", runTime.timeName(), mesh, Foam::IOobject::NO_READ, Foam::IOobject::NO_WRITE ), rho*U),
  rhoE(Foam::IOobject("rhoE", runTime.timeName(), mesh, Foam::IOobject::NO_READ, Foam::IOobject::NO_WRITE ), rho*(e + 0.5*magSqr(U))),
  pos(Foam::IOobject("pos", runTime.timeName(), mesh), mesh, Foam::dimensionedScalar("pos", Foam::dimless, 1.0)),
  neg(Foam::IOobject("neg", runTime.timeName(), mesh), mesh, Foam::dimensionedScalar("neg", Foam::dimless, -1.0)),
  phi("phi", Foam::fvc::flux(rhoU)),
  turbulence(Foam::compressible::turbulenceModel::New(rho, U, phi, thermo)),
  v_zero(Foam::dimVolume/Foam::dimTime, Foam::Zero),
  p(thermo.p()),
  T(thermo.T()),
  psi(thermo.psi()),
  mu(thermo.mu()),
  _rank(rank)
  {
    if (skipRank()) { return; }
    initialiseOpenfoam();
  }


  int initialiseOpenfoam(){
    using namespace Foam;
    int argc{1};
    char* dummy_args[] = { "rhoCentralFoam", NULL };
    char** argv = dummy_args;
    #include "postProcess.H"
    #include "addCheckCaseOptions.H"
    // #include "setRootCaseLists.H" // No initialisation
    turbulence->validate();
    if (max(mu.primitiveField()) > 0.0) // copied from #include "createFieldRefs.H"
    {
      inviscid = false;
    }
    if (mesh.schemesDict().readIfPresent("fluxScheme", fluxScheme)) // copied from #include "readFluxScheme.H"
    {
      if (!(fluxScheme == "Tadmor") && !(fluxScheme == "Kurganov"))
      {
        FatalErrorInFunction << "fluxScheme: " << fluxScheme << " is not a valid choice. "
          << "Options are: Tadmor, Kurganov" << abort(FatalError);
      }
    }
    return 0;
  }

  virtual ~RhoCentralInterface4Evaporation() {
    if (skipRank()) {
      return;
    }
  }

  /* @brief run one time step and advance all the fields */
  void advance() {
    using namespace Foam;
    #include "readTimeControls.H"

    if (!LTS)
    {
      #include "setDeltaT.H"
      ++runTime;
      // Do any mesh changes
      mesh.update();
    }

    // --- Directed interpolation of primitive fields onto faces
    surfaceScalarField rho_pos(interpolate(rho, pos));
    surfaceScalarField rho_neg(interpolate(rho, neg));
    surfaceVectorField rhoU_pos(interpolate(rhoU, pos, U.name()));
    surfaceVectorField rhoU_neg(interpolate(rhoU, neg, U.name()));
    volScalarField rPsi("rPsi", 1.0/psi);
    surfaceScalarField rPsi_pos(interpolate(rPsi, pos, T.name()));
    surfaceScalarField rPsi_neg(interpolate(rPsi, neg, T.name()));
    surfaceScalarField e_pos(interpolate(e, pos, T.name()));
    surfaceScalarField e_neg(interpolate(e, neg, T.name()));
    surfaceVectorField U_pos("U_pos", rhoU_pos/rho_pos);
    surfaceVectorField U_neg("U_neg", rhoU_neg/rho_neg);
    surfaceScalarField p_pos("p_pos", rho_pos*rPsi_pos);
    surfaceScalarField p_neg("p_neg", rho_neg*rPsi_neg);
    surfaceScalarField phiv_pos("phiv_pos", U_pos & mesh.Sf());

    // Note: extracted out the orientation so becomes unoriented
    phiv_pos.setOriented(false);
    surfaceScalarField phiv_neg("phiv_neg", U_neg & mesh.Sf());
    phiv_neg.setOriented(false);

    // Make fluxes relative to mesh-motion
    if (mesh.moving())
    {
      surfaceScalarField meshPhi(mesh.phi());
      meshPhi.setOriented(false);
      phiv_pos -= meshPhi;
      phiv_neg -= meshPhi;
    }

    volScalarField c("c", sqrt(thermo.Cp()/thermo.Cv()*rPsi));

    surfaceScalarField cSf_pos("cSf_pos", interpolate(c, pos, T.name())*mesh.magSf());

    surfaceScalarField cSf_neg("cSf_neg", interpolate(c, neg, T.name())*mesh.magSf());

    surfaceScalarField ap("ap", max(max(phiv_pos + cSf_pos, phiv_neg + cSf_neg), v_zero));

    surfaceScalarField am("am",min(min(phiv_pos - cSf_pos, phiv_neg - cSf_neg), v_zero));

    surfaceScalarField a_pos("a_pos", ap/(ap - am));

    surfaceScalarField amaxSf("amaxSf", max(mag(am), mag(ap)));

    surfaceScalarField aSf("aSf", am*a_pos);

    if (fluxScheme == "Tadmor")
    {
      aSf = -0.5*amaxSf;
      a_pos = 0.5;
    }

    surfaceScalarField a_neg("a_neg", 1.0 - a_pos);

    phiv_pos *= a_pos;
    phiv_neg *= a_neg;

    surfaceScalarField aphiv_pos("aphiv_pos", phiv_pos - aSf);
    surfaceScalarField aphiv_neg("aphiv_neg", phiv_neg + aSf);

    // Reuse amaxSf for the maximum positive and negative fluxes estimated by the central scheme
    amaxSf = max(mag(aphiv_pos), mag(aphiv_neg));
    #include "centralCourantNo.H"

    if (LTS)
    {
      #include "setRDeltaT.H"
      ++runTime;
    }

    Info<< "Time = " << runTime.timeName() << nl << endl;

    phi = aphiv_pos*rho_pos + aphiv_neg*rho_neg;

    surfaceVectorField phiU(aphiv_pos*rhoU_pos + aphiv_neg*rhoU_neg);
    // Note: reassembled orientation from the pos and neg parts so becomes oriented
    phiU.setOriented(true);

    surfaceVectorField phiUp(phiU + (a_pos*p_pos + a_neg*p_neg)*mesh.Sf());
    surfaceScalarField phiEp("phiEp", aphiv_pos*(rho_pos*(e_pos + 0.5*magSqr(U_pos)) + p_pos)
      + aphiv_neg*(rho_neg*(e_neg + 0.5*magSqr(U_neg)) + p_neg) + aSf*p_pos - aSf*p_neg);

    // Make flux for pressure-work absolute
    if (mesh.moving())
    {
      surfaceScalarField meshPhi(mesh.phi());
      meshPhi.setOriented(false);
      phiEp += meshPhi*(a_pos*p_pos + a_neg*p_neg);
    }

    volScalarField muEff("muEff", turbulence->muEff());
    volTensorField tauMC("tauMC", muEff*dev2(Foam::T(fvc::grad(U))));

    // --- Solve density
    solve(fvm::ddt(rho) + fvc::div(phi));

    // --- Solve momentum
    solve(fvm::ddt(rhoU) + fvc::div(phiUp));

    U.ref() = rhoU()/rho();
    U.correctBoundaryConditions();
    rhoU.boundaryFieldRef() == rho.boundaryField()*U.boundaryField();

    if (!inviscid)
    {
      solve(fvm::ddt(rho, U) - fvc::ddt(rho, U) - fvm::laplacian(muEff, U) - fvc::div(tauMC));
      rhoU = rho*U;
    }

    // --- Solve energy
    surfaceScalarField sigmaDotU("sigmaDotU", (fvc::interpolate(muEff)*mesh.magSf()*fvc::snGrad(U)
      + fvc::dotInterpolate(mesh.Sf(), tauMC)) & (a_pos*U_pos + a_neg*U_neg));
    solve(fvm::ddt(rhoE) + fvc::div(phiEp) - fvc::div(sigmaDotU));
    e = rhoE/rho - 0.5*magSqr(U);
    e.correctBoundaryConditions();
    thermo.correct();
    rhoE.boundaryFieldRef() == rho.boundaryField()*(e.boundaryField() + 0.5*magSqr(U.boundaryField()));

    if (!inviscid)
    {
      solve(fvm::ddt(rho, e) - fvc::ddt(rho, e) - fvm::laplacian(turbulence->alphaEff(), e));
      thermo.correct();
      rhoE = rho*(e + 0.5*magSqr(U));
    }

    p.ref() = rho()/psi();
    p.correctBoundaryConditions();
    rho.boundaryFieldRef() == psi.boundaryField()*p.boundaryField();

    turbulence->correct();

    runTime.write();

    runTime.printExecutionTime(Info);
  }

  /* @brief returns the velocity at a position for the current time step
   * @param pos vector of the position to get the velocity */
  tarch::la::Vector<3, double> getVelocity(tarch::la::Vector<3, double> pos) const {
    using namespace Foam;
    const Foam::vector foamPosition(pos[0], pos[1], pos[2]);
    const int foamIndice = U.mesh().findCell(foamPosition);
    if (foamIndice > 0) {
      Info << "velocity " ;
      Info << U[foamIndice] << endl;
      return tarch::la::Vector<3, double>(U[foamIndice][0], U[foamIndice][1], U[foamIndice][2]);
    } else {
      std::cout << "ERROR EvaporationTest::RhoCentralInterface4Evaporation::getVelocity(): pos " <<
      pos[0] << " " << pos[1] << " " << pos[2] <<" can't be found in the OpenFOAM mesh" << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  /* @brief returns the density at a position for the current time step
   * @param pos vector of the position to get the density */
  double getDensity(tarch::la::Vector<3, double> pos) const {
    using namespace Foam;
    const Foam::vector foamPosition(pos[0], pos[1], pos[2]);
    Info << "vector ";
    Info << foamPosition << endl;
    const int foamIndice = rho.mesh().findCell(foamPosition);
    // Foam::vector tmp = rho.mesh()[foamIndice];
    Info << "findCell ";
    Info << rho.mesh().C()[foamIndice] << endl;
    if (foamIndice > 0) {
      Info << "density " << rho[foamIndice] << endl;
      return rho[foamIndice];

    } else {
      std::cout << "ERROR EvaporationTest::RhoCentralInterface4Evaporation::getDensity(): pos " <<
      pos[0] << " " << pos[1] << " " << pos[2] <<" can't be found in the OpenFOAM mesh" << std::endl;
      exit(EXIT_FAILURE);
    }
  }

private:

  void plottxt() const {
    // if (_plotEveryTimestep < 1 || _timestepCounter % _plotEveryTimestep > 0)
    //   return;
    // std::stringstream ss;
    // ss << "velocity_" << _timestepCounter << ".txt";
    // std::ofstream file(ss.str().c_str());
    // if (!file.is_open()) {
    //   std::cout << "ERROR NumericalSolver::plottxt(): Could not open file " << ss.str() << "!" << std::endl;
    //   exit(EXIT_FAILURE);
    // }
    // std::stringstream velocity;
    //
    // // loop over domain (incl. boundary)
    // double y = _channelheight / 2;
    // double x = _channelheight / 2;
    // for (double z = _dx / 2; z < _channelheight; z = z + _dx) {
    //   const int foamIndice = U.mesh().findCell(Foam::vector(x, y, z));
    //   // write information to streams
    //   if (foamIndice > 0) {
    //     velocity << U[foamIndice][0] << ", " << U[foamIndice][1] << ", " << U[foamIndice][2] << std::endl;
    //   }
    // }
    // file << velocity.str() << std::endl;
    // file.close();
  }

  // The solver runs sequentially on rank 0. Therefore the function checks if the acutal rank is zero.
  bool skipRank() const { return !(_rank == 0); }

  /** @brief */
  Foam::Time runTime;
  /** @brief */
  Foam::autoPtr<Foam::dynamicFvMesh> meshPtr;
  /** @brief */
  Foam::dynamicFvMesh& mesh;
  /** @brief */
  bool adjustTimeStep;
  /** @brief */
  Foam::scalar maxCo;
  /** @brief */
  Foam::scalar maxDeltaT;
  /** @brief */
  Foam::autoPtr<Foam::psiThermo> pThermo;
  /** @brief */
  Foam::psiThermo& thermo;
  /** @brief */
  Foam::volScalarField& e;
  /** @brief */
  Foam::volVectorField U;
  /** @brief */
  Foam::volScalarField rho;
  /** @brief */
  Foam::volVectorField rhoU;
  /** @brief */
  Foam::volScalarField rhoE;
  /** @brief */
  Foam::surfaceScalarField pos;
  /** @brief */
  Foam::surfaceScalarField neg;
  /** @brief */
  Foam::surfaceScalarField phi;
  /** @brief */
  Foam::autoPtr<Foam::compressible::turbulenceModel> turbulence;
  /** @brief */
  Foam::tmp<Foam::volScalarField> trDeltaT;
  /** @brief */
  Foam::dimensionedScalar v_zero;
  /** @brief */
  Foam::volScalarField& p;
  /** @brief */
  const Foam::volScalarField& T;
  /** @brief */
  const Foam::volScalarField& psi;
  /** @brief */
  const Foam::volScalarField& mu;
  /** @brief */
  Foam::word fluxScheme{"Kurganov"};
  /** @brief */
  bool inviscid{true};
  /** @brief */
  Foam::scalar CoNum{0.0};
  /** @brief */
  Foam::scalar meanCoNum{0.0};
  /** @brief */
  const bool LTS{false};
  /** @brief */
  const int _rank{0};
  /** @brief */
  int _plotEveryTimestep{0};
  /** @brief */
  int _timestepCounter{0};
  /** @brief */
  double _channelheight{200};
  /** @brief */
  double _dx{2.5};
};

#endif // _COUPLING_SOLVERS_RHOCENTRALINTERFACE4EVAPORATION_H_
