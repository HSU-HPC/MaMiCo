// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _COUPLING_TESTS_FoamTest_H_
#define _COUPLING_TESTS_FoamTest_H_

#include "tarch/utils/MultiMDService.h"
#include "simplemd/configurations/MolecularDynamicsConfiguration.h"
#include "coupling/tests/Test.h"
#include "coupling/CouplingMDDefinitions.h"
#include "tarch/configuration/ParseConfiguration.h"
#include "coupling/solvers/CouetteSolverInterface.h"
#include "coupling/interface/MDSimulationFactory.h"
#include "coupling/interface/impl/SimpleMD/SimpleMDSolverInterface.h"
#include "coupling/configurations/MaMiCoConfiguration.h"
#include "coupling/services/MultiMDCellService.h"
#if (COUPLING_MD_PARALLEL==COUPLING_MD_YES)
#include <mpi.h>
#endif
#include <sys/time.h>
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

/**
 *
 *  @author Helene Wittenberg
 */
class FoamTest: public Test {
public:
  FoamTest(): Test("FoamTest") {
    Foam::setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);
    mesh.setFluxRequired(p.name());
  }
  virtual ~FoamTest(){}

  virtual void run(){
    init();
    for (int cycle = 0; cycle < _cfg.couplingCycles; cycle++){
      advanceMacro();
      advanceMicro(cycle);
      twoWayCoupling(cycle);}
    shutdown();
  }

private:
  void init(){
    initMPI();
    parseConfigurations();
    initSolvers();
  }

  void initMPI(){
    _rank = 0;
    #if (COUPLING_MD_PARALLEL==COUPLING_MD_YES)
    MPI_Comm_rank(MPI_COMM_WORLD,&_rank);
    #endif
  }

  void parseConfigurations(){
    tarch::configuration::ParseConfiguration::parseConfiguration<coupling::configurations::MaMiCoConfiguration<3> >("couette_mamico.xml","mamico",_mamicoConfig);
    if (!_mamicoConfig.isValid()){ std::cout << "ERROR FoamTest: Invalid MaMiCo config!" << std::endl; exit(EXIT_FAILURE); }

    tarch::configuration::ParseConfiguration::parseConfiguration<simplemd::configurations::MolecularDynamicsConfiguration>("couette_simplemd.xml","molecular-dynamics",_simpleMDConfig);
    if (!_simpleMDConfig.isValid()){std::cout << "ERROR FoamTest: Invalid SimpleMD config!" << std::endl; exit(EXIT_FAILURE);}

    parseFoamTestConfiguration();
  }

  void parseFoamTestConfiguration(){
    tinyxml2::XMLDocument conffile;
    tinyxml2::XMLElement *node = NULL;
    conffile.LoadFile("couette.xml");
    node = conffile.FirstChildElement("couette-test");
    if (node == NULL){
      std::cout << "Could not read input file couette.xml: missing element <couette-test>" << std::endl;
      exit(EXIT_FAILURE);
    }
    tinyxml2::XMLElement *n2 = node->NextSiblingElement();
    if(n2 != NULL){
      std::cout << "Could not read input file couette.xml: unknown element " << n2->Name() << std::endl;
      exit(EXIT_FAILURE);
    }

    tinyxml2::XMLElement* subtag = node->FirstChildElement("coupling");
    if (subtag == NULL){
      std::cout << "Could not read input file couette.xml: Missing subtag: coupling" << std::endl;
      exit(EXIT_FAILURE);
    }
    tarch::configuration::ParseConfiguration::readIntMandatory(_cfg.couplingCycles,subtag,"coupling-cycles");
    tarch::configuration::ParseConfiguration::readBoolMandatory(_cfg.twoWayCoupling,subtag,"two-way-coupling");
    tarch::configuration::ParseConfiguration::readBoolMandatory(_cfg.md2Macro,subtag,"send-from-md-to-macro");
    tarch::configuration::ParseConfiguration::readBoolMandatory(_cfg.macro2Md,subtag,"send-from-macro-to-md");
    tarch::configuration::ParseConfiguration::readIntMandatory(_cfg.csvEveryTimestep,subtag,"write-csv-every-timestep");

    subtag = node->FirstChildElement("microscopic-solver");
    if (subtag == NULL){
      std::cout << "Could not read input file couette.xml: Missing subtag: microscopic-solver" << std::endl;
      exit(EXIT_FAILURE);
    }
    std::string type;
    tarch::configuration::ParseConfiguration::readStringMandatory(type,subtag,"type");
    if(type == "md"){
      tarch::configuration::ParseConfiguration::readDoubleMandatory(_cfg.temp,subtag,"temperature");
      tarch::configuration::ParseConfiguration::readIntMandatory(_cfg.equSteps,subtag,"equilibration-steps");
      tarch::configuration::ParseConfiguration::readIntMandatory(_cfg.totalNumberMDSimulations,subtag,"number-md-simulations");
      if(_cfg.totalNumberMDSimulations < 1){
        std::cout << "Could not read input file couette.xml: number-md-simulations < 1" << std::endl;
        exit(EXIT_FAILURE);
      }
      tarch::configuration::ParseConfiguration::readDoubleMandatory(_cfg.density,subtag,"density");
    }
    else{
      std::cout << "Could not read input file couette.xml: Unknown microscopic solver type!" << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  void initSolvers(){
    // for timing measurements
    _tv.micro = 0;
    _tv.macro = 0;
    _tv.filter = 0;

    _multiMDService = new tarch::utils::MultiMDService<3>(_simpleMDConfig.getMPIConfiguration().getNumberOfProcesses(),_cfg.totalNumberMDSimulations);
    _localMDInstances = _multiMDService->getLocalNumberOfMDSimulations();

    for (unsigned int i = 0; i < _localMDInstances; i++){
      _simpleMD.push_back(coupling::interface::SimulationAndInterfaceFactory::getInstance().getMDSimulation(
        _simpleMDConfig,_mamicoConfig
        #if (COUPLING_MD_PARALLEL==COUPLING_MD_YES)
        , _multiMDService->getLocalCommunicator()
        #endif
      ));
      if (_simpleMD[i]==NULL){
        std::cout << "ERROR FoamTest: _simpleMD[" << i << "]==NULL!" << std::endl;
        exit(EXIT_FAILURE);
      }
    }

    _mdStepCounter = 0;
    if (_rank == 0){ gettimeofday(&_tv.start,NULL); }

    for (unsigned int i = 0; i < _localMDInstances; i++){
      _simpleMD[i]->init(*_multiMDService,_multiMDService->getGlobalNumberOfLocalMDSimulation(i));
    }

    // equilibrate MD
    for (unsigned int i = 0; i < _localMDInstances; i++){
      _simpleMD[i]->switchOffCoupling();
      _simpleMD[i]->simulateTimesteps(_cfg.equSteps,_mdStepCounter);
    }

    _mdStepCounter += _cfg.equSteps;

    // allocate coupling interfaces
    for (unsigned int i = 0; i < _localMDInstances; i++){
      _simpleMD[i]->switchOnCoupling();
      _mdSolverInterface.push_back(coupling::interface::SimulationAndInterfaceFactory::getInstance().
        getMDSolverInterface(_simpleMDConfig, _mamicoConfig, _simpleMD[i]));
      if (_mdSolverInterface[i] == NULL){
        std::cout << "ERROR FoamTest: mdSolverInterface[" << i << "] == NULL!" << std::endl;
        exit(EXIT_FAILURE);
      }
    }

    coupling::interface::MacroscopicSolverInterface<3>* couetteSolverInterface = new coupling::solvers::CouetteSolverInterface<3>(
      getGlobalNumberMacroscopicCells(_simpleMDConfig,_mamicoConfig),_mamicoConfig.getMomentumInsertionConfiguration().getInnerOverlap());

    // initialise macroscopic cell service for multi-MD case and set single cell services in each MD simulation
    _multiMDCellService = new coupling::services::MultiMDCellService<MY_LINKEDCELL,3>(
      _mdSolverInterface,couetteSolverInterface, _simpleMDConfig.getMPIConfiguration().getNumberOfProcesses(), (unsigned int) _rank, _cfg.totalNumberMDSimulations,
      _mamicoConfig.getParticleInsertionConfiguration(), _mamicoConfig.getMomentumInsertionConfiguration(), _mamicoConfig.getBoundaryForceConfiguration(),
      _mamicoConfig.getTransferStrategyConfiguration(), _mamicoConfig.getNoiseReductionConfiguration(), _mamicoConfig.getParallelTopologyConfiguration(), _simpleMDConfig.getSimulationConfiguration().getNumberOfTimesteps(),
      _mamicoConfig.getMacroscopicCellConfiguration(), *_multiMDService
    );

      // set couette solver interface in MamicoInterfaceProvider
    coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL,3>::getInstance().setMacroscopicSolverInterface(couetteSolverInterface);

    for (unsigned int i = 0; i < _localMDInstances; i++){
      _simpleMD[i]->setMacroscopicCellService(&(_multiMDCellService->getMacroscopicCellService(i)));
      // compute and store temperature in macroscopic cells (temp=1.1 everywhere)
      _multiMDCellService->getMacroscopicCellService(i).computeAndStoreTemperature(_cfg.temp);
    }

    // allocate buffers for send/recv operations
    allocateSendBuffer(_multiMDCellService->getMacroscopicCellService(0).getIndexConversion(),*couetteSolverInterface);
    allocateRecvBuffer(_multiMDCellService->getMacroscopicCellService(0).getIndexConversion(),*couetteSolverInterface);
    _buf.foamCellIndices4SendBuffer = globalIndices2FoamIndices(_buf.globalCellIndices4SendBuffer, _buf.sendBuffer.size(), _multiMDCellService->getMacroscopicCellService(0).getIndexConversion());
    findPointsInFoamBoundary(_buf.globalCellIndices4RecvBuffer, _buf.recvBuffer.size(), _multiMDCellService->getMacroscopicCellService(0).getIndexConversion());
    // finish time measurement for initialisation
    if(_rank == 0){
      gettimeofday(&_tv.end,NULL);
      double runtime = (_tv.end.tv_sec - _tv.start.tv_sec)*1000000 + (_tv.end.tv_usec - _tv.start.tv_usec);
      std::cout << "Initialization: " << (int)(runtime/1000) << "ms" << std::endl;
    }

    if (_rank == 0){ gettimeofday(&_tv.start_total,NULL); }
    std::cout << "Finish FoamTest::initSolvers() " << std::endl;
  }

  void advanceMacro(){
    if (_rank==0){
      gettimeofday(&_tv.start,NULL);
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
      runTime.write();
      plottxt();
      gettimeofday(&_tv.end,NULL);
      _tv.macro += (_tv.end.tv_sec - _tv.start.tv_sec)*1000000 + (_tv.end.tv_usec - _tv.start.tv_usec);}
      //extract data from couette solver and send them to MD (can take any index-conversion object)
      fillSendBuffer(_cfg.density,_multiMDCellService->getMacroscopicCellService(0).getIndexConversion(),_buf.sendBuffer);
      std::cout << "Finished MacroTimestep" << std::endl;
    if(_cfg.macro2Md){
      _multiMDCellService->sendFromMacro2MD(_buf.sendBuffer,_buf.globalCellIndices4SendBuffer);
      std::cout << "Finish _multiMDCellService->sendFromMacro2MD " << std::endl;
    }
  }

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

  void advanceMicro(int cycle){
    if (_rank==0){ gettimeofday(&_tv.start,NULL); }
    std::cout << "Start MDTimpestep" << std::endl;
    // run MD instances
    for (unsigned int i = 0; i < _localMDInstances; i++){
      // set macroscopic cell service and interfaces in MamicoInterfaceProvider
      coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL,3>::getInstance().setMacroscopicCellService(&(_multiMDCellService->getMacroscopicCellService(i)));
      coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL,3>::getInstance().setMDSolverInterface(_mdSolverInterface[i]);
      _simpleMD[i]->simulateTimesteps(_simpleMDConfig.getSimulationConfiguration().getNumberOfTimesteps(),_mdStepCounter);
      //std::cout << "Finish _simpleMD[i]->simulateTimesteps " << std::endl;
      // plot macroscopic time step info in multi md service
      _multiMDCellService->getMacroscopicCellService(i).plotEveryMacroscopicTimestep(cycle);
    }
    _mdStepCounter += _simpleMDConfig.getSimulationConfiguration().getNumberOfTimesteps();
    if (_rank==0){
      gettimeofday(&_tv.end,NULL);
      _tv.micro += (_tv.end.tv_sec - _tv.start.tv_sec)*1000000 + (_tv.end.tv_usec - _tv.start.tv_usec);
    }

    // send back data from MD instances and merge it
    if(_cfg.md2Macro){
      _tv.filter += _multiMDCellService->sendFromMD2Macro(_buf.recvBuffer,_buf.globalCellIndices4RecvBuffer);
      //std::cout << "Finish _multiMDCellService->sendFromMD2Macro " << std::endl;
    }
    std::cout << "Finished MDTimestep" << std::endl;
  }

  void twoWayCoupling(int cycle){ // Function to set the values from MD within the macro solver
    if ( _cfg.twoWayCoupling & (cycle>10)){writeMDvalues2Foam(); std::cout << "Done" << std::endl;}
    //write data to csv-compatible file for evaluation
    write2CSV(_buf.recvBuffer,_buf.globalCellIndices4RecvBuffer,_multiMDCellService->getMacroscopicCellService(0).getIndexConversion(),cycle);
  }

  void writeMDvalues2Foam(){
    if(_rank==0){
      // std::stringstream ss; ss << "recvBuffer_" << runTime.timeOutputValue()/0.25 << ".txt";
      // std::ofstream file(ss.str().c_str());
      // if (!file.is_open()){std::cout << "ERROR NumericalSolver::plottxt(): Could not open file " << ss.str() << "!" << std::endl; exit(EXIT_FAILURE);}
      // std::stringstream Buffer;
      for(unsigned int i=0; i < _buf.numberFoamBoundaryPoints; i++){
        unsigned int j = _buf.foam2RecvBufferIndices[i];
        tarch::la::Vector<3,double> localVel( (1.0/_buf.recvBuffer[j]->getMacroscopicMass())*_buf.recvBuffer[j]->getMacroscopicMomentum() );
        _buf.foamBoundaryIndices[i]->x() = localVel[0];
        _buf.foamBoundaryIndices[i]->y() = localVel[1];
        _buf.foamBoundaryIndices[i]->z() = localVel[2];
        //std::cout << localVel[0] << ", " << localVel[1] << ", " << localVel[2] << ", "<< _buf.foamBoundaryIndices[i][0] << ", " << _buf.foamBoundaryIndices[i][1] << ", " << _buf.foamBoundaryIndices[i][2] << ", " << std::endl;
      }
      // for (unsigned int i = 6; i < 12; i++){
      //   for (unsigned int j = 0; j < 36; j++){
      //     Foam::vectorField FoamCoord = U.boundaryFieldRef()[i].patch().Cf()[j]+(U.boundaryFieldRef()[i].patch().nf()*1.25);
      //     Buffer << FoamCoord[0][0] << ", " << FoamCoord[0][1] << ", " << FoamCoord[0][2] << ", " << U.boundaryFieldRef()[i][j][0] << ", " << U.boundaryFieldRef()[i][j][1] << ", " << U.boundaryFieldRef()[i][j][2] << std::endl;
      //   }
      // }
      // file << Buffer.str() << std::endl;
      // file.close();
    }
  }

  void shutdown(){
    // finish time measurement for coupled simulation
    if (_rank==0){
      gettimeofday(&_tv.end,NULL);
      double time_total = (_tv.end.tv_sec - _tv.start_total.tv_sec)*1000000 + (_tv.end.tv_usec - _tv.start_total.tv_usec);
      std::cout << "Finished all coupling cycles after " << time_total/1000000 << " s" << std::endl;
      std::cout << "Time percentages Micro, Macro, Filter, Other: " << std::endl;
      std::cout << _tv.micro/time_total*100 << ", " << _tv.macro/time_total*100 << ",  "
                << _tv.filter/time_total*100 << ", " << (1-(_tv.micro+_tv.macro+_tv.filter)/time_total)*100 << std::endl;
    }

    // free buffers/arrays
    deleteBuffer(_buf.sendBuffer);
    if (_buf.globalCellIndices4SendBuffer != NULL){
      delete [] _buf.globalCellIndices4SendBuffer;
      _buf.globalCellIndices4SendBuffer = NULL;
    }
    deleteBuffer(_buf.recvBuffer);
    if (_buf.globalCellIndices4RecvBuffer != NULL){
      delete [] _buf.globalCellIndices4RecvBuffer;
      _buf.globalCellIndices4RecvBuffer = NULL;
    }

    // shutdown MD simulation
    for (unsigned int i = 0; i < _localMDInstances; i++) {
      // the shutdown operation may also delete the md solver interface; therefore, we update the MD solver interface in the vector _mdSolverInteface after the shutdown is completed
      coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL,3>::getInstance().setMDSolverInterface(_mdSolverInterface[i]);
      _simpleMD[i]->shutdown(); delete _simpleMD[i]; _simpleMD[i] = NULL;
      _mdSolverInterface[i] = coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL,3>::getInstance().getMDSolverInterface();
    }
    _simpleMD.clear();
    for (unsigned int i = 0; i < _localMDInstances; i++) {
      if (_mdSolverInterface[i] != NULL){
        delete _mdSolverInterface[i];
        _mdSolverInterface[i] = NULL;
      }
    }
    _mdSolverInterface.clear();

    coupling::interface::MacroscopicSolverInterface<3>* couetteSolverInterface =
      coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL,3>::getInstance().getMacroscopicSolverInterface();

    if( _multiMDService != NULL){delete _multiMDService; _multiMDService = NULL;}
    if (couetteSolverInterface != NULL){delete couetteSolverInterface; couetteSolverInterface = NULL;}
    if(_multiMDCellService != NULL){delete _multiMDCellService; _multiMDCellService=NULL;}

    std::cout << "Finish FoamTest::shutdown() " << std::endl;
  }

  /** computes global number of macroscopic cells from configs. Required by couette solver interface before MacroscopicCellService is initialised! */
  tarch::la::Vector<3,unsigned int> getGlobalNumberMacroscopicCells(
    const simplemd::configurations::MolecularDynamicsConfiguration &simpleMDConfig, const coupling::configurations::MaMiCoConfiguration<3> &mamicoConfig) const {
    tarch::la::Vector<3,double> domainSize(simpleMDConfig.getDomainConfiguration().getGlobalDomainSize());
    tarch::la::Vector<3,double> dx(mamicoConfig.getMacroscopicCellConfiguration().getMacroscopicCellSize());
    tarch::la::Vector<3,unsigned int> globalNumberMacroscopicCells(0);
    for (unsigned int d = 0;  d < 3; d++){
     int buf = floor(domainSize[d]/dx[d]+0.5);
     globalNumberMacroscopicCells[d] = (unsigned int) buf;
    }
    return globalNumberMacroscopicCells;
  }

  /** allocates the send buffer (with values for all macroscopic cells) and returns indices. This is only done on rank 0. */
  void allocateSendBuffer(const coupling::IndexConversion<3>& indexConversion,
  coupling::interface::MacroscopicSolverInterface<3> &couetteSolverInterface) {
    // determine global number of cells
    const tarch::la::Vector<3,unsigned int> cells(indexConversion.getGlobalNumberMacroscopicCells()+tarch::la::Vector<3,unsigned int>(2));
    const unsigned int num = cells[0]*cells[1]*cells[2];

    // delete all potential entries of sendBuffer
    deleteBuffer(_buf.sendBuffer);

    // count number of cells to be sent from this process; therefore, loop over all global macroscopic cells...
    unsigned int numCellsSent=0;
    for (unsigned int i =0; i < num; i++){
       // ... and find out, if the current cell should be send to MD from this couette solver process
       if(couetteSolverInterface.sendMacroscopicQuantityToMDSolver(indexConversion.getGlobalVectorCellIndex(i))){
         numCellsSent++;
      }
    }

    // allocate array for cell indices
    unsigned int* indices = new unsigned int [numCellsSent];
    if (indices==NULL){std::cout << "ERROR FoamTest::allocateSendBuffer(): indices==NULL!" << std::endl; exit(EXIT_FAILURE); }

    // allocate sendBuffer and initialise all entries, incl. indices
    for (unsigned int i = 0; i < num; i++){
      if (couetteSolverInterface.sendMacroscopicQuantityToMDSolver(indexConversion.getGlobalVectorCellIndex(i))){
         std::vector<unsigned int> ranks = couetteSolverInterface.getSourceRanks(indexConversion.getGlobalVectorCellIndex(i));
         bool containsThisRank=false;
         for (unsigned int k = 0; k < ranks.size(); k++){
           containsThisRank = containsThisRank || (ranks[k]==(unsigned int)_rank);
         }
         if (containsThisRank){
           _buf.sendBuffer.push_back(new coupling::datastructures::MacroscopicCell<3>());
           if (_buf.sendBuffer[_buf.sendBuffer.size()-1]==NULL){std::cout << "ERROR FoamTest::allocateSendBuffer: sendBuffer[" << _buf.sendBuffer.size()-1 << "]==NULL!" << std::endl; exit(EXIT_FAILURE);}
           indices[_buf.sendBuffer.size()-1] = i;
         }
      }
    }

    #if (COUPLING_MD_DEBUG==COUPLING_MD_YES)
    for (unsigned int i = 0; i < numCellsSent; i++){
      std::vector<unsigned int> ranks = couetteSolverInterface.getSourceRanks(indexConversion.getGlobalVectorCellIndex(indices[i]));
      std::cout << "Current rank= " << _rank << ", Send cell " << indexConversion.getGlobalVectorCellIndex(indices[i]) << ", ranks=";
      for (unsigned int j = 0; j < ranks.size(); j++){ std::cout << " " << ranks[j];}
      std::cout << std::endl;
    }
    #endif
    _buf.globalCellIndices4SendBuffer = indices;
  }

  /** allocates the recv-buffer. This buffer contains all global inner macroscopic cells, but only on rank 0. On all other ranks, no cells are stored and a NULL ptr is returned */
  void allocateRecvBuffer(const coupling::IndexConversion<3>& indexConversion,
  coupling::interface::MacroscopicSolverInterface<3> &couetteSolverInterface) {

    // determine global number of cells
    const tarch::la::Vector<3,unsigned int> cells(indexConversion.getGlobalNumberMacroscopicCells()+tarch::la::Vector<3,unsigned int>(2));
    const unsigned int num = cells[0]*cells[1]*cells[2];

    // delete all potential entries of sendBuffer
    deleteBuffer(_buf.recvBuffer);

    // determine number of cells that should be received
    unsigned int numCellsRecv = 0;
    for (unsigned int i = 0; i < num; i++){
      if(couetteSolverInterface.receiveMacroscopicQuantityFromMDSolver(indexConversion.getGlobalVectorCellIndex(i))){
        std::vector<unsigned int> ranks = couetteSolverInterface.getTargetRanks(indexConversion.getGlobalVectorCellIndex(i));
        bool containsThisRank=false;
        for (unsigned int k = 0; k < ranks.size(); k++){
          containsThisRank = containsThisRank || (ranks[k]==(unsigned int)_rank);
        }
        if (containsThisRank){ numCellsRecv++; }
      }
    }
    // allocate array for cell indices
    unsigned int* indices = new unsigned int [numCellsRecv];
    if (indices==NULL){std::cout << "ERROR FoamTest::allocateRecvBuffer(): indices==NULL!" << std::endl; exit(EXIT_FAILURE); }

    // allocate recvBuffer and initialise all entries, incl. indices
    for (unsigned int i = 0; i < num; i++){
      if (couetteSolverInterface.receiveMacroscopicQuantityFromMDSolver(indexConversion.getGlobalVectorCellIndex(i))){
        std::vector<unsigned int> ranks = couetteSolverInterface.getTargetRanks(indexConversion.getGlobalVectorCellIndex(i));
        bool containsThisRank=false;
        for (unsigned int k = 0; k < ranks.size(); k++){
          containsThisRank = containsThisRank || (ranks[k]==(unsigned int)_rank);
        }
        if (containsThisRank){
          _buf.recvBuffer.push_back(new coupling::datastructures::MacroscopicCell<3>());
          if (_buf.recvBuffer[_buf.recvBuffer.size()-1]==NULL){std::cout << "ERROR FoamTest::allocateRecvBuffer: recvBuffer[" << _buf.recvBuffer.size()-1 << "]==NULL!" << std::endl; exit(EXIT_FAILURE);}
          // set linearized index
          indices[_buf.recvBuffer.size()-1] = i;
        }
      }
    }

    #if (COUPLING_MD_DEBUG==COUPLING_MD_YES)
    for (unsigned int i = 0; i < numCellsRecv; i++){
      std::vector<unsigned int> ranks = couetteSolverInterface.getTargetRanks(indexConversion.getGlobalVectorCellIndex(indices[i]));
      std::cout << "Current rank= " << _rank << ", Recv cell " << indexConversion.getGlobalVectorCellIndex(indices[i]) << ", ranks=";
      for (unsigned int j = 0; j < ranks.size(); j++){ std::cout << " " << ranks[j];}
      std::cout << std::endl;
    }
    #endif
    _buf.globalCellIndices4RecvBuffer = indices;
  }

  unsigned int* globalIndices2FoamIndices(const unsigned int* const Index, unsigned int size,
    const coupling::IndexConversion<3>& indexConversion){
    unsigned int* foamIndice = new unsigned int [size];
    const tarch::la::Vector<3,double> domainOffset(indexConversion.getGlobalMDDomainOffset());
    const tarch::la::Vector<3,double> cellSize(2.5);

    for(unsigned int i=0; i<size;i++){
      const tarch::la::Vector<3,unsigned int> globalIndex(indexConversion.getGlobalVectorCellIndex(Index[i]));
      tarch::la::Vector<3,double> cellMidPoint(domainOffset-0.5*cellSize);
      for (unsigned int d = 0; d < 3; d++){ cellMidPoint[d] = cellMidPoint[d] + ((double)globalIndex[d])*2.5; }
      const Foam::vector foamPosition(cellMidPoint[0],cellMidPoint[1],cellMidPoint[2]);
      foamIndice[i] = U.mesh().findCell(foamPosition) > 0 ? mesh.findCell(foamPosition) : 0;
    }
    return foamIndice;
  }

  void findPointsInFoamBoundary(const unsigned int* const globalIndice, unsigned int size,
    const coupling::IndexConversion<3>& indexConversion){
    _buf.numberFoamBoundaryPoints = 6*36;
    _buf.foam2RecvBufferIndices = new unsigned int [_buf.numberFoamBoundaryPoints];
    _buf.foamBoundaryIndices = new Foam::vector* [_buf.numberFoamBoundaryPoints];
    unsigned int counter = 0;

    // const tarch::la::Vector<3,double> domainOffset(indexConversion.getGlobalMDDomainOffset());
    // const tarch::la::Vector<3,double> cellSize(2.5);

    // std::stringstream ss; ss << "FoamCoodinates.txt";
    // std::ofstream file(ss.str().c_str());
    // if (!file.is_open()){std::cout << "ERROR NumericalSolver::plottxt(): Could not open file " << ss.str() << "!" << std::endl; exit(EXIT_FAILURE);}
    // std::stringstream coord;
    for (unsigned int i = 6; i < 12; i++){
      for (unsigned int j = 0; j < 36; j++){
        Foam::vectorField FoamCoord = U.boundaryFieldRef()[i].patch().Cf()[j]+(U.boundaryFieldRef()[i].patch().nf()*1.25);
        _buf.foamBoundaryIndices[counter] = &(U.boundaryFieldRef()[i][j]);
        const unsigned int globalIndex = indexConversion.getGlobalCellIndex(indexConversion.getGlobalVectorCellIndex(tarch::la::Vector<3,double>(FoamCoord[0][0],FoamCoord[0][1],FoamCoord[0][2])));
        for(unsigned int k = 0; k < size; k++){
          if(globalIndex==globalIndice[k]){
            _buf.foam2RecvBufferIndices[counter] = k;
            // const tarch::la::Vector<3,unsigned int> globalIndex1(indexConversion.getGlobalVectorCellIndex(globalIndex));
            // tarch::la::Vector<3,double> cellMidPoint(domainOffset-0.5*cellSize);
            // for (unsigned int d = 0; d < 3; d++){ cellMidPoint[d] = cellMidPoint[d] + ((double)globalIndex1[d])*2.5; }
            // coord << FoamCoord[0][0] << ", " << FoamCoord[0][1] << ", " << FoamCoord[0][2] << ", ";
            // coord << cellMidPoint[0] << ", " << cellMidPoint[1] << ", " << cellMidPoint[2] << ", " ;
            // coord << globalIndex << ", " << i << ", " << j << ", " << _buf.foamNumberIndices4RecvBuffer[k] << std::endl;
            goto endloop;
          }
        }
      endloop:
      counter++;
      }
    }
    // file << coord.str();
    // file.close();
  }

  /** write cells that have been received from MD to csv file */
  void write2CSV(
    std::vector<coupling::datastructures::MacroscopicCell<3>* >& recvBuffer,const unsigned int * const recvIndices,
    const coupling::IndexConversion<3>& indexConversion, int couplingCycle
  ) const {
    if(_cfg.csvEveryTimestep < 1 || couplingCycle % _cfg.csvEveryTimestep > 0) return;
    // form file name and open file
    std::stringstream ss;
    ss << "CouetteAvgMultiMDCells_" << _rank << "_" << couplingCycle << ".csv";
    std::ofstream file(ss.str().c_str());
    if (!file.is_open()){std::cout << "ERROR FoamTest::write2CSV(): Could not open file " << ss.str() << "!" << std::endl; exit(EXIT_FAILURE);}

    // loop over received cells; read macroscopic mass+momentum buffers and write cell index, mass and velocity to one line in the csv-file
    const unsigned int numCellsRecv = recvBuffer.size();
    for (unsigned int i = 0; i < numCellsRecv; i++){
      tarch::la::Vector<3,double> vel(recvBuffer[i]->getMacroscopicMomentum());
      if (recvBuffer[i]->getMacroscopicMass()!=0.0){ vel = (1.0/recvBuffer[i]->getMacroscopicMass())*vel; }
      const tarch::la::Vector<3,unsigned int> counter(indexConversion.getGlobalVectorCellIndex(recvIndices[i]));
      file << counter[0] << " ; " << counter[1] << " ; " << counter[2] << " ; " << vel[0] << " ; " << vel[1] << " ; " << vel[2] << " ; " << recvBuffer[i]->getMacroscopicMass() << ";";
      file << std::endl;
    }

    // close file
    file.close();
  }


  /** deletes the send buffer */
  void deleteBuffer(std::vector<coupling::datastructures::MacroscopicCell<3>* >& sendBuffer) const {
    // delete all potential entries of sendBuffer
    for (unsigned int i = 0; i < sendBuffer.size(); i++){if (sendBuffer[i]!=NULL){ delete sendBuffer[i]; sendBuffer[i]=NULL;}}
    sendBuffer.clear();
  }

  /** fills send buffer with data from couette solver */
  void fillSendBuffer(
    const double density, const coupling::IndexConversion<3>& indexConversion,
    std::vector<coupling::datastructures::MacroscopicCell<3>* >& sendBuffer
  ) const {
    const unsigned int size = sendBuffer.size();
    const tarch::la::Vector<3,double> macroscopicCellSize(indexConversion.getMacroscopicCellSize());
    double mass = density * macroscopicCellSize[0]*macroscopicCellSize[1]*macroscopicCellSize[2];

    // std::stringstream ss; ss << "buffer_" << runTime.timeName() << ".txt";
    // std::ofstream file(ss.str().c_str());
    // if (!file.is_open()){std::cout << "ERROR NumericalSolver::plottxt(): Could not open file " << ss.str() << "!" << std::endl; exit(EXIT_FAILURE);}
    // std::stringstream output;

    for (unsigned int i = 0; i < size; i++){
      if(_buf.foamCellIndices4SendBuffer[i]>0){
        const tarch::la::Vector<3,double> velocity(U[_buf.foamCellIndices4SendBuffer[i]][0], U[_buf.foamCellIndices4SendBuffer[i]][1], U[_buf.foamCellIndices4SendBuffer[i]][2]);
        //const tarch::la::Vector<3,double> momentum(mass*U[_buf.foamCellIndices4SendBuffer[i]]);
        sendBuffer[i]->setMicroscopicMass(mass);
        sendBuffer[i]->setMicroscopicMomentum(mass*velocity);
        // write information to streams
        // output << velocity[0] << ", " << velocity[1] << ", " << velocity[2] << ", " << mass*velocity[0] << ", " << mass*velocity[1] << ", " << mass*velocity[2] << std::endl;
      }
    }
    // file << output.str() << std::endl;
    // output.str("");
    // file.close();
  }

    struct CouetteConfig{
    // number of coupling cycles, that is continuum time steps; MD/DPD: 1000
    int couplingCycles;
    bool md2Macro, macro2Md, twoWayCoupling;
    int csvEveryTimestep;
    double density;
    // only for LB couette solver: VTK plotting per time step
    int plotEveryTimestep;
    double temp;
    int equSteps;
    int totalNumberMDSimulations;
  };

  struct CouplingBuffer{
    std::vector<coupling::datastructures::MacroscopicCell<3>* > sendBuffer;
    unsigned int *globalCellIndices4SendBuffer;
    unsigned int *foamCellIndices4SendBuffer;
    std::vector<coupling::datastructures::MacroscopicCell<3>* > recvBuffer;
    unsigned int *globalCellIndices4RecvBuffer;
    Foam::vector** foamBoundaryIndices;
    unsigned int *foam2RecvBufferIndices;
    unsigned int numberFoamBoundaryPoints;
  };

  struct TimingValues{
    timeval start_total;
    timeval start;
    timeval end;
    double micro;
    double macro;
    double filter;
  };

  int _rank;
  simplemd::configurations::MolecularDynamicsConfiguration _simpleMDConfig;
  coupling::configurations::MaMiCoConfiguration<3> _mamicoConfig;
  CouetteConfig _cfg;
  unsigned int _mdStepCounter;
  tarch::utils::MultiMDService<3>* _multiMDService;
  coupling::services::MultiMDCellService<MY_LINKEDCELL,3> *_multiMDCellService;
  CouplingBuffer _buf;
  unsigned int _localMDInstances;
  std::vector<coupling::interface::MDSolverInterface<MY_LINKEDCELL,3>* > _mdSolverInterface;
  std::vector<coupling::interface::MDSimulation*> _simpleMD;
  TimingValues _tv;
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
#endif // _COUPLING_TESTS_FOAMTEST_H_