// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#include "simplemd/MolecularDynamicsSimulation.h"
#include "simplemd/cell-mappings/VaryCheckpointMapping.h"

simplemd::MolecularDynamicsSimulation::MolecularDynamicsSimulation(const simplemd::configurations::MolecularDynamicsConfiguration& configuration)
    : _configuration(configuration), _timeIntegrator(NULL), _updateLinkedCellListsMapping(NULL), _vtkMoleculeWriter(NULL),
#if BUILD_WITH_ADIOS2
      _Adios2Writer(NULL),
#endif
      _lennardJonesForce(NULL), _emptyLinkedListsMapping(NULL), _rdfMapping(NULL), _boundaryTreatment(NULL), _localMDSimulation(0), _profilePlotter(NULL),
      _parallelTopologyService(NULL), _moleculeService(NULL), _linkedCellService(NULL), _molecularPropertiesService(NULL),
      // initialise external forces
      _externalForceService(configuration.getExternalForceConfigurations()) {
}

double simplemd::MolecularDynamicsSimulation::getNumberDensity(unsigned int numberMolecules, const tarch::la::Vector<MD_DIM, double>& domainSize) const {
  double density = 1.0;
  for (unsigned int d = 0; d < MD_DIM; d++) {
    density = density / domainSize[d];
  }
  return (density * numberMolecules);
}

// TODO: fix duplicate copy-paste code in both versions of initServices
void simplemd::MolecularDynamicsSimulation::initServices() {
  // set vtk file stem and checkpoint filestem -> only one MD simulation runs
  _localMDSimulation = 0;
  _vtkFilestem = _configuration.getVTKConfiguration().getFilename();
  _checkpointFilestem = _configuration.getCheckpointConfiguration().getFilename();
  // initialise local variable with global information first (they are adapted
  // later on)
  tarch::la::Vector<MD_DIM, double> localDomainSize(_configuration.getDomainConfiguration().getGlobalDomainSize());
  tarch::la::Vector<MD_DIM, double> localDomainOffset(_configuration.getDomainConfiguration().getGlobalDomainOffset());
  tarch::la::Vector<MD_DIM, unsigned int> moleculesPerDirection(_configuration.getDomainConfiguration().getMoleculesPerDirection());
  tarch::la::Vector<MD_DIM, double> realMeshWidth(0.0);
  tarch::la::Vector<MD_DIM, unsigned int> processCoordinates(0);
  double linkedCellVolume = 1.0;
  simplemd::moleculemappings::InitialPositionAndForceUpdate initialPositionAndForceUpdate(_configuration.getSimulationConfiguration().getDt(),
                                                                                          _configuration.getMoleculeConfiguration().getMass());
  // initialise services -> initialise ParallelTopologyService first (also for
  // serial case!) and then, initialise other services, if this process
  // contributes to the simulation in some way Note: initBuffers must also be
  // called on non-idle processors, but only after MoleculeService has been
  // initialised.
  _parallelTopologyService = new simplemd::services::ParallelTopologyService(
      _configuration.getDomainConfiguration().getGlobalDomainSize(), _configuration.getDomainConfiguration().getGlobalDomainOffset(),
      _configuration.getDomainConfiguration().getMeshWidth(), _configuration.getMPIConfiguration().getNumberOfProcesses(),
      _configuration.getDomainConfiguration().getBoundary()
#if (MD_PARALLEL == MD_YES)
          ,
      MPI_COMM_WORLD
#endif
  );
  if (_parallelTopologyService == NULL) {
    std::cout << "ERROR simplemd::MolecularDynamicsSimulation::initServices(): "
                 "_parallelTopologyService==NULL!"
              << std::endl;
    exit(EXIT_FAILURE);
  }

  if (!_parallelTopologyService->isIdle()) {

    // set local periodic boundary information in periodicBoundary from here on
    _localBoundary = _parallelTopologyService->getLocalBoundaryInformation();
    // get process coordinates for the local process
    processCoordinates = _parallelTopologyService->getProcessCoordinates();
    // get local meshwidth defined by ParallelTopologyService
    realMeshWidth = _parallelTopologyService->getMeshWidth();
    // determine local domain size and local domain offset and local number of
    // molecules
    for (unsigned int d = 0; d < MD_DIM; d++) {
      localDomainSize[d] = localDomainSize[d] / (_configuration.getMPIConfiguration().getNumberOfProcesses()[d]);
      localDomainOffset[d] = localDomainOffset[d] + (processCoordinates[d] * _parallelTopologyService->getLocalNumberOfCells()[d] * realMeshWidth[d]);
      int missingMolecules = moleculesPerDirection[d] % (_configuration.getMPIConfiguration().getNumberOfProcesses()[d]);
      moleculesPerDirection[d] = moleculesPerDirection[d] / (_configuration.getMPIConfiguration().getNumberOfProcesses()[d]);
      for (int i = 0; i < missingMolecules; i++) {
        if (processCoordinates[d] == (unsigned int)(i * (1.0 * _configuration.getMPIConfiguration().getNumberOfProcesses()[d]) / missingMolecules))
          moleculesPerDirection[d]++;
      }
    }
#if (MD_DEBUG == MD_YES)
    if (_parallelTopologyService->getProcessCoordinates() == tarch::la::Vector<MD_DIM, unsigned int>(0)) {
      std::cout << "Local number of cells: " << _parallelTopologyService->getLocalNumberOfCells() << std::endl;
    }
    std::cout << "Rank: " << _parallelTopologyService->getRank() << " Boundary: " << _localBoundary << ", process coordinates: " << processCoordinates
              << std::endl;
    std::cout << " Local Domain size: " << localDomainSize << ", local domain offset: " << localDomainOffset << std::endl;
#endif

    tarch::utils::RandomNumberService::getInstance().init(_configuration.getSimulationConfiguration().fixSeed());
    _molecularPropertiesService = new simplemd::services::MolecularPropertiesService(
        _configuration.getMoleculeConfiguration().getMass(), _configuration.getMoleculeConfiguration().getEpsilon(),
        _configuration.getMoleculeConfiguration().getSigma(), _configuration.getDomainConfiguration().getCutoffRadius(),
        _configuration.getDomainConfiguration().getKB());
    if (_molecularPropertiesService == NULL) {
      std::cout << "ERROR MolecularDynamicsSimulation::initServices(): "
                   "_molecularPropertiesService==NULL!"
                << std::endl;
      exit(EXIT_FAILURE);
    }

    // either initialise from checkpoint data or via a certain number of
    // molecules per direction
    if (_configuration.getDomainConfiguration().initFromCheckpoint()) {
      _moleculeService =
          new simplemd::services::MoleculeService(localDomainSize, localDomainOffset, _configuration.getDomainConfiguration().getCheckpointFilestem(),
                                                  _configuration.getDomainConfiguration().getBlockSize(), *_parallelTopologyService);
    } else if (_configuration.getDomainConfiguration().initFromSequentialCheckpoint()) {
      _moleculeService =
          new simplemd::services::MoleculeService(localDomainSize, localDomainOffset, _configuration.getDomainConfiguration().getCheckpointFilestem(),
                                                  _configuration.getDomainConfiguration().getBlockSize());
    } else {
      _moleculeService = new simplemd::services::MoleculeService(
          localDomainSize, localDomainOffset, moleculesPerDirection, _configuration.getMoleculeConfiguration().getMeanVelocity(),
          _configuration.getDomainConfiguration().getKB(), _configuration.getMoleculeConfiguration().getTemperature(),
          _configuration.getDomainConfiguration().getBlockSize(), *_molecularPropertiesService);
    }
    if (_moleculeService == NULL) {
      std::cout << "ERROR MolecularDynamicsSimulation::initServices(): "
                   "_moleculeService==NULL!"
                << std::endl;
      exit(EXIT_FAILURE);
    }

    // initialise buffers. After this call, the ParallelTopologyService
    // initialisation is complete
    _parallelTopologyService->initBuffers(_moleculeService->getNumberMolecules());
    _linkedCellService = new simplemd::services::LinkedCellService(localDomainSize, localDomainOffset, *_parallelTopologyService, *_moleculeService);
    if (_linkedCellService == NULL) {
      std::cout << "ERROR simplemd::MolecularDynamicsSimulation::initServices(): "
                   "_linkedCellService==NULL!"
                << std::endl;
      exit(EXIT_FAILURE);
    }

    // init boundary treatment
    _boundaryTreatment = new simplemd::BoundaryTreatment(*_parallelTopologyService, *_moleculeService, *_linkedCellService);
    if (_boundaryTreatment == NULL) {
      std::cout << "ERROR simplemd::MolecularDynamicsSimulation::initServices(): "
                   "_boundaryTreatment==NULL!"
                << std::endl;
      exit(EXIT_FAILURE);
    }

    // init all mappings here
    _timeIntegrator = new simplemd::moleculemappings::VelocityStoermerVerletMapping(
        _configuration.getDomainConfiguration().getKB(), _configuration.getSimulationConfiguration().getDt(),
        _configuration.getMoleculeConfiguration().getMass(), _configuration.getDomainConfiguration().getBoundary(),
        _configuration.getDomainConfiguration().getGlobalDomainOffset(), _configuration.getDomainConfiguration().getGlobalDomainSize());
    if (_timeIntegrator == NULL) {
      std::cout << "ERROR simplemd::MolecularDynamicsSimulation::initServices(): "
                   "_timeIntegrator==NULL!"
                << std::endl;
      exit(EXIT_FAILURE);
    }
    _updateLinkedCellListsMapping = new simplemd::moleculemappings::UpdateLinkedCellListsMapping(*_parallelTopologyService, *_linkedCellService);
    if (_updateLinkedCellListsMapping == NULL) {
      std::cout << "ERROR simplemd::MolecularDynamicsSimulation::initServices(): "
                   "_updateLinkedCellListsMapping==NULL!"
                << std::endl;
      exit(EXIT_FAILURE);
    }
    _vtkMoleculeWriter = new simplemd::moleculemappings::VTKMoleculeWriter(*_parallelTopologyService, *_moleculeService, _vtkFilestem);
    if (_vtkMoleculeWriter == NULL) {
      std::cout << "ERROR simplemd::MolecularDynamicsSimulation::initServices(): "
                   "_vtkMoleculeWriter==NULL!"
                << std::endl;
      exit(EXIT_FAILURE);
    }

#if BUILD_WITH_ADIOS2
    _Adios2Writer = new simplemd::moleculemappings::Adios2Writer(*_parallelTopologyService, *_moleculeService, _configuration
#if (MD_PARALLEL == MD_YES)
                                                                 ,
                                                                 MPI_COMM_WORLD
#endif
    );
    if (_Adios2Writer == NULL) {
      std::cout << "ERROR simplemd::MolecularDynamicsSimulation::initServices(): "
                   "_Adios2Writer==NULL!"
                << std::endl;
      exit(EXIT_FAILURE);
    }
#endif

    // cell mappings
    _lennardJonesForce = new simplemd::cellmappings::LennardJonesForceMapping(_externalForceService, *_molecularPropertiesService);
    if (_lennardJonesForce == NULL) {
      std::cout << "ERROR simplemd::MolecularDynamicsSimulation::initServices(): "
                   "_lennardJonesForce==NULL!"
                << std::endl;
      exit(EXIT_FAILURE);
    }
    _emptyLinkedListsMapping = new simplemd::cellmappings::EmptyLinkedListsMapping();
    if (_emptyLinkedListsMapping == NULL) {
      std::cout << "ERROR simplemd::MolecularDynamicsSimulation::initServices(): "
                   "_emptyLinkedListsMapping==NULL!"
                << std::endl;
      exit(EXIT_FAILURE);
    }
    _rdfMapping =
        new simplemd::cellmappings::RDFMapping(*_parallelTopologyService, *_linkedCellService, _configuration.getDomainConfiguration().getCutoffRadius(),
                                               _configuration.getRDFConfiguration().getNumberOfPoints());
    if (_rdfMapping == NULL) {
      std::cout << "ERROR simplemd::MolecularDynamicsSimulation::initServices(): "
                   "_rdfMapping==NULL!"
                << std::endl;
      exit(EXIT_FAILURE);
    }

    // initialise profile plotter
    if (_profilePlotter != NULL) {
      delete _profilePlotter;
      _profilePlotter = NULL;
    }
    linkedCellVolume = 1.0;
    for (unsigned int d = 0; d < MD_DIM; d++) {
      linkedCellVolume = linkedCellVolume * _linkedCellService->getMeshWidth()[d];
    }
    _profilePlotter = new simplemd::ProfilePlotter(_configuration.getProfilePlotterConfigurations(), *_parallelTopologyService, *_linkedCellService,
                                                   linkedCellVolume, _localMDSimulation);
    if (_profilePlotter == NULL) {
      std::cout << "ERROR simplemd::MolecularDynamicsSimulation::initService(): "
                   "_profilePlotter==NULL!"
                << std::endl;
      exit(EXIT_FAILURE);
    }

    // -------------- do initial force computations and position update
    // ----------
    _boundaryTreatment->fillBoundaryCells(_localBoundary, *_parallelTopologyService);
    // compute forces between molecules.
    // After this step, each molecule has received all force contributions from
    // its neighbors.
    _linkedCellService->iterateCellPairs(*_lennardJonesForce, false);
    _boundaryTreatment->emptyGhostBoundaryCells();
    _linkedCellService->iterateCells(*_emptyLinkedListsMapping, false);
    _moleculeService->iterateMolecules(initialPositionAndForceUpdate, false);

    // sort molecules into linked cells
    _moleculeService->iterateMolecules(*_updateLinkedCellListsMapping, false);

    // -------------- do initial force computations and position update (end)
    // ----------
  } // end is process not idle
}

void simplemd::MolecularDynamicsSimulation::initServices(const tarch::utils::MultiMDService<MD_DIM>& multiMDService, unsigned int localMDSimulation) {
  // set vtk file stem and checkpoint filestem and adios2 filestem;
  _localMDSimulation = localMDSimulation;
  std::stringstream filestems;
  filestems << _configuration.getVTKConfiguration().getFilename() << "_" << localMDSimulation << "_";
  _vtkFilestem = filestems.str();
  filestems.str("");
  filestems << _configuration.getCheckpointConfiguration().getFilename() << "_" << localMDSimulation << "_";
  _checkpointFilestem = filestems.str();
  filestems.str("");
  // initialise local variable with global information first (they are adapted
  // later on)
  tarch::la::Vector<MD_DIM, double> localDomainSize(_configuration.getDomainConfiguration().getGlobalDomainSize());
  tarch::la::Vector<MD_DIM, double> localDomainOffset(_configuration.getDomainConfiguration().getGlobalDomainOffset());
  tarch::la::Vector<MD_DIM, unsigned int> moleculesPerDirection(_configuration.getDomainConfiguration().getMoleculesPerDirection());
  tarch::la::Vector<MD_DIM, double> realMeshWidth(0.0);
  tarch::la::Vector<MD_DIM, unsigned int> processCoordinates(0);
  double linkedCellVolume = 1.0;
  simplemd::moleculemappings::InitialPositionAndForceUpdate initialPositionAndForceUpdate(_configuration.getSimulationConfiguration().getDt(),
                                                                                          _configuration.getMoleculeConfiguration().getMass());
  // initialise services -> initialise ParallelTopologyService first (also for
  // serial case!) and then, initialise other services, if this process
  // contributes to the simulation in some way Note: initBuffers must also be
  // called on non-idle processors, but only after MoleculeService has been
  // initialised.
  _parallelTopologyService = new simplemd::services::ParallelTopologyService(
      _configuration.getDomainConfiguration().getGlobalDomainSize(), _configuration.getDomainConfiguration().getGlobalDomainOffset(),
      _configuration.getDomainConfiguration().getMeshWidth(), _configuration.getMPIConfiguration().getNumberOfProcesses(),
      _configuration.getDomainConfiguration().getBoundary()
// DIFFERENCE TO initServices(): use communicator of MultiMDService
#if (MD_PARALLEL == MD_YES)
          ,
      multiMDService.getLocalCommunicator()
#endif
  );
  if (_parallelTopologyService == NULL) {
    std::cout << "ERROR simplemd::MolecularDynamicsSimulation::initServices(): "
                 "_parallelTopologyService==NULL!"
              << std::endl;
    exit(EXIT_FAILURE);
  }

  if (!_parallelTopologyService->isIdle()) {

    // set local periodic boundary information in periodicBoundary from here on
    _localBoundary = _parallelTopologyService->getLocalBoundaryInformation();
    // get process coordinates for the local process
    processCoordinates = _parallelTopologyService->getProcessCoordinates();
    // get local meshwidth defined by ParallelTopologyService
    realMeshWidth = _parallelTopologyService->getMeshWidth();
    // determine local domain size and local domain offset and local number of
    // molecules
    for (unsigned int d = 0; d < MD_DIM; d++) {
      localDomainSize[d] = localDomainSize[d] / (_configuration.getMPIConfiguration().getNumberOfProcesses()[d]);
      localDomainOffset[d] = localDomainOffset[d] + (processCoordinates[d] * _parallelTopologyService->getLocalNumberOfCells()[d] * realMeshWidth[d]);
      int missingMolecules = moleculesPerDirection[d] % (_configuration.getMPIConfiguration().getNumberOfProcesses()[d]);
      moleculesPerDirection[d] = moleculesPerDirection[d] / (_configuration.getMPIConfiguration().getNumberOfProcesses()[d]);
      for (int i = 0; i < missingMolecules; i++) {
        if (processCoordinates[d] == (unsigned int)(i * (1.0 * _configuration.getMPIConfiguration().getNumberOfProcesses()[d]) / missingMolecules))
          moleculesPerDirection[d]++;
      }
    }
#if (MD_DEBUG == MD_YES)
    if (_parallelTopologyService->getProcessCoordinates() == tarch::la::Vector<MD_DIM, unsigned int>(0)) {
      std::cout << "Local number of cells: " << _parallelTopologyService->getLocalNumberOfCells() << std::endl;
    }
    std::cout << "Rank: " << _parallelTopologyService->getRank() << " Boundary: " << _localBoundary << ", process coordinates: " << processCoordinates
              << std::endl;
    std::cout << " Local Domain size: " << localDomainSize << ", local domain offset: " << localDomainOffset << std::endl;
#endif

    tarch::utils::RandomNumberService::getInstance().init(_configuration.getSimulationConfiguration().fixSeed());
    _molecularPropertiesService = new simplemd::services::MolecularPropertiesService(
        _configuration.getMoleculeConfiguration().getMass(), _configuration.getMoleculeConfiguration().getEpsilon(),
        _configuration.getMoleculeConfiguration().getSigma(), _configuration.getDomainConfiguration().getCutoffRadius(),
        _configuration.getDomainConfiguration().getKB());
    if (_molecularPropertiesService == NULL) {
      std::cout << "ERROR MolecularDynamicsSimulation::initServices(): "
                   "_molecularPropertiesService==NULL!"
                << std::endl;
      exit(EXIT_FAILURE);
    }

    // either initialise from checkpoint data or via a certain number of
    // molecules per direction
    if (_configuration.getDomainConfiguration().initFromCheckpoint()) {
      _moleculeService =
          new simplemd::services::MoleculeService(localDomainSize, localDomainOffset, _configuration.getDomainConfiguration().getCheckpointFilestem(),
                                                  _configuration.getDomainConfiguration().getBlockSize(), *_parallelTopologyService);
    } else if (_configuration.getDomainConfiguration().initFromSequentialCheckpoint()) {
      _moleculeService =
          new simplemd::services::MoleculeService(localDomainSize, localDomainOffset, _configuration.getDomainConfiguration().getCheckpointFilestem(),
                                                  _configuration.getDomainConfiguration().getBlockSize());
    } else {
      _moleculeService = new simplemd::services::MoleculeService(
          localDomainSize, localDomainOffset, moleculesPerDirection, _configuration.getMoleculeConfiguration().getMeanVelocity(),
          _configuration.getDomainConfiguration().getKB(), _configuration.getMoleculeConfiguration().getTemperature(),
          _configuration.getDomainConfiguration().getBlockSize(), *_molecularPropertiesService);
    }
    if (_moleculeService == NULL) {
      std::cout << "ERROR MolecularDynamicsSimulation::initServices(): "
                   "_moleculeService==NULL!"
                << std::endl;
      exit(EXIT_FAILURE);
    }

    // initialise buffers. After this call, the ParallelTopologyService
    // initialisation is complete
    _parallelTopologyService->initBuffers(_moleculeService->getNumberMolecules());
    _linkedCellService = new simplemd::services::LinkedCellService(localDomainSize, localDomainOffset, *_parallelTopologyService, *_moleculeService);
    if (_linkedCellService == NULL) {
      std::cout << "ERROR simplemd::MolecularDynamicsSimulation::initServices(): "
                   "_linkedCellService==NULL!"
                << std::endl;
      exit(EXIT_FAILURE);
    }

    // init boundary treatment
    _boundaryTreatment = new simplemd::BoundaryTreatment(*_parallelTopologyService, *_moleculeService, *_linkedCellService);
    if (_boundaryTreatment == NULL) {
      std::cout << "ERROR simplemd::MolecularDynamicsSimulation::initServices(): "
                   "_boundaryTreatment==NULL!"
                << std::endl;
      exit(EXIT_FAILURE);
    }

    // init all mappings here
    _timeIntegrator = new simplemd::moleculemappings::VelocityStoermerVerletMapping(
        _configuration.getDomainConfiguration().getKB(), _configuration.getSimulationConfiguration().getDt(),
        _configuration.getMoleculeConfiguration().getMass(), _configuration.getDomainConfiguration().getBoundary(),
        _configuration.getDomainConfiguration().getGlobalDomainOffset(), _configuration.getDomainConfiguration().getGlobalDomainSize());
    if (_timeIntegrator == NULL) {
      std::cout << "ERROR simplemd::MolecularDynamicsSimulation::initServices(): "
                   "_timeIntegrator==NULL!"
                << std::endl;
      exit(EXIT_FAILURE);
    }
    _updateLinkedCellListsMapping = new simplemd::moleculemappings::UpdateLinkedCellListsMapping(*_parallelTopologyService, *_linkedCellService);
    if (_updateLinkedCellListsMapping == NULL) {
      std::cout << "ERROR simplemd::MolecularDynamicsSimulation::initServices(): "
                   "_updateLinkedCellListsMapping==NULL!"
                << std::endl;
      exit(EXIT_FAILURE);
    }
    _vtkMoleculeWriter = new simplemd::moleculemappings::VTKMoleculeWriter(*_parallelTopologyService, *_moleculeService, _vtkFilestem);
    if (_vtkMoleculeWriter == NULL) {
      std::cout << "ERROR simplemd::MolecularDynamicsSimulation::initServices(): "
                   "_vtkMoleculeWriter==NULL!"
                << std::endl;
      exit(EXIT_FAILURE);
    }

#if BUILD_WITH_ADIOS2
    _Adios2Writer = new simplemd::moleculemappings::Adios2Writer(*_parallelTopologyService, *_moleculeService, _configuration
#if (MD_PARALLEL == MD_YES)
                                                                 ,
                                                                 multiMDService.getLocalCommunicator()
#endif
    );
    if (_Adios2Writer == NULL) {
      std::cout << "ERROR simplemd::MolecularDynamicsSimulation::initServices(): "
                   "_Adios2Writer==NULL!"
                << std::endl;
      exit(EXIT_FAILURE);
    }
#endif

    // cell mappings
    _lennardJonesForce = new simplemd::cellmappings::LennardJonesForceMapping(_externalForceService, *_molecularPropertiesService);
    if (_lennardJonesForce == NULL) {
      std::cout << "ERROR simplemd::MolecularDynamicsSimulation::initServices(): "
                   "_lennardJonesForce==NULL!"
                << std::endl;
      exit(EXIT_FAILURE);
    }
    _emptyLinkedListsMapping = new simplemd::cellmappings::EmptyLinkedListsMapping();
    if (_emptyLinkedListsMapping == NULL) {
      std::cout << "ERROR simplemd::MolecularDynamicsSimulation::initServices(): "
                   "_emptyLinkedListsMapping==NULL!"
                << std::endl;
      exit(EXIT_FAILURE);
    }
    _rdfMapping =
        new simplemd::cellmappings::RDFMapping(*_parallelTopologyService, *_linkedCellService, _configuration.getDomainConfiguration().getCutoffRadius(),
                                               _configuration.getRDFConfiguration().getNumberOfPoints());
    if (_rdfMapping == NULL) {
      std::cout << "ERROR simplemd::MolecularDynamicsSimulation::initServices(): "
                   "_rdfMapping==NULL!"
                << std::endl;
      exit(EXIT_FAILURE);
    }

    // initialise profile plotter
    if (_profilePlotter != NULL) {
      delete _profilePlotter;
      _profilePlotter = NULL;
    }
    linkedCellVolume = 1.0;
    for (unsigned int d = 0; d < MD_DIM; d++) {
      linkedCellVolume = linkedCellVolume * _linkedCellService->getMeshWidth()[d];
    }
    _profilePlotter = new simplemd::ProfilePlotter(_configuration.getProfilePlotterConfigurations(), *_parallelTopologyService, *_linkedCellService,
                                                   linkedCellVolume, _localMDSimulation);
    if (_profilePlotter == NULL) {
      std::cout << "ERROR simplemd::MolecularDynamicsSimulation::initService(): "
                   "_profilePlotter==NULL!"
                << std::endl;
      exit(EXIT_FAILURE);
    }

    // -------------- do initial force computations and position update
    // ----------
    _boundaryTreatment->fillBoundaryCells(_localBoundary, *_parallelTopologyService);
    // compute forces between molecules.
    // After this step, each molecule has received all force contributions from
    // its neighbors.
    _linkedCellService->iterateCellPairs(*_lennardJonesForce, false);
    _boundaryTreatment->emptyGhostBoundaryCells();
    _linkedCellService->iterateCells(*_emptyLinkedListsMapping, false);
    _moleculeService->iterateMolecules(initialPositionAndForceUpdate, false);
    // sort molecules into linked cells
    _moleculeService->iterateMolecules(*_updateLinkedCellListsMapping, false);
    // -------------- do initial force computations and position update (end)
    // ----------
  } // end is process not idle
}

void simplemd::MolecularDynamicsSimulation::shutdownServices() {
  // shutdown services
  if (!_parallelTopologyService->isIdle()) {
    _linkedCellService->shutdown();
    _moleculeService->shutdown();
    _molecularPropertiesService->shutdown();
    tarch::utils::RandomNumberService::getInstance().shutdown();
  }
  _parallelTopologyService->shutdown();

  // delete mappings
  if (_timeIntegrator != NULL) {
    delete _timeIntegrator;
    _timeIntegrator = NULL;
  }
  if (_updateLinkedCellListsMapping != NULL) {
    delete _updateLinkedCellListsMapping;
    _updateLinkedCellListsMapping = NULL;
  }
  if (_vtkMoleculeWriter != NULL) {
    delete _vtkMoleculeWriter;
    _vtkMoleculeWriter = NULL;
  }
#if BUILD_WITH_ADIOS2
  if (_Adios2Writer != NULL) {
    delete _Adios2Writer;
    _Adios2Writer = NULL;
  }
#endif
  if (_lennardJonesForce != NULL) {
    delete _lennardJonesForce;
    _lennardJonesForce = NULL;
  }
  if (_emptyLinkedListsMapping != NULL) {
    delete _emptyLinkedListsMapping;
    _emptyLinkedListsMapping = NULL;
  }
  if (_rdfMapping != NULL) {
    delete _rdfMapping;
    _rdfMapping = NULL;
  }

  // delete profile plotter and parallel topology service, and boundary
  // treatment, and molecule service
  if (_boundaryTreatment != NULL) {
    delete _boundaryTreatment;
    _boundaryTreatment = NULL;
  }
  if (_profilePlotter != NULL) {
    delete _profilePlotter;
    _profilePlotter = NULL;
  }
  if (_linkedCellService != NULL) {
    delete _linkedCellService;
    _linkedCellService = NULL;
  }
  if (_moleculeService != NULL) {
    delete _moleculeService;
    _moleculeService = NULL;
  }
  if (_molecularPropertiesService != NULL) {
    delete _molecularPropertiesService;
    _molecularPropertiesService = NULL;
  }
  if (_parallelTopologyService != NULL) {
    delete _parallelTopologyService;
    _parallelTopologyService = NULL;
  }
}

void simplemd::MolecularDynamicsSimulation::simulateOneTimestep(const unsigned int& t) {
  // nop for idle processes
  if (_parallelTopologyService->isIdle()) {
    return;
  }

  // Before this step, the ghost layer contains exactly those particles, which
  // have left the local part of the domain during the previous timestep. The
  // algorithm proceeds as follows:
  //
  // 1. put particles from ghost cells into the correct cells again:
  //   + If the ghost cell is part of a parallel boundary, molecules in the cell
  //   are sent to the respective neighboring process
  //        and sorted into the corresponding inner cell upon unpacking of
  //        receive buffer.
  //   + If the ghost cell is part of a local periodic boundary, molecules in
  //   ghost cells are replicated with adapted position with
  //        respect to periodic boundary conditions, and stored in the local
  //        buffer. Upon unpacking the local buffer, molecules are resorted into
  //        the correct cells.
  //
  // 2. fill ghost cells (for periodic boundaries or parallel boundaries).
  //   After this step, all periodic/ parallel ghost layer cells are
  //   populated with molecules from the respective periodic/ parallel neighbour
  //   cell. In the parallel case, all inner, but near boundary, cells are
  //   broadcasted to each respective ghost layer cell.
  //
  // Depending on whether overlapping is switched on, we enter the respective
  // methods:
  //   + If overlapping is off, perform described steps 1. and 2. via
  //   putBoundaryParticlesToInnerCellsAndFillBoundaryCells and then
  //     compute forces.
  //   + If overlapping is on, perform described steps 1. and 2., compute forces
  //   on inner part of domain, wait for communication
  //     buffers to arrive, then compute forces near boundaries, all within
  //     putBoundaryParticlesToInnerCellsFillBoundaryCellsAndOverlapWithForceComputations.
  if (!_configuration.getSimulationConfiguration().useOverlappingCommunicationWithForceComputation()) {
    _boundaryTreatment->putBoundaryParticlesToInnerCellsAndFillBoundaryCells(_localBoundary, *_parallelTopologyService);
    // compute forces between molecules.
    _linkedCellService->iterateCellPairs(*_lennardJonesForce, false);
  } else {
    _boundaryTreatment->putBoundaryParticlesToInnerCellsFillBoundaryCellsAndOverlapWithForceComputations(_localBoundary, *_parallelTopologyService,
                                                                                                         *_lennardJonesForce, false);
  }

  evaluateStatistics(t);

  _boundaryTreatment->emptyGhostBoundaryCells();

  // plot VTK output
  if ((_configuration.getVTKConfiguration().getWriteEveryTimestep() > 0) && (t % _configuration.getVTKConfiguration().getWriteEveryTimestep() == 0)) {
    _vtkMoleculeWriter->setTimestep(t);
    _moleculeService->iterateMolecules(*_vtkMoleculeWriter, false);
  }

#if BUILD_WITH_ADIOS2
  // plot Adios2 output
  if ((_configuration.getAdios2Configuration().getWriteEveryTimestep() > 0) && (t % _configuration.getAdios2Configuration().getWriteEveryTimestep() == 0)) {
    _Adios2Writer->setTimestep(t);
    _moleculeService->iterateMolecules(*_Adios2Writer, false);
  }
#endif

  // write checkpoint
  if ((_configuration.getCheckpointConfiguration().getWriteEveryTimestep() != 0) &&
      (t % _configuration.getCheckpointConfiguration().getWriteEveryTimestep() == 0)) {
    _moleculeService->writeCheckPoint(*_parallelTopologyService, _checkpointFilestem, t);
  }

  // reorganise memory if needed
  if ((_configuration.getSimulationConfiguration().getReorganiseMemoryEveryTimestep() != 0) &&
      (t % _configuration.getSimulationConfiguration().getReorganiseMemoryEveryTimestep() == 0)) {
    _moleculeService->reorganiseMemory(*_parallelTopologyService, *_linkedCellService);
  }

  if (t % 500 == 0) {
    cellmappings::VaryCheckpointMapping varyCheckpointMapping(
        _configuration.getMoleculeConfiguration().getMass(), _configuration.getDomainConfiguration().getKB(),
        _configuration.getMoleculeConfiguration().getTemperature(), _configuration.getMoleculeConfiguration().getSigma(),
        _configuration.getDomainConfiguration().getMeshWidth());
    _linkedCellService->iterateCells(varyCheckpointMapping, false);
  }

  // empty linked lists
  _linkedCellService->iterateCells(*_emptyLinkedListsMapping, false);

  // time integration. After this step, the velocities and the positions of the
  // molecules have been updated.
  _moleculeService->iterateMolecules(*_timeIntegrator, false);

  // sort molecules into linked cells
  _moleculeService->iterateMolecules(*_updateLinkedCellListsMapping, false);

  if (_parallelTopologyService->getProcessCoordinates() == tarch::la::Vector<MD_DIM, unsigned int>(0)) {
    if (t % 50 == 0 && _localMDSimulation == 0)
      std::cout << "Finish MD timestep " << t << "..." << std::endl;
  }
}

void simplemd::MolecularDynamicsSimulation::runSimulation() {
  // time loop
  for (unsigned int t = 0; t < _configuration.getSimulationConfiguration().getNumberOfTimesteps(); t++) {
    // simulate one timestep
    simulateOneTimestep(t);
  }
  if (_parallelTopologyService->getProcessCoordinates() == tarch::la::Vector<MD_DIM, unsigned int>(0)) {
    std::cout << "Finish MD simulation." << std::endl;
  }
}

void simplemd::MolecularDynamicsSimulation::evaluateStatistics(const unsigned int& t) {
  const unsigned int timeInterval = _configuration.getSimulationConfiguration().computeMacroscopicQuantitiesEveryTimestep();

  if (_configuration.getRDFConfiguration().isDefined()) {
    if (t >= _configuration.getRDFConfiguration().getStartAtTimestep()) {
      if ((t - _configuration.getRDFConfiguration().getStartAtTimestep()) % _configuration.getRDFConfiguration().getEvaluateEveryTimestep() == 0) {
        _linkedCellService->iterateCellPairs(*_rdfMapping, false);
        if ((t - _configuration.getRDFConfiguration().getStartAtTimestep()) % _configuration.getRDFConfiguration().getWriteEveryTimestep() == 0) {
          _rdfMapping->evaluateRDF(_localMDSimulation);
        }
      }
    }
  }

  if ((timeInterval != 0) && (t % timeInterval == 0)) {
    // compute average velocity
    simplemd::cellmappings::ComputeMeanVelocityMapping computeMeanVelocityMapping(*_parallelTopologyService, _localMDSimulation);
    _linkedCellService->iterateCells(computeMeanVelocityMapping, false);
    // compute average temperature
    simplemd::cellmappings::ComputeTemperatureMapping computeTemperatureMapping(*_parallelTopologyService, *_molecularPropertiesService,
                                                                                computeMeanVelocityMapping.getMeanVelocity(), _localMDSimulation);
    _linkedCellService->iterateCells(computeTemperatureMapping, false);
  }

  // trigger profile plotting
  _profilePlotter->accumulateAndPlotInformation(t);
}
