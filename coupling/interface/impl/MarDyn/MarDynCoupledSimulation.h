// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef MARDYNCOUPLEDSIMULATION_H_
#define MARDYNCOUPLEDSIMULATION_H_

//MarDyn
#include "Simulation.h"
#include "ensemble/CanonicalEnsemble.h"
#include "ensemble/EnsembleBase.h"
#include "integrators/Integrator.h"
#include "integrators/Leapfrog.h"
#include "io/GeneratorFactory.h"
#include "io/InputOldstyle.h"
#include "io/io.h"
#include "particleContainer/handlerInterfaces/ParticlePairsHandler.h"
#include "utils/Logger.h"
#include "utils/Timer.h"

#ifdef ENABLE_MPI
#include <mpi.h>
#include "parallel/DomainDecomposition.h"
#endif

//MaMico
#include "coupling/services/MacroscopicCellService.h"
#include "coupling/interface/MamicoInterfaceProvider.h"
#include "coupling/interface/Molecule.h"
#include "coupling/interface/MoleculeIterator.h"
#include "coupling/interface/MDSolverInterface.h"
#include "coupling/interface/impl/MarDyn/LinkedCellsForCoupling.h"
#include "coupling/interface/impl/MarDyn/MarDynCell.h"

#include <list>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <limits>

using Log::global_log;

/*	Extends Mardyns Simulation class in order to be used in coupled simulation.
 * 	- initConfigOldStyle: extended for LinkedCellsForCoupling, new DomainDecomposition and removed some none-relevant parts
 * 	- setInitValuesForCoupling(): additional initialization for coupled simulation
 * 	- simulateOneMDTimestep: adapted from the simulate() method to simulate one MD timestep
 * 	- added some get()-methods for Simulation parameters
 * 	- method to remove a molecule from the molecule container
 * 	@author Hanno Flohr
 */
class MarDynCoupledSimulation : public Simulation {

public:
	//	Constructor that sets the needed mamico cell values for LinkedCellsForCoupling and the MPI grid dimensions
	MarDynCoupledSimulation(tarch::la::Vector<3,double> mamicoCellSize, tarch::la::Vector<3,unsigned int> linkedCellsPerMacroscopicCell,
								tarch::la::Vector<3,unsigned int> mpiGridDimensions) :
		Simulation(), _mamicoCellSize(mamicoCellSize), _linkedCellsPerMacroscopicCell(linkedCellsPerMacroscopicCell),
		_macroscopicCellService(NULL), _couplingOn(true), _mpiGridDimensions(mpiGridDimensions), ensemble(NULL)
	{
		global_simulation = this;
	}

	//same as other constructor but without the MPI grid dimensions as input
	MarDynCoupledSimulation(tarch::la::Vector<3,double> mamicoCellSize, tarch::la::Vector<3,unsigned int> linkedCellsPerMacroscopicCell) :
		Simulation(), _mamicoCellSize(mamicoCellSize), _linkedCellsPerMacroscopicCell(linkedCellsPerMacroscopicCell),
		_macroscopicCellService(NULL), _couplingOn(true), _mpiGridDimensions(tarch::la::Vector<3,unsigned int>(0)), ensemble(NULL)
	{
		global_simulation = this;
	}

	virtual ~MarDynCoupledSimulation() {}

	//returns the timestep length used in the simulation
	double getTimestepLength() { return this->_integrator->getTimestepLength(); }

	//returns the type of the particle container
	ParticleContainerType getContainerType() { return _particleContainerType; }

	//returns a pointer to the cell processor that is used in the simulation
	CellProcessor* getCellProcessor() { return _cellProcessor; }

	//returns a pointer to the particle pairs handler that is used in the simulation
	ParticlePairsHandler* getParticlePairsHandler() { return _particlePairsHandler; }

	/* adapted from MarDyn Simulation.cpp, enhanced with search for 'LinkedCellsForCoupling'
	 * and new domain decomposition initialization; some parts not relevant for the coupling removed
	 */
	void initConfigOldstyle(const std::string& inputfilename) {
		global_log->info() << "init oldstyle config file: " << inputfilename << std::endl;

		// open filestream to the input file
		std::ifstream inputfilestream(inputfilename.c_str());
		if (!inputfilestream.is_open()) {
			global_log->error() << "Could not open file " << inputfilename << std::endl;
			exit(1);
		}

		// used to store one token of the inputfilestream
		std::string token;

		double timestepLength;

		// The first line of the config file has to contain the token "MDProjectConfig"
		inputfilestream >> token;
		if ((token != "mardynconfig") && (token != "MDProjectConfig")) {
			global_log->error() << "Not a mardynconfig file! First token: "
					<< token << std::endl;
			exit(1);
		}

		while (inputfilestream) {
			token.clear();
			inputfilestream >> token;
			global_log->debug() << " [[" << token << "]]" << std::endl;

			if (token.substr(0, 1) == "#") {
				inputfilestream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
				continue;
			}
			if (token == "phaseSpaceFile") {
				std::string phaseSpaceFileFormat;
				inputfilestream >> phaseSpaceFileFormat;

				if (timestepLength == 0.0) {
					global_log->error() << "timestep missing." << std::endl;
					exit(1);
				}
				if (phaseSpaceFileFormat == "OldStyle") {
					std::string phaseSpaceFileName;
					inputfilestream >> phaseSpaceFileName;
					_inputReader = (InputBase*) new InputOldstyle();
					_inputReader->setPhaseSpaceFile(phaseSpaceFileName);
					_inputReader->setPhaseSpaceHeaderFile(phaseSpaceFileName);
					_inputReader->readPhaseSpaceHeader(_domain, timestepLength);
				} else if (phaseSpaceFileFormat == "Generator") {
					global_log->info() << "phaseSpaceFileFormat is Generator!"
							<< std::endl;
					std::string generatorName; // name of the library to load
					std::string inputFile; // name of the input file for the generator

					std::string line;
					getline(inputfilestream, line);
					std::stringstream lineStream(line);
					lineStream >> generatorName >> inputFile;
					_inputReader = GeneratorFactory::loadGenerator(generatorName,
							inputFile);
					_inputReader->readPhaseSpaceHeader(_domain, timestepLength);
				} else {
					global_log->error() << "Don't recognize phasespaceFile reader "
							<< phaseSpaceFileFormat << std::endl;
					exit(1);
				}
				if (_LJCutoffRadius == 0.0)
					_LJCutoffRadius = _cutoffRadius;
				_domain->initParameterStreams(_cutoffRadius, _LJCutoffRadius);
			} else if (token == "timestepLength") {
				inputfilestream >> timestepLength;
			} else if (token == "cutoffRadius") {
				inputfilestream >> _cutoffRadius;
			} else if (token == "LJCutoffRadius") {
				inputfilestream >> _LJCutoffRadius;
			} else if ((token == "parallelization") || (token == "parallelisation")) {
	#ifndef ENABLE_MPI
				global_log->warning()
						<< "Input file demands parallelization, but the current compilation doesn't\n\tsupport parallel execution.\n"
						<< std::endl;
				inputfilestream >> token;
	#else
				inputfilestream >> token;
				if (token=="DomainDecomposition") {
					//needed for coupled simulation to assure that Mardyn works with the same grid dimensions as MaMiCo
					if(_mpiGridDimensions[0]>0 && _mpiGridDimensions[1]>0 && _mpiGridDimensions[2]>0) {
						int gridDims[3] {(int)_mpiGridDimensions[0], (int)_mpiGridDimensions[1], (int)_mpiGridDimensions[2]};
						delete _domainDecomposition;
						_domainDecomposition = (DomainDecompBase*) new DomainDecomposition(gridDims);
					}
				}
	#endif
			} else if (token == "datastructure") {

				if (_domainDecomposition == NULL) {
					global_log->error()
							<< "_domainDecomposition is NULL! Probably you compiled for MPI, but didn't specify line \"parallelization\" before line \"datastructure\"!"
							<< std::endl;
					exit(1);
				}

				inputfilestream >> token;
				if (token == "LinkedCells") {
					_particleContainerType = LINKED_CELL;
					int cellsInCutoffRadius;
					inputfilestream >> cellsInCutoffRadius;
					double bBoxMin[3];
					double bBoxMax[3];
					for (int i = 0; i < 3; i++) {
						bBoxMin[i] = _domainDecomposition->getBoundingBoxMin(i,_domain);
						bBoxMax[i] = _domainDecomposition->getBoundingBoxMax(i,_domain);
					}
					if (this->_LJCutoffRadius == 0.0) _LJCutoffRadius = this->_cutoffRadius;
					_moleculeContainer = new LinkedCells(bBoxMin, bBoxMax, _cutoffRadius, _LJCutoffRadius, cellsInCutoffRadius);
				}
				else if(token == "LinkedCellsForCoupling") {
//LINKED CELLS FOR COUPLING
					//check if mamico values are set, abort otherwise
					if(_mamicoCellSize[0]==0.0 || _mamicoCellSize[1]==0.0 || _mamicoCellSize[2]==0.0 ||
						_linkedCellsPerMacroscopicCell[0]==0 || _linkedCellsPerMacroscopicCell[1]==0 || _linkedCellsPerMacroscopicCell[2]==0) {
						global_log->error() << "Input file demands 'LinkedCellsForCoupling' as particle container, "
								"but mamico values are not set. Use the right constructor!" << std::endl;
						exit(1);
					}

					_particleContainerType = LINKED_CELL;
					double bBoxMin[3];
					double bBoxMax[3];
					for(int d=0; d<3; d++) {
						bBoxMin[d] = _domainDecomposition->getBoundingBoxMin(d,_domain);
						bBoxMax[d] = _domainDecomposition->getBoundingBoxMax(d,_domain);
					}
					if(this->_LJCutoffRadius == 0.0) _LJCutoffRadius = this->_cutoffRadius;
					_moleculeContainer = new LinkedCellsForCoupling(_mamicoCellSize, _linkedCellsPerMacroscopicCell, bBoxMin, bBoxMax, _cutoffRadius);
//LINKED CELLS FOR COUPLING
				} else {
					global_log->error() << "UNKOWN DATASTRUCTURE: " << token
							<< std::endl;
					exit(1);
				}
			} else if (token == "output") {
				inputfilestream >> token;
				if (token == "ResultWriter") {
					unsigned long writeFrequency;
					std::string outputPathAndPrefix;
					inputfilestream >> writeFrequency >> outputPathAndPrefix;
					_outputPlugins.push_back(new ResultWriter(writeFrequency,
							outputPathAndPrefix));
					global_log->debug() << "ResultWriter '" << outputPathAndPrefix
							<< "'.\n";
				} else if (token == "XyzWriter") {
					unsigned long writeFrequency;
					std::string outputPathAndPrefix;
					inputfilestream >> writeFrequency >> outputPathAndPrefix;
					_outputPlugins.push_back(new XyzWriter(writeFrequency,
							outputPathAndPrefix, true));
					global_log->debug() << "XyzWriter " << writeFrequency << " '"
							<< outputPathAndPrefix << "'.\n";
				} else if (token == "CheckpointWriter") {
					unsigned long writeFrequency;
					std::string outputPathAndPrefix;
					inputfilestream >> writeFrequency >> outputPathAndPrefix;
					_outputPlugins.push_back(new CheckpointWriter(writeFrequency,
							outputPathAndPrefix, true));
					global_log->debug() << "CheckpointWriter " << writeFrequency
							<< " '" << outputPathAndPrefix << "'.\n";
				} else if (token == "PovWriter") {
					unsigned long writeFrequency;
					std::string outputPathAndPrefix;
					inputfilestream >> writeFrequency >> outputPathAndPrefix;
					_outputPlugins.push_back(new PovWriter(writeFrequency,
							outputPathAndPrefix, true));
					global_log->debug() << "POVWriter " << writeFrequency << " '"
							<< outputPathAndPrefix << "'.\n";
				} else if (token == "DecompWriter") {
					unsigned long writeFrequency;
					std::string mode;
					std::string outputPathAndPrefix;
					inputfilestream >> writeFrequency >> mode
							>> outputPathAndPrefix;
					_outputPlugins.push_back(new DecompWriter(writeFrequency, mode,
							outputPathAndPrefix, true));
					global_log->debug() << "DecompWriter " << writeFrequency
							<< " '" << outputPathAndPrefix << "'.\n";
				} else if ((token == "VisittWriter") || (token == "VISWriter")) {
					unsigned long writeFrequency;
					std::string outputPathAndPrefix;
					inputfilestream >> writeFrequency >> outputPathAndPrefix;
					_outputPlugins.push_back(new VISWriter(writeFrequency,
							outputPathAndPrefix, true));
					global_log->debug() << "VISWriter " << writeFrequency << " '"
							<< outputPathAndPrefix << "'.\n";
				} else if (token == "VTKWriter") {
	#ifdef VTK
					unsigned long writeFrequency = 0;
					std::string outputPathAndPrefix;
					inputfilestream >> writeFrequency >> outputPathAndPrefix;
					_outputPlugins.push_back(new VTKMoleculeWriter(writeFrequency,
							outputPathAndPrefix));
					global_log->debug() << "VTKWriter " << writeFrequency << " '"
							<< outputPathAndPrefix << "'.\n";
	#else
					Log::global_log->error() << std::endl << "VKT-Plotting demanded, but programme compiled without -DVTK!" << std::endl << std::endl;
	#endif
				} else if (token == "VTKGridWriter") {
	#ifdef VTK
					unsigned long writeFrequency = 0;
					std::string outputPathAndPrefix;
					inputfilestream >> writeFrequency >> outputPathAndPrefix;

					if (_particleContainerType == LINKED_CELL) {
						_outputPlugins.push_back(new VTKGridWriter(writeFrequency,
								outputPathAndPrefix));
						global_log->debug() << "VTKGridWriter " << writeFrequency
								<< " '" << outputPathAndPrefix << "'.\n";
					} else {
						global_log->warning()
								<< "VTKGridWriter only supported with LinkedCells!"
								<< std::endl;
						global_log->warning()
								<< "Generating no VTK output for the grid!"
								<< std::endl;
					}
	#else
					Log::global_log->error() << std::endl << "VKT-Plotting demanded, but programme compiled without -DVTK!" << std::endl << std::endl;
	#endif
				} else {
					global_log->warning() << "Unknown output plugin " << token << std::endl;
				}
			} else if (token == "cutoffRadius") {
				double rc;
				inputfilestream >> rc;
				this->setcutoffRadius(rc);
			} else if (token == "LJCutoffRadius") {
				double rc;
				inputfilestream >> rc;
				this->setLJCutoff(rc);
			}
			else {
				if (token != "")
					global_log->warning() << "Did not process unknown token "
							<< token << std::endl;
			}
		}

		// read particle data
		maxid = _inputReader->readPhaseSpace(_moleculeContainer, &_lmu, _domain, _domainDecomposition);

		if (this->_LJCutoffRadius == 0.0)
			_LJCutoffRadius = this->_cutoffRadius;
		_domain->initFarFieldCorr(_cutoffRadius, _LJCutoffRadius);

		// initialize the Leapfrog integrator
		_integrator = new Leapfrog(timestepLength);

		// test new Decomposition
		_moleculeContainer->update();
		_moleculeContainer->deleteOuterParticles();
	}

	//disables the execution of the coupling steps
	void switchOffCoupling() { _couplingOn = false; }

	//enables the execution of the coupling steps
	void switchOnCoupling() { _couplingOn = true; }

	//additional initialization for the coupled step by step simulation
	void setInitValuesForCoupling() {
		//The system is equilibrated at a higher temperature until _simstep==_initCanonical
		//by setting it to zero, no initial temperature increase should happen
		_initCanonical = 0;

		_initSimulation = 0;
		_initStatistics = 0;

		//set macroscopic cell service for coupling actions
		if(_macroscopicCellService==NULL) {
			_macroscopicCellService = (coupling::services::MacroscopicCellServiceImpl<MarDynCell,3>*)
					coupling::interface::MamicoInterfaceProvider<MarDynCell,3>::getInstance().getMacroscopicCellService();
		}

		//set the Ensemble and update global variables
		ensemble = new CanonicalEnsemble(_moleculeContainer, global_simulation->getEnsemble()->components() );
		ensemble->updateGlobalVariable(NUM_PARTICLES);
		global_log->debug() << "Number of particles in the ensemble: " << ensemble->N() << std::endl;
		ensemble->updateGlobalVariable(ENERGY);
		global_log->debug() << "Kinetic energy in the ensemble: " << ensemble->E() << std::endl;
		ensemble->updateGlobalVariable(TEMPERATURE);
		global_log->debug() << "Temperature of the ensemble: " << ensemble->T() << std::endl;

		//turn off Mardyn's thermostat
		_domain->thermostatOff();
	}

	//simulates one MD timestep
	void simulateOneMDTimestep(const unsigned int &timestep) {
		_simstep = timestep;

		global_log->debug() << "timestep: " << getSimulationStep() << std::endl;
		global_log->debug() << "simulation time: " << getSimulationTime() << std::endl;

		//exectue MaMiCo coupling steps (if coupling is turned on)
		if(_couplingOn) executeCouplingSteps();

		//inform integrator about new timestep & first step of leapfrog integrator
		_integrator->eventNewTimestep(_moleculeContainer, _domain);

		//ensure that all particles are in the right cells and exchange particles
		global_log->debug() << "Updating container and decomposition.." << std::endl;
		updateParticleContainerAndDecomposition();

		//force calculation
		global_log->debug() << "Traversing pairs for force calculation.." << std::endl;
		_moleculeContainer->traverseCells(*_cellProcessor);

		//clear halo
		global_log->debug() << "Delete outer particles / clearing halo.." << std::endl;
		_moleculeContainer->deleteOuterParticles();

		//inform integrator about calculated forces & execute steps between force calculation and end of timestep
		global_log->debug() << "Inform integrator about calculated forces.." << std::endl;
		_integrator->eventForcesCalculated(_moleculeContainer, _domain);

		//calculate global macroscopic values from the local values
		global_log->debug() << "Calculate macroscopic values.." << std::endl;
		_domain->calculateGlobalValues(_domainDecomposition, _moleculeContainer, (!(_simstep % _collectThermostatDirectedVelocity)), Tfactor(_simstep) );

		advanceSimulationTime(_integrator->getTimestepLength());

		/* BEGIN PHYSICAL SECTION
		 * the system is in a consistent state, so the global variables can be updated
		 */
		ensemble->updateGlobalVariable(NUM_PARTICLES);
		global_log->debug() << "Number of particles in the ensemble: " << ensemble->N() << std::endl;
		ensemble->updateGlobalVariable(ENERGY);
		global_log->debug() << "Kinetic energy in the ensemble: " << ensemble->E() << std::endl;
		ensemble->updateGlobalVariable(TEMPERATURE);
		global_log->debug() << "Temperature of the ensemble: " << ensemble->T() << std::endl;
		/* END PHYSICAL SECTION */

		//create output for the current timestep, if needed
		output(timestep);
	}


	//removes the molecule with id 'id' from the molecule container
	bool removeMoleculeFromContainer(unsigned long id) {
		bool found = false;
		Molecule* curMolecule;

		//loop through molecules in the particle container
		for(curMolecule = _moleculeContainer->begin(); curMolecule != _moleculeContainer->end(); curMolecule = _moleculeContainer->next()) {
			//if the desired d is found, delete that molecule
			if(curMolecule->id() == id) {
				_moleculeContainer->deleteCurrent();
				found = true;
				break;
			}
		}
		return found;
	}

protected:
	//calls the three coupling methods of the MacroscopicCellService
	void executeCouplingSteps() {
		if(_macroscopicCellService==NULL) {
			_macroscopicCellService = (coupling::services::MacroscopicCellServiceImpl<MarDynCell,3>*)
					coupling::interface::MamicoInterfaceProvider<MarDynCell,3>::getInstance().getMacroscopicCellService();
		}
		_macroscopicCellService->processInnerMacroscopicCellAfterMDTimestep();
                // TODO move distribute momentum after force accumulation (if required)
		_macroscopicCellService->distributeMass(_simstep);
                _macroscopicCellService->distributeMomentum(_simstep);
		_macroscopicCellService->applyTemperatureToMolecules(_simstep);
	}

	//MaMiCo values needed for the initialization of LinkedCellsForCoupling
	tarch::la::Vector<3,double> _mamicoCellSize;
	tarch::la::Vector<3,unsigned int> 	_linkedCellsPerMacroscopicCell;
	//the macroscopic cell service to execute the coupling steps
	coupling::services::MacroscopicCellServiceImpl<MarDynCell,3>* _macroscopicCellService;
	//used to check if coupling steps shall be executed
	bool _couplingOn;
	//grid dimensions for the MPI initialization
	tarch::la::Vector<3,unsigned int> _mpiGridDimensions;

	CanonicalEnsemble* ensemble;
};

#endif /* MARDYNCOUPLEDSIMULATION_H_ */
