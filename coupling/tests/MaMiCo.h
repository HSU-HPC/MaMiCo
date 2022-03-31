// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MAMICO_H
#define _MAMICO_H

#include <iostream>
#include <string>
#include <mpi.h>

#include "coupling/configurations/MaMiCoConfiguration.h"
#include "coupling/interface/impl/Espresso/EspressoMDSolverInterface.h"
#include "coupling/services/MacroscopicCellService.h"
#include "coupling/solvers/DummySolver.h"
#include "coupling/solvers/DummySolverInterfaceService.h"
#include "tarch/configuration/ParseConfiguration.h"

#include "particle_data.hpp"
#include "domain_decomposition.hpp"
#include "cells.hpp"

static coupling::interface::EspressoMDSolverInterface *
    _espressoMDSolverInterface;
static coupling::services::MacroscopicCellService<3> *_macroscopicCellService;
static DummySolver _dummySolver(18, 18, 18, 0.48);
//static DummySolverInterfaceService *_dummySolverInterfaceService;

void initialize() {
  // Initialize an instance of MDSolverInterface
  if (_espressoMDSolverInterface != NULL) {
    delete _espressoMDSolverInterface;
    _espressoMDSolverInterface = NULL;
  }
  _espressoMDSolverInterface =
      new coupling::interface::EspressoMDSolverInterface();

  // Initialize an instance of MacroscopicSolverInterfaceService
  tarch::la::Vector<3, unsigned int> numberProcesses(2, 2, 1);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  tarch::la::Vector<3, double> globalMDDomainSize(0.0);
  tarch::la::Vector<3, double> globalMDDomainOffset(0.0);
  tarch::la::Vector<3, double> macroscopicCellSize(0.0);
  tarch::la::Vector<3, unsigned int> linkedCellsPerMacroscopicCell(2);
  for (unsigned int d = 0; d < 3; d++) {
    globalMDDomainSize[d] = box_l[d];
    //globalMDDomainOffset(d) = my_left[d];
    globalMDDomainOffset[d] = 0.0;
    macroscopicCellSize[d] = dd.cell_size[d] * linkedCellsPerMacroscopicCell[d];
  }
  tarch::la::Vector<3, unsigned int> nP(1);
  DummySolverInterfaceService::getInstance().init(
      numberProcesses, rank, globalMDDomainSize, globalMDDomainOffset,
      macroscopicCellSize);

  // Initialize configurations for MacroscopicCellService
  const std::string filenameMamico = "mamico_espresso_test_configuration.xml";
  coupling::configurations::MaMiCoConfiguration<3> _configurationMamico;
  tarch::configuration::ParseConfiguration::parseConfiguration<
      coupling::configurations::MaMiCoConfiguration<3> >(
      filenameMamico, "mamico", _configurationMamico);

  unsigned int numberOfTimesteps = 20;

  // Initialize an instance of MacroscopicCellService
  if (_macroscopicCellService != NULL) {
    delete _macroscopicCellService;
    _macroscopicCellService = NULL;
  }
  _macroscopicCellService =
      new coupling::services::MacroscopicCellServiceImpl<ParticleList, 3>(
          0, _espressoMDSolverInterface,
          DummySolverInterfaceService::getInstance().getInterface(),
          numberProcesses, rank,
          _configurationMamico.getParticleInsertionConfiguration(),
          _configurationMamico.getMomentumInsertionConfiguration(),
          _configurationMamico.getBoundaryForceConfiguration(),
          _configurationMamico.getTransferStrategyConfiguration(),
          _configurationMamico.getParallelTopologyConfiguration(),
          numberOfTimesteps,
          _configurationMamico.getMacroscopicCellConfiguration());

  // Read Temperature from the configuration file and set the temperature for
  // the coupled simulation
  _macroscopicCellService->computeAndStoreTemperature(1.8);

}

void shutdown() {
  // shut down interfaces and MD simulation
  if (_macroscopicCellService != NULL) {
    delete _macroscopicCellService;
    _macroscopicCellService = NULL;
  }
  if (_espressoMDSolverInterface != NULL) {
    delete _espressoMDSolverInterface;
    _espressoMDSolverInterface = NULL;
  }
  DummySolverInterfaceService::getInstance().shutdown();
}

// Store Dummy Solver information to send buffers
void storeDummySolverDataInSendBuffer() {
  tarch::la::Vector<3, unsigned int> loop(2);
  unsigned int sendCounter = 0;
  for (loop[2] = 2; loop[2] < 16; loop[2]++) {
    for (loop[1] = 2; loop[1] < 16; loop[1]++) {
      for (loop[0] = 2; loop[0] < 16; loop[0]++) {
        const tarch::la::Vector<3, unsigned int> index =
            coupling::initDimVector<3>(loop);
        const double mass =
            _dummySolver.getDensity(index[0], index[1], index[2]);
        const tarch::la::Vector<3, double> momentum =
            _dummySolver.getVelocity(index[0], index[1], index[2]);
        bool flagsend =
            DummySolverInterfaceService::getInstance().addToSendBuffer(
                mass, momentum, index);
        if (flagsend == true) {
          sendCounter++;
        }
      }
    }
  }

#ifdef MAMICO_DEBUG
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::cout << "On rank " << rank << std::endl;
  std::cout << sendCounter << std::endl;
  const std::vector<coupling::datastructures::MacroscopicCell<3> *>
      _sendBuffer = DummySolverInterfaceService::getInstance().getSendBuffer();
  std::cout << _sendBuffer.size() << std::endl;
  const unsigned int *_globalIndices4SendBuffer =
      DummySolverInterfaceService::getInstance()
          .getGlobalCellIndices4SendBuffer();
  for (unsigned int i = 0; i < 2232; i++) {
    std::cout << _globalIndices4SendBuffer[i] << std::endl;
    std::cout << _sendBuffer[i]->getMicroscopicMass() << " "
              << _sendBuffer[i]->getMicroscopicMomentum() << std::endl;
  }
#endif
}

// Write back data to the dummy solver from the receive buffer
void writeReceiveBufferDataToDummySolver() {
  tarch::la::Vector<3, unsigned int> loop(2);
  unsigned int recvCounter = 0;

#ifdef MAMICO_DEBUG
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::cout << "On rank " << rank << std::endl;
  const std::vector<coupling::datastructures::MacroscopicCell<3> *>
      _recvBuffer =
          DummySolverInterfaceService::getInstance().getReceiveBuffer();
  std::cout << _recvBuffer.size() << std::endl;
  unsigned int *_globalIndices4ReceiveBuffer =
      DummySolverInterfaceService::getInstance()
          .getGlobalCellIndices4ReceiveBuffer();
  for (unsigned int i = 0; i < 512; i++) {
    std::cout << _globalIndices4ReceiveBuffer[i] << std::endl;
    std::cout << _recvBuffer[i]->getMacroscopicMass() << " "
              << _recvBuffer[511]->getMacroscopicMomentum() << std::endl;
  }
#endif

  for (loop[2] = 2; loop[2] < 16; loop[2]++) {
    for (loop[1] = 2; loop[1] < 16; loop[1]++) {
      for (loop[0] = 2; loop[0] < 16; loop[0]++) {
        const tarch::la::Vector<3, unsigned int> index =
            coupling::initDimVector<3>(loop);
        double density;
        tarch::la::Vector<3, double> velocity;
        bool flagrecv =
            DummySolverInterfaceService::getInstance().getFromReceiveBuffer(
                density, velocity, index);
        if (flagrecv == true) {
          recvCounter++;
          _dummySolver.setDensity(density, index[0], index[1], index[2]);
          _dummySolver.setVelocity(velocity, index[0], index[1], index[2]);
        }
      }
    }
  }
}

#endif //_MAMICO_H
