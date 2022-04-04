// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#include "simplemd/services/ParallelAndLocalBufferService.h"

/* Public methods: */

simplemd::services::ParallelAndLocalBufferService::SimpleBuffer::SimpleBuffer(){
  _values = NULL;
  _length = 0;
  _capacity = 0;
}

simplemd::services::ParallelAndLocalBufferService::SimpleBuffer::~SimpleBuffer(){
  /* memory freed during shutdown */
}

bool simplemd::services::ParallelAndLocalBufferService::SimpleBuffer::initialise(
  const unsigned int doublesPerMolecule,
  const unsigned int upperBoundOnNumberOfMolecules)
{
  _capacity = 3*upperBoundOnNumberOfMolecules * doublesPerMolecule;

  _values = (double *) malloc(_capacity * sizeof(double));
  if(_values == NULL) {
    std::cout << "_values of SimpleBuffer, part of the ParallelAndLocalBufferService, could not be allocated. Terminating..." << std::endl;
    return false;
  }
  _length = 0;
  return true;
}

void simplemd::services::ParallelAndLocalBufferService::SimpleBuffer::shutdown()
{
  free(_values);
  _length = 0;
  _capacity = 0;
  _values = NULL;
}

bool simplemd::services::ParallelAndLocalBufferService::SimpleBuffer::pushData(
  const tarch::la::Vector<MD_DIM, double> position,
  const tarch::la::Vector<MD_DIM, double> velocity,
  const tarch::la::Vector<MD_DIM, double> force,
  const double isFixed,
  const bool permitReallocation)
{
  // check if buffer has enough storage for new data
  // we want to add position, velocity and forceOld, so we need 3 * MD_DIM more places

  // check is not in MD_ERROR, because the local buffer can be reallocated and if send/receive buffers are exhausted, we need to terminate in any case!
  if(_length + 3 * MD_DIM + 1 > _capacity) {
    if(permitReallocation == true) {
      //increase buffer capacity via realloc and proceed
      reallocate();
    }
    else {
      std::cout << "Capacity of buffer was exceeded when reallocation is not permitted. Terminating..." << std::endl;
      #if (MD_ERROR == MD_NO) // if MD_ERROR is off, exit here; if it is on, more checks are made, giving info on which processor failed
        exit(EXIT_FAILURE);
      #endif
      return false;
    }
  }

  // push position, velocity, force in this order from first dimension to last
  // unroll loops manually:
  #if (MD_DIM == 1)
    _values[_length]     = position[0];
    _values[_length + 1] = velocity[0];
    _values[_length + 2] = force   [0];
    _values[_length + 3] = isFixed    ;
  #endif
  #if (MD_DIM == 2)
    _values[_length]     = position[0];
    _values[_length + 1] = position[1];
    _values[_length + 2] = velocity[0];
    _values[_length + 3] = velocity[1];
    _values[_length + 4] = force   [0];
    _values[_length + 5] = force   [1];
    _values[_length + 6] = isFixed    ;
  #endif
  #if (MD_DIM == 3)
    _values[_length]     = position[0];
    _values[_length + 1] = position[1];
    _values[_length + 2] = position[2];
    _values[_length + 3] = velocity[0];
    _values[_length + 4] = velocity[1];
    _values[_length + 5] = velocity[2];
    _values[_length + 6] = force   [0];
    _values[_length + 7] = force   [1];
    _values[_length + 8] = force   [2];
    _values[_length + 9] = isFixed    ;
  #endif
  _length += 3 * MD_DIM + 1;

  return true;
}


/* Private methods: */

bool simplemd::services::ParallelAndLocalBufferService::SimpleBuffer::reallocate()
{
  double* temp = NULL;

  unsigned int newsize = _capacity * 2;
  temp = (double*) realloc (_values, newsize * sizeof(double));
  if(temp != NULL) {
    _values = temp;
    _capacity = newsize;
    #if (MD_DEBUG == MD_YES)
      std::cout << "Local buffer reallocated. Buffer capacity is now " << _capacity << " doubles "<< std::endl;
    #endif
    return true;
  }
  else {
    std::cout << "Simple buffer could not be reallocated. Terminating..." << std::endl;
    return false;
  }
}



/* Methods of ParallelAndLocalBufferService */

/* Public methods: */

bool simplemd::services::ParallelAndLocalBufferService::initialise(
  const unsigned int numUniqueNeighbours,
  const unsigned int numCellsPerBuffer[],
  const double avMoleculesPerCell)
{
  bool isOk = true;
  unsigned int doublesPerMolecule = MD_DIM * 3 + 1;
  /* Reallocation of local buffer is permitted, so initialize it with a small value, say 20 molecules */
  unsigned int buffUpperBound = 20;

  isOk = _localBuffer.initialise(doublesPerMolecule, buffUpperBound);
  if(!isOk) {
    std::cout << "Allocation of local buffer with " << buffUpperBound * doublesPerMolecule << " doubles failed. Terminating..." << std::endl;
    return false;
  }

  #if (MD_PARALLEL==MD_YES)
    _numberActiveParallelBuffers = numUniqueNeighbours;
    unsigned int i_buffer;
    for(i_buffer = 0; i_buffer < _numberActiveParallelBuffers; i_buffer ++) {
      buffUpperBound = computeBufferUpperBound(numCellsPerBuffer[i_buffer], avMoleculesPerCell);

      isOk = _sendBuffers[i_buffer].initialise(doublesPerMolecule, buffUpperBound);
      if(!isOk) {
        std::cout << "Allocation of send buffer " << i_buffer << " with " << buffUpperBound  * doublesPerMolecule << " doubles failed. Terminating..." << std::endl;
        return false;
      }
      #if (MD_DEBUG == MD_YES)
        std::cout << "Send buffer " << i_buffer << " was successfully allocated with an upper bound of " << buffUpperBound << std::endl;
      #endif

      isOk = _receiveBuffers[i_buffer].initialise(doublesPerMolecule, buffUpperBound);
      if(!isOk) {
        std::cout << "Allocation of receive buffer " << i_buffer << " with " << buffUpperBound  * doublesPerMolecule << " doubles failed. Terminating..." << std::endl;
        return false;
      }
      #if (MD_DEBUG == MD_YES)
        std::cout << "Receive buffer " << i_buffer << " was successfully allocated with an upper bound of " << buffUpperBound << std::endl;
      #endif
    }
  #endif
  return true;
}

void simplemd::services::ParallelAndLocalBufferService::shutdown()
{
  _localBuffer.shutdown();

  #if (MD_PARALLEL==MD_YES)
  unsigned int i_buffer;
    for(i_buffer = 0; i_buffer < _numberActiveParallelBuffers; i_buffer ++) {
      _sendBuffers[i_buffer].shutdown();
      _receiveBuffers[i_buffer].shutdown();
    }
    _numberActiveParallelBuffers = 0;
  #endif
}

bool simplemd::services::ParallelAndLocalBufferService::pushMoleculeToLocalBuffer(
  const tarch::la::Vector<MD_DIM, double> & position,
  const Molecule* mol)
{
  const bool permitReallocation = true;
  #if (MD_ERROR == MD_YES)
  bool isOk;
  isOk =
  #endif
   _localBuffer.pushData(position, mol->getConstVelocity(), mol->getConstForceOld(), (double)mol->isFixed(), permitReallocation);
  #if (MD_ERROR == MD_YES)
    if(!isOk) {
      std::cout << "Pushing molecule to local buffer failed. Terminating..." << std::endl;
      return false;
    }
  #endif
  return true;
}

#if (MD_PARALLEL == MD_YES)
bool simplemd::services::ParallelAndLocalBufferService::pushMoleculeToSendBuffer(
  const tarch::la::Vector<MD_DIM, double> & position,
  const Molecule* mol,
  const unsigned int i_buffer)
{
  const bool permitReallocation = false;
  #if (MD_ERROR == MD_YES)
  bool isOk;
  isOk =
  #endif
   _sendBuffers[i_buffer].pushData(position, mol->getConstVelocity(), mol->getConstForceOld(), (double)mol->isFixed(), permitReallocation);
  #if (MD_ERROR == MD_YES)
    if(!isOk) {
      std::cout << "Pushing molecule to send buffer " << i_buffer << " failed. Terminating..." << std::endl;
      return false;
    }
  #endif
  return true;
}
#endif


/* Private methods */

unsigned int simplemd::services::ParallelAndLocalBufferService::computeBufferUpperBound(
  const unsigned int numCells,
  const double avMoleculesPerCell) const
{
  // modeling a function of the type A(x) = c * x^-alpha + d
  return (unsigned int)std::ceil( (double)numCells * avMoleculesPerCell * (4.0/std::sqrt((double)numCells) + 1.5) );
}
