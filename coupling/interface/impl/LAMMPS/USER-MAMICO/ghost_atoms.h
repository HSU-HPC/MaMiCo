// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef LMP_MAMICO_GHOSTATOMS_H
#define LMP_MAMICO_GHOSTATOMS_H

#include "lammps.h"
#include "atom.h"
#include "memory.h"

#include "coupling/CouplingMDDefinitions.h"

namespace LAMMPS_NS {

/** stores the ghost atom positions in a separate buffer and gives access to that buffer.
 *  @author Philipp Neumann
 */
template<unsigned int dim>
class GhostAtoms {
  public:
    GhostAtoms(LAMMPS_NS::LAMMPS *lmp): _lmp(lmp),_ghostX(NULL),_nghost(0),_nghostBufferSize(0){}
    ~GhostAtoms(){
      if (_ghostX!=NULL) _lmp->memory->destroy(_ghostX);
      _lmp = NULL;
      _nghost = 0;
      _nghostBufferSize = 0;
    }

    /** access the ghost atom buffer */
    double ** const getGhostAtomPositions() const {return _ghostX;}
    /** returns the number of ghost atoms */
    unsigned int getNGhost() const {return _nghost;}

    /** reads all ghost atoms from _lmp and stores them in a separate buffer */
    void extractGhostAtoms(){
      const int nall = _lmp->atom->nlocal + _lmp->atom->nghost;
      const int nlocal = _lmp->atom->nlocal;
      #if (COUPLING_MD_DEBUG==COUPLING_MD_YES)
      int rank; MPI_Comm_rank(MPI_COMM_WORLD,&rank);
      #endif

      // read the new number of ghost atoms
      _nghost = _lmp->atom->nghost;

      // reallocate buffers, if we have more atoms than available in buffer
      if (_nghost > _nghostBufferSize){
        _ghostX = _lmp->memory->grow(_ghostX,_nghost,dim,"mamico:ghostX");
        if (_ghostX==NULL){ std::cout << "ERROR Sorting::extractGhostAtoms: Could not reallocate buffer for " << _nghost << " atoms!" << std::endl; exit(EXIT_FAILURE); }
        // set buffer size
        _nghostBufferSize = _nghost;
      }

      // extract positions of ghost atoms and write them to buffer
      for (int i = nlocal; i < nall; i++){
        for (int d = 0; d < dim; d++){ _ghostX[i-nlocal][d] = _lmp->atom->x[i][d]; }
        #if (COUPLING_MD_DEBUG==COUPLING_MD_YES)
        std::cout << "Rank " << rank << ": Ghost particle "; for (int d = 0; d < dim; d++){ std::cout << _ghostX[i-nlocal][d] << " ";}std::cout << std::endl;
        #endif
      }
    }

  private:
    LAMMPS_NS::LAMMPS* _lmp;        // pointer to lammps instance
    double **_ghostX;               // stores the positions of the ghost atoms; allocated using Lammps memory class
    unsigned int _nghost;           // number of ghost atoms in the buffer
    unsigned int _nghostBufferSize; // size of the buffer; this value is always >= _nghost
};

}
#endif // LMP_MAMICO_GHOSTATOMS_H

