// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef LMP_MAMICO_SORTING_H_
#define LMP_MAMICO_SORTING_H_

#include "mamico_cell.h"
#include "ghost_atoms.h"

#include "coupling/IndexConversion.h"
#include "coupling/CouplingMDDefinitions.h"
#include "tarch/la/Vector.h"

#include "lammps.h"
#include "memory.h"
#include "atom.h"


namespace LAMMPS_NS {


/** this class is used for sorting particles into the MamicoCell structures and to provide access to those.
 *  @author Philipp Neumann
 */
template<unsigned int dim>
class Sorting {
  public:
    enum PrintType {PRINT_ALL_CELLS=0,PRINT_INNER_CELLS=1,PRINT_GHOST_CELLS=2};


    /** initialise "numberCells" mamico cells. This number must be big enough to represent all local Mamico-macroscopic cells (incl.
     *  ghost cells) on every process.
     */
    Sorting(int numberCells,LAMMPS_NS::LAMMPS* lmp):
    _lmp(lmp),
    _numberCells((unsigned int) numberCells),
    _mamicoCells(new MamicoCell[numberCells]),
    _ghostAtoms(lmp)
    {
      if (_mamicoCells == NULL){std::cout << "ERROR Sorting: _mamicoCells==NULL!" << std::endl; exit(EXIT_FAILURE);}
    }

    ~Sorting(){
      if (_mamicoCells != NULL) delete [] _mamicoCells;
      _numberCells = 0;
      _lmp = NULL;
    }


    /** prints the molecules in the mamico cells for debugging purposes */
    void printMolecules(const coupling::IndexConversion<dim>& indexConversion, LAMMPS_NS::Sorting<dim>::PrintType printType) {
      tarch::la::Vector<3,unsigned int> loop(0);
      const tarch::la::Vector<3,unsigned int> end = coupling::initRange<dim>(
        indexConversion.getLocalNumberMacroscopicCells()+tarch::la::Vector<dim,unsigned int>(2)
      );

      for (loop[2] = 0; loop[2] < end[2]; loop[2]++){
      for (loop[1] = 0; loop[1] < end[1]; loop[1]++){
      for (loop[0] = 0; loop[0] < end[0]; loop[0]++){
        //decied whether to output this line or not
        bool decide = true;
        switch(printType){
          case PRINT_ALL_CELLS:
            decide=true; break;
          case PRINT_INNER_CELLS:
            decide=true; for (int d = 0; d < dim; d++){decide=decide&&(loop[d]>0)&&(loop[d]<end[d]-1);} break;
          case PRINT_GHOST_CELLS:
            decide=true; for (int d = 0; d < dim; d++){decide=decide&&(loop[d]>0)&&(loop[d]<end[d]-1);} decide=!decide; break;
          default:
            std::cout << "ERROR printMolecules(): This case should never be reached!" << std::endl; exit(EXIT_FAILURE); break;
        }

        if (decide){
          const unsigned int cellIndex = indexConversion.getLocalCellIndex(coupling::initDimVector<dim>(loop));
          coupling::interface::MoleculeIterator<LAMMPS_NS::MamicoCell,dim> *it = coupling::interface::
            MamicoInterfaceProvider<LAMMPS_NS::MamicoCell,dim>::getInstance().getMDSolverInterface()->getMoleculeIterator(_mamicoCells[cellIndex]);

          for (it->begin(); it->continueIteration(); it->next()){
            const coupling::interface::Molecule<dim>& molecule = it->getConst();
            std::cout << "Rank " << indexConversion.getThisRank() << ", cell " << coupling::initDimVector<dim>(loop) << ", molecule " << molecule.getPosition() << std::endl;
          }
          delete it;
        }
      }
      }
      }
    }



    /** returns the cell at index "index" */
    MamicoCell& getMamicoCell(unsigned int index){ return _mamicoCells[index]; }


    void updateAllCells(const coupling::IndexConversion<dim>& indexConversion){
      // remove all atoms from the cell lists
      #if (COUPLING_MD_DEBUG==COUPLING_MD_YES)
      std::cout << "Flag and reset cells..." << std::endl;
      #endif
      flagAndResetCells(indexConversion);
      // extract ghost atoms from lammps
      #if (COUPLING_MD_DEBUG==COUPLING_MD_YES)
      std::cout << "Extract ghost atoms..." << std::endl;
      #endif
      _ghostAtoms.extractGhostAtoms();
      // update all non-ghost cells
      #if (COUPLING_MD_DEBUG==COUPLING_MD_YES)
      std::cout << "Update non-ghost cells..." << std::endl;
      #endif
      updateNonGhostCells(indexConversion,false);
      #if (COUPLING_MD_DEBUG==COUPLING_MD_YES)
      printMolecules(indexConversion,PRINT_INNER_CELLS);
      #endif
      // update all ghost cells
      #if (COUPLING_MD_DEBUG==COUPLING_MD_YES)
      std::cout << "Update ghost cells..." << std::endl;
      #endif
      updateGhostCells(indexConversion);
    }



    /** removes all atoms from the non-ghost cells and sorts the "nlocal" atoms from this process into the local cells.
      *  If "clearCellLists is set "false", the cell lists are not cleared before the update is carried out; by default, the
      *  cells are emptied.
      */
    void updateNonGhostCells(
      const coupling::IndexConversion<dim>& indexConversion,
      bool clearCellLists=true
    ) {
      #if (COUPLING_MD_DEBUG==COUPLING_MD_YES)
      int rank; MPI_Comm_rank(MPI_COMM_WORLD,&rank);
      #endif
      const tarch::la::Vector<dim,double> localOffset = getLocalOffset(indexConversion);
      const tarch::la::Vector<dim,double> localSize   = getLocalSize(indexConversion);
      // reset all non-ghost cells

      if (clearCellLists){
        const tarch::la::Vector<dim,unsigned int> localCells = indexConversion.getLocalNumberMacroscopicCells()+tarch::la::Vector<dim,unsigned int>(1);
        tarch::la::Vector<3,unsigned int> end = coupling::initRange<dim>(localCells);
        tarch::la::Vector<3,unsigned int> start(0); for (unsigned int d = 0; d < dim; d++) start[d] = 1;
        tarch::la::Vector<3,unsigned int> loop(0);
        for (loop[2] = start[2]; loop[2] < end[2]; loop[2]++){
        for (loop[1] = start[1]; loop[1] < end[1]; loop[1]++){
        for (loop[0] = start[0]; loop[0] < end[0]; loop[0]++){
          tarch::la::Vector<dim,unsigned int> index = coupling::initDimVector<dim>(loop);
          _mamicoCells[indexConversion.getLocalCellIndex(index)].clear();
        }
        }
        }
      }


      // loop over all local atoms and sort them into the cells
      const int nLocalAtoms = _lmp->atom->nlocal;
      for (int n = 0; n < nLocalAtoms; n++){
        // extract atom position
        tarch::la::Vector<dim,double> position(0.0);
        for (unsigned int d = 0; d < dim; d++){ position[d] = _lmp->atom->x[n][d]; }

        // determine global cell index for this atom; convert global to local vector index
        if (isInLocalDomain(position,localOffset,localSize)){
          tarch::la::Vector<dim,unsigned int> vectorIndex = indexConversion.getGlobalVectorCellIndex(position);
          #if (COUPLING_MD_DEBUG==COUPLING_MD_YES)
          std::cout << "Rank " << rank << ": Sort molecule at position " << position << " into global cell index " << vectorIndex << std::endl;
          // check if this is a global non-ghost cell and throw error otherwise
          bool isInnerCell=true;
          for (unsigned int d = 0; d < dim; d++){
            isInnerCell = isInnerCell && (vectorIndex[d]>0) && (vectorIndex[d]<indexConversion.getGlobalNumberMacroscopicCells()[d]+1);
          }
          if (!isInnerCell){ std::cout << "ERROR Sorting::updateNonGhostCells: Molecule is not sorted into an inner cell!" << std::endl;exit(EXIT_FAILURE);}
          #endif
          vectorIndex  = indexConversion.convertGlobalToLocalVectorCellIndex(vectorIndex);
          // add atom to cell
          #if (COUPLING_MD_DEBUG==COUPLING_MD_YES)
          std::cout << "Sort molecule into local cell " << vectorIndex << ", corresponding to " << indexConversion.getLocalCellIndex(vectorIndex) << std::endl;
          #endif
          _mamicoCells[indexConversion.getLocalCellIndex(vectorIndex)].addAtom(n);
        }
      }
    }

    /** forward call to _ghostAtoms */
    double ** const getGhostAtomPositions() const {return _ghostAtoms.getGhostAtomPositions(); }


  private:

    /** sets the ghost flag in all local mamico cells and removes all particle ids from the cells */
    void flagAndResetCells(const coupling::IndexConversion<dim>& indexConversion){
      const tarch::la::Vector<dim,unsigned int> localNumberCells = tarch::la::Vector<dim,unsigned int>(2) + indexConversion.getLocalNumberMacroscopicCells();
      tarch::la::Vector<3,unsigned int> loop(0);
      tarch::la::Vector<3,unsigned int> end = coupling::initRange<dim>(localNumberCells);

      // loop over all local mamico cells
      unsigned int counter=0;
      for (loop[2] = 0; loop[2] < end[2]; loop[2]++){
      for (loop[1] = 0; loop[1] < end[1]; loop[1]++){
      for (loop[0] = 0; loop[0] < end[0]; loop[0]++){
        // determine ghost flag of cell
        bool isInnerCell=true;
        for (unsigned int d = 0; d < dim; d++){ isInnerCell = isInnerCell && (loop[d]>0) && (loop[d] < end[d]-1); }
        #if (COUPLING_MD_DEBUG==COUPLING_MD_YES)
        std::cout << "Flag cell " << loop << " to be ghost cell: " << isInnerCell << std::endl;
        #endif
        _mamicoCells[counter].setGhostCell(!isInnerCell);
        // reset vector with ids in cell
        _mamicoCells[counter].clear();
        counter++;        
      }
      }
      }
    }



    /** sorts the ghost atoms from ghostX into the ghost cells */
    void updateGhostCells(const coupling::IndexConversion<dim>& indexConversion){
      #if (COUPLING_MD_DEBUG==COUPLING_MD_YES)
      int rank; MPI_Comm_rank(MPI_COMM_WORLD,&rank);
      #endif
      // compute bounding box for ghost particles -> if atoms are outside this bounding box, we do not want to consider them anymore
      const tarch::la::Vector<dim,double> localOffset = getLocalOffset(indexConversion);
      const tarch::la::Vector<dim,double> localSize   = getLocalSize(indexConversion);

      double **ghostX = _ghostAtoms.getGhostAtomPositions();
      const int nghost = _ghostAtoms.getNGhost();

      for (int n = 0; n < nghost; n++){
        // extract atom position
        tarch::la::Vector<dim,double> position(0.0);
        for (unsigned int d = 0; d < dim; d++){ position[d] = ghostX[n][d]; }

        // if the current position is in the region of interest, compute cell index and add it to mamico cell
        if (isInLocalDomain(position,localOffset,localSize)){
          // determine global cell index for this atom; convert global to local vector index
          tarch::la::Vector<dim,unsigned int> vectorIndex = indexConversion.convertGlobalToLocalVectorCellIndex(indexConversion.getGlobalVectorCellIndex(position));
          #if (COUPLING_MD_DEBUG==COUPLING_MD_YES)
          std::cout << "Rank " << rank << ": Sort molecule at position " << position << " into global cell index " << indexConversion.convertLocalToGlobalVectorCellIndex(vectorIndex) << std::endl;
          // further check if this is a ghost cell and throw and error otherwise
          bool isGhostCell=false;
          for (int d = 0; d < dim; d++){
            isGhostCell = isGhostCell || (vectorIndex[d]==0) || (vectorIndex[d]==indexConversion.getLocalNumberMacroscopicCells()[d]+1);
          }
          if (!isGhostCell){ std::cout << "ERROR Sorting::updateGhostCells: Molecule is not sorted into a ghost cell!" << std::endl; exit(EXIT_FAILURE);}
          #endif
          // add atom to cell
          _mamicoCells[indexConversion.getLocalCellIndex(vectorIndex)].addAtom(n);
        }
      }
    }

    /** returns the offset of the local MD domain, incl. a ghost layer of mamico cells */
    tarch::la::Vector<dim,double> getLocalOffset(const coupling::IndexConversion<dim>& indexConversion) const {
      // init local offset to global MD offset (very lower left of ghost layer)
      const tarch::la::Vector<dim,double> meshsize = indexConversion.getMacroscopicCellSize();
      tarch::la::Vector<dim,double> localOffset    = indexConversion.getGlobalMDDomainOffset()-meshsize;

      // shift lower offset to the correct process
      const tarch::la::Vector<dim,unsigned int> thisProcess    = indexConversion.getThisProcess();
      const tarch::la::Vector<dim,unsigned int> avgNumberCells = indexConversion.getAverageLocalNumberMacroscopicCells();
      for (unsigned int d = 0; d < dim; d++){ localOffset[d] += meshsize[d]*thisProcess[d]*avgNumberCells[d]; }

      return localOffset;
    }


    /** returns the local domain size of this process, including a ghost layer of mamico cells */
    tarch::la::Vector<dim,double> getLocalSize(const coupling::IndexConversion<dim>& indexConversion) const {
      const tarch::la::Vector<dim,double> meshsize               = indexConversion.getMacroscopicCellSize();
      const tarch::la::Vector<dim,unsigned int> localNumberCells = indexConversion.getLocalNumberMacroscopicCells();
      tarch::la::Vector<dim,double> localSize(0.0);
      for (unsigned int d = 0; d < dim; d++){ localSize[d] = meshsize[d]*(localNumberCells[d]+2); }
      return localSize; 
    }


    /** returns true if the position vector lies inside the local domain, specified by localOffset and localSize */
    bool isInLocalDomain(const tarch::la::Vector<dim,double>& position, const tarch::la::Vector<dim,double>& localOffset,
      const tarch::la::Vector<dim,double>& localSize) const {
      bool isInside=true;
      const tarch::la::Vector<dim,double> upperOffset = localSize+localOffset;
      for (unsigned int d = 0; d < dim; d++){ isInside = isInside && (position[d] >= localOffset[d]) && (position[d] < upperOffset[d]); }
      return isInside;
    }


    LAMMPS_NS::LAMMPS *_lmp;        // points to our lammps instance
    unsigned int _numberCells;      // stores the number of allocated mamico cells
    MamicoCell *_mamicoCells;       // stores all mamico cells
    GhostAtoms<dim> _ghostAtoms;   // stores and manages ghost atoms
};

}

#endif // LMP_MAMICO_SORTING_H_

