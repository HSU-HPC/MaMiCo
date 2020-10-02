// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_SENDRECV_DATAEXCHANGEFROMMD2MACRO_H_
#define _MOLECULARDYNAMICS_COUPLING_SENDRECV_DATAEXCHANGEFROMMD2MACRO_H_

#include "coupling/CouplingMDDefinitions.h"
#include "coupling/sendrecv/DataExchange.h"
#include "coupling/datastructures/MacroscopicCell.h"
#include "coupling/interface/MacroscopicSolverInterface.h"

namespace coupling {
  namespace sendrecv{
    template<unsigned int dim>
    class DataExchangeFromMD2Macro;
  }
}


/** data exchange from MD, that is the coupling tool, to the macroscopic solver. We transfer the buffers macroscopicMass
 *  and macroscopicMomentum of the MacroscopicCell.
 *  The target ranks are determined by the getRanks()-method of the macroscopic solver interface: the macroscopic
 *  solver interface thus needs to know on which ranks information from the MD simulation are required as input.
 *  The source ranks arise from the unique rank determination within the IndexConversion. We only allow transfer of
 *  non-ghost macroscopic cells to the macroscopic solver, i.e. cells which are completely embedded into the MD domain.
 *  @author Philipp Neumann
 */
template<unsigned int dim>
class coupling::sendrecv::DataExchangeFromMD2Macro:
public coupling::sendrecv::DataExchange<coupling::datastructures::MacroscopicCell<dim>,dim> {

  public:
    DataExchangeFromMD2Macro(
      coupling::interface::MacroscopicSolverInterface<dim>* interface,
      coupling::IndexConversion<dim> * indexConversion,
      unsigned int tagoffset=0
    ):
    coupling::sendrecv::DataExchange<coupling::datastructures::MacroscopicCell<dim>,dim>(TAG_FROM_MD2MACRO+tagoffset),
    _interface(interface), 
    _indexConversion(indexConversion)
    {
      #if (COUPLING_MD_DEBUG==COUPLING_MD_YES)
      std::cout << "DataExchangeFromMD2Macro initialised..." << std::endl;
      #endif
    }
    virtual ~DataExchangeFromMD2Macro(){}


    /** returns the ranks to which a particular cell (at index globalCellIndex) should be sent. */
    virtual std::vector<unsigned int> getTargetRanks(tarch::la::Vector<dim,unsigned int> globalCellIndex){
      // if we need information on macroscopic solver side, return the respective ranks via interface
      if (_interface->receiveMacroscopicQuantityFromMDSolver(globalCellIndex)){
        return _interface->getTargetRanks(globalCellIndex);
      // otherwise return empty vector
      } else {
        return std::vector<unsigned int>();
      }
    }


    /** returns all ranks from which a particular cell (at index globalCellIndex) is sent. */
    virtual std::vector<unsigned int> getSourceRanks(tarch::la::Vector<dim,unsigned int> globalCellIndex){
      std::vector<unsigned int> sourceRanks;
      if (_interface->receiveMacroscopicQuantityFromMDSolver(globalCellIndex)){
        sourceRanks.push_back(_indexConversion->getUniqueRankForMacroscopicCell(globalCellIndex));
      }
      return sourceRanks;
    }


    /** local rule to read from macroscopic cell and write data to (e.g. send) buffer.
     *  We only send the macroscopic mass and macroscopic momentum from MD to the macroscopic solver.
     */
    virtual void readFromCell(double * const buffer, const coupling::datastructures::MacroscopicCell<dim> & cell){
      buffer[0] = cell.getMacroscopicMass();
      for (unsigned int d = 0; d < dim; d++){ buffer[d+1] = cell.getMacroscopicMomentum()[d]; }
    }


    /** local rule to read from receive buffer and write data to macroscopic cell */
    virtual void writeToCell(const double * const buffer, coupling::datastructures::MacroscopicCell<dim>& cell){
      tarch::la::Vector<dim,double> macroscopicMomentum(0.0);
      for (unsigned int d = 0; d < dim; d++){ macroscopicMomentum[d] = buffer[1+d]; }
      cell.setMacroscopicMomentum(macroscopicMomentum);
      cell.setMacroscopicMass(buffer[0]);
    }


    /** returns the number of doubles that are sent per macroscopic cell. */
    virtual unsigned int getDoublesPerCell() const {
      // 1 double: macroscopic mass; dim doubles: macroscopic momentum
      return 1+dim;
    }

    void setIndexConversion(coupling::IndexConversion<dim> * indexConversion) {
      _indexConversion = indexConversion;
    }

  private:
    coupling::interface::MacroscopicSolverInterface<dim> *_interface;
    coupling::IndexConversion<dim> * _indexConversion;
};
#endif // _MOLECULARDYNAMICS_COUPLING_SENDRECV_DATAEXCHANGEFROMMD2MACRO_H_
