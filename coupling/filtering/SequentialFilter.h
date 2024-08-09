// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#pragma once
#include "coupling/CouplingMDDefinitions.h"
#include "coupling/filtering/interfaces/FilterInterface.h"
#include "coupling/indexing/IndexingService.h"
#include <mpi.h>

#define FILTER_SEQUENTIAL true
#define FILTER_PARALLEL false

// #define DEBUG_SEQ_FILTER
// #define DEBUG_SEQ_FILTER_VERBOSE

namespace coupling {
namespace filtering {
template <class Container_T, class CellIndex_T, unsigned int dim> class SequentialFilter;
}
} // namespace coupling

/*
 * Implementation of FilterInterface.h for filters which operate in a sequential
 * manner, i.e. process data on one master rank. For such filters, operator()()
 * will
 * 		- contribute to one dedicated processing rank: by calling
 * contribute()
 * 		- process on master rank only and then scatter data correctly:
 * by calling process()
 *
 * Meant to be used as a wrapper class, i.e. take a pointer to a filter object
 * and then sequentualizes it. Used in
 * coupling::FilterSequence<dim>::loadFiltersFromXML.
 *
 * Disclaimer: Sequential filters be added to a sequence via XML only. This
 * implies FFF with e.g. a python function is not compatible.
 *
 * @author Felix Maurer
 */
template <class Container_T, class CellIndex_T, unsigned int dim> class coupling::filtering::SequentialFilter : public coupling::filtering::FilterInterface<Container_T, dim> {
public:
  SequentialFilter(coupling::filtering::FilterInterface<coupling::datastructures::CellContainer<CellIndex_T, dim>, dim>* filter, const MPI_Comm comm);

  ~SequentialFilter() {
    delete _filter;
    for (auto cell : _inputCells_Global)
      delete cell;
    for (auto cell : _outputCells_Global)
      delete cell;
  }

  /*
   * Implements FilterInterface's requirement of having a ()-operand defined.
   */
  virtual void operator()();

private:
  /*
   * When sequentialized, all ranks call this function.
   */
  virtual void contribute();

  /*
   * When sequentialized, only the processing rank calls this function. It acts
   * as a wrapper of _filter's operator() member function.
   */
  virtual void process(bool sequential);

  /*
   * Auxilliary functions providing an interface between low-level double
   * buffers used by MPI and Coupling Cells.
   */
  void couplingCellToBuffer(std::vector<double>& buf, const coupling::datastructures::CouplingCell<dim>* cell);

  void bufferToCouplingCell(const std::vector<double>& buf, coupling::datastructures::CouplingCell<dim>* cell);

  // The sequentialized Filter
  coupling::filtering::FilterInterface<coupling::datastructures::CellContainer<CellIndex_T, dim>, dim>* _filter;

  // MPI related stuff
  const MPI_Comm _comm;
  int _commSize;
  int _processingRank;
  int _myRank;

  // Globalized variants of cell data structures (i.e spanning across all cells
  // of the global domain). Only the master rank uses these.
  coupling::datastructures::CellContainer<I13, dim> _inputCells_Global;
  coupling::datastructures::CellContainer<I13, dim> _outputCells_Global;

  // Used by the processing rank to remember from which rank it received cells
  // located at which global index TODO: use CellIndex?
  std::vector<unsigned int> _cellRanks;

  // Buffers coupling cells for MPI communication
  std::vector<double> _cellbuf;

  // Used by the processing rank to remember its local domain
  Container_T _inputCells_Local;
  Container_T _outputCells_Local;
};

#include "SequentialFilter.cpph"
