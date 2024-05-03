// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#pragma once

// #define DEBUG_FILTER_PIPELINE
#define POST_MULTI_INSTANCE_FILTERING_YES true
#define POST_MULTI_INSTANCE_FILTERING_NO false

// include dependencies
#include "coupling/filtering/sequencing/AsymmetricalFilterJunction.h"
#include "coupling/filtering/sequencing/FilterJunction.h"
#include "coupling/filtering/sequencing/FilterSequence.h"
#include "tarch/configuration/ParseConfiguration.h"
#include "tarch/tinyxml2/tinyxml2.h"

using tarch::configuration::ParseConfiguration;

namespace coupling {
namespace filtering {
template <unsigned int dim> class FilterPipeline;

/*
 * Used as member of FilterPipeline. Displays where that FP is used.
 * per instance:		apply filtering for each MD instance
 * individually, before merging instances post multi instance:	apply filtering
 * after MD instances have been merged
 */
enum class Scope { perInstance, postMultiInstance };
} // namespace filtering
} // namespace coupling

/*
 * Manages different branches of filtering sequences.
 * These filtering sequences may be interdependant by using another's sequences
 * input or completely isolated. As this entire filtering process is applied
 * during MD to Macro communication, it uses the MD simulation's output
 * Macro-Cells as input and output. All configuration is made using an
 * XML-config file and does not require recompilation when modified.
 *
 * @author Felix Maurer
 */
template <unsigned int dim> class coupling::filtering::FilterPipeline {
public:
  FilterPipeline(std::vector<coupling::datastructures::MacroscopicCell<dim>*> inputCells, const coupling::filtering::Scope scope,
                 const tarch::utils::MultiMDService<dim>& multiMDService, const char* cfgpath);

  ~FilterPipeline() {
    for (auto sequence : _sequences)
      delete sequence;
      // TODO: do i have to delete the _...cells as well?

#ifdef DEBUG_FILTER_PIPELINE
    std::cout << "FP: FilterPipeline deconstructed." << std::endl;
#endif
  }

  /*
   * Applies each FilterSequence in order of their appearance in the config
   * file. Ouput of the specified output-FilterSequence will be written to
   * _md2MacroCells.
   */
  void operator()();

  /*
   * Getters for FilterSequences.
   * Not that Junction is a subtype of Sequence, so this is how to get Junctions
   * as well.
   */
  coupling::filtering::FilterSequence<dim>* getSequence(const char* identifier) const;
  std::vector<coupling::filtering::FilterSequence<dim>*> getAllSequences() const { return _sequences; }

#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
  /*
   * Get MPI communicator used by all parallel filters.
   */
  MPI_Comm getFilteringCommunicator() { return _comm; };
#endif

private:
  /*
   * Detects errors in XML config file.
   */
  bool configIsValid(ParseConfiguration::XMLConfiguration& xmlConfig);

  /*
   * Interprets configuration of sequences and intializes them. Parameters
   * known:
   *   -"input": Name of another FilterSequence previously defined (optional,
   * uses MD output (i.e. _md2MacroCells) by default)
   *
   * Also detects which sequence will be used as output to this FilterPipeline.
   */
  void loadSequencesFromXML(tinyxml2::XMLElement* metaNode);

  /*
   * Input cells within the local, md2macro, ghost layer excluding domain
   */
  std::vector<coupling::datastructures::MacroscopicCell<dim>*> _md2MacroCells;
  /*
   * Input cells that do not match the criteria to be in _md2MacroCells.
   */
  std::vector<coupling::datastructures::MacroscopicCell<dim>*> _outerCells;

  ParseConfiguration::XMLConfiguration _config;

  /*
   * Scope in which this FilterPipeline is applied. Cf. coupling::Scope
   * definition.
   */
  const coupling::filtering::Scope _scope;

  std::vector<coupling::filtering::FilterSequence<dim>*> _sequences;

#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
  MPI_Comm _comm;
#endif
};

// include implementation of header
#include "FilterPipeline.cpph"
