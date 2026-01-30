// This file is part of the MaMiCo project. For conditions of distribution
// and use, please see the copyright notice in MaMiCo's main folder

#include <limits>
#include <mpi.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <stdexcept>
#include <string>

#include "coupling/CouplingMDDefinitions.h"
#include "coupling/configurations/MaMiCoConfiguration.h"
#include "coupling/indexing/IndexingService.h"
#include "coupling/interface/MDSimulationFactory.h"
#include "coupling/services/MultiMDCellService.h"
#include "coupling/solvers/CouetteSolverInterface.h"
#include "simplemd/configurations/MolecularDynamicsConfiguration.h"
#include "tarch/configuration/ParseConfiguration.h"
#include "tarch/utils/MultiMDService.h"
#include "tarch/utils/Utils.h"

// for debugging purposes only
#include <iostream>

// #define PYBIND_USE_DMALLOC
//  used for memory debugging
#ifdef PYBIND_USE_DMALLOC
#include "dmalloc.h"
#endif

#if defined(LS1_MARDYN)
#include "coupling/interface/impl/ls1/LS1MDSolverInterface.h"
#include "coupling/interface/impl/ls1/LS1StaticCommData.h"
#include "utils/Logger.h"
using Log::global_log;
#endif

/**
 *  Python bindings for MaMiCo framework
 *
 *  Features:
 *   - Submodule hierarchy structures MaMiCo on the python side
 *   - Interactive runtime access to all XML configuration values
 *   - Automatic bidirectional conversion between python tuple/list and
 * tarch:la::Vector
 *   - MPI parallel (MPI communication on C++ side, rank is returned to python
 * side)
 *   - Interactive runtime management of coupling services and interface classes
 *       (using python wrappers around C++ objects)
 *   - Automatic garbage collection is responsible for deleting objects that are
 *       created / managed on the python side. (Manual memory management on C++
 * side.)
 *
 *  Limitations:
 *   - Implemented for dim=3 only
 *   - Openmpi version minimum required: 3.0
 *
 *  @author Piet Jarmatz
 */

namespace py = pybind11;
using namespace pybind11::literals;

// custom pybind11 headers
#include "coupling/python-binding/conversion.h"

// Helper function for allocation of configuration object and XML parsing
template <class T> T* makeConfiguration(const std::string filename, const std::string topleveltag) {
  T* cfg = new T{}; // corresponding delete will be called by python, because of
                    // return_value_policy::take_ownership
  tarch::configuration::ParseConfiguration::parseConfiguration<T>(filename, topleveltag, *cfg);
  return cfg;
}

///////////////////////////////////////////////////////////////////////////////////
//
// Type casters for automatic runtime conversion between
// PySequence, e.g. tuple or list, on the python side
// tarch:la::Vector on the MaMiCo C++ side
//
///////////////////////////////////////////////////////////////////////////////////

using Vec3ui = tarch::la::Vector<3, unsigned int>;
namespace pybind11 {
namespace detail {
template <> struct type_caster<Vec3ui> {
public:
  PYBIND11_TYPE_CASTER(Vec3ui, _("tarch::la::Vector<3,unsigned int>"));

  /** Conversion part 1 (Python->C++): */
  bool load(handle src, bool implicit) {
    PyObject* source = src.ptr();
    if (!PySequence_Check(source))
      return false;
    if (PySequence_Length(source) != 3)
      return false;

    for (unsigned int index = 0; index < 3; index++) {
      PyObject* itm = PySequence_GetItem(source, index);
      if (!itm)
        return false;

      PyObject* tmp = PyNumber_Long(itm);
      Py_DECREF(itm);
      if (!tmp)
        return false;

      unsigned long val = PyLong_AsUnsignedLong(tmp);
      Py_DECREF(tmp);
      if (PyErr_Occurred())
        return false;

      if (val > std::numeric_limits<unsigned int>::max())
        return false;
      value[index] = (unsigned int)val;
    }
    return true;
  }

  /*** Conversion part 2 (C++ -> Python)  */
  static handle cast(Vec3ui src, return_value_policy /* policy */, handle /* parent */) {
    return PyTuple_Pack(3, PyLong_FromUnsignedLong(src[0]), PyLong_FromUnsignedLong(src[1]), PyLong_FromUnsignedLong(src[2]));
  }
};
} // namespace detail
} // namespace pybind11

using Vec3d = tarch::la::Vector<3, double>;
namespace pybind11 {
namespace detail {
template <> struct type_caster<Vec3d> {
public:
  PYBIND11_TYPE_CASTER(Vec3d, _("tarch::la::Vector<3,double>"));

  /** Conversion part 1 (Python->C++): */
  bool load(handle src, bool implicit) {
    PyObject* source = src.ptr();
    if (!PySequence_Check(source))
      return false;
    if (PySequence_Length(source) != 3)
      return false;

    for (unsigned int index = 0; index < 3; index++) {
      PyObject* itm = PySequence_GetItem(source, index);
      if (!itm)
        return false;
      value[index] = PyFloat_AsDouble(itm);
      Py_DECREF(itm);
      if (PyErr_Occurred())
        return false;
    }
    return true;
  }

  /*** Conversion part 2 (C++ -> Python)  */
  static handle cast(Vec3d src, return_value_policy /* policy */, handle /* parent */) {
    return PyTuple_Pack(3, PyFloat_FromDouble(src[0]), PyFloat_FromDouble(src[1]), PyFloat_FromDouble(src[2]));
  }
};
} // namespace detail
} // namespace pybind11

///////////////////////////////////////////////////////////////////////////////////
//   End of type casters
//   //////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

// Helper function, calls MPI_Init and returns this rank
int initMPI() {
  // Get argc and argv
  py::list argv_list = py::module::import("sys").attr("argv");
  char** argv = new char*[argv_list.size()];
  int argc = 0;
  for (auto arg : argv_list)
    // (todo mini memleak here)
    argv[argc++] = (new std::string(PyUnicode_AsUTF8(arg.ptr())))->data();

  int rank = 0;
  MPI_Init(&argc, &argv);
  delete[] argv;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#if defined(LS1_MARDYN)
  Log::global_log = std::make_unique<Log::Logger>(Log::Error); // Log::Info
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
  global_log->set_mpi_output_root(0);
#endif
#endif
  return rank;
}

////////////////////////////////////////////////////////////////
// Coupling Buffer management //////////////////////////////////
////////////////////////////////////////////////////////////////

class CouplingBuffer {
public:
  CouplingBuffer(coupling::interface::MacroscopicSolverInterface<3>& macroscopicSolverInterface, unsigned int rank, unsigned int outerRegion)
      : _macro2MDBuffer(), _md2macroBuffer(), _outerRegion(outerRegion) {
    allocateMacro2MDBuffer(macroscopicSolverInterface, rank);
    allocateMD2macroBuffer(macroscopicSolverInterface, rank);
  }

  ~CouplingBuffer() {
    deleteBuffer(_macro2MDBuffer);
    deleteBuffer(_md2macroBuffer);
  }

  void recv2CSV(const std::string filename) const {
    std::ofstream file(filename.c_str());
    if (!file.is_open())
      throw std::runtime_error("Could not open file " + filename);
    // loop over received cells; read macroscopic mass+momentum buffers and
    // write cell index, mass and velocity to one line in the csv-file
    I01 idx;
    coupling::datastructures::CouplingCell<3>* cell;
    for (auto pair : _md2macroBuffer) {
      std::tie(cell, idx) = pair;
      tarch::la::Vector<3, double> vel(cell->getMacroscopicMomentum());
      if (cell->getMacroscopicMass() != 0.0) {
        vel = (1.0 / cell->getMacroscopicMass()) * vel;
      }
      file << idx << " ; " << vel[0] << " ; " << vel[1] << " ; " << vel[2] << " ; " << cell->getMacroscopicMass() << ";";
      file << std::endl;
    }
    file.close();
  }

  void store2send(const double cellmass, py::array_t<double> velocity, py::array_t<double> density) {
    // Some sanity checks to see if dimension and shape of
    // velocity and density numpy arrays actually match sendBuffer
    if (velocity.ndim() != 4)
      throw std::runtime_error("velocity.ndim() != 4");
    if (density.ndim() != 3)
      throw std::runtime_error("density.ndim() != 3");
    auto numcells = I09::numberCellsInDomain;
    for (unsigned int d = 0; d < 3; d++) {
      if (velocity.shape(d) != numcells[d])
        throw std::runtime_error("velocity.shape(" + std::to_string(d) + ") should be " + std::to_string(numcells[d]) + " but is " +
                                 std::to_string(velocity.shape(d)));
      if (density.shape(d) != numcells[d])
        throw std::runtime_error("density.shape(" + std::to_string(d) + ") should be " + std::to_string(numcells[d]) + " but is " +
                                 std::to_string(density.shape(d)));
    }
    if (velocity.shape(3) != 3)
      throw std::runtime_error("velocity.shape(3) != 3");

    // Create unchecked proxy objects for direct element access
    // without internal checking of dimensions and bounds
    auto velocity_raw = velocity.unchecked<4>();
    auto density_raw = density.unchecked<3>();

    I01 idx;
    coupling::datastructures::CouplingCell<3>* cell;
    for (auto pair : _macro2MDBuffer) {
      std::tie(cell, idx) = pair;
      const tarch::la::Vector<3, unsigned int> globalIndex{idx.get()};
      double mass = cellmass * density_raw(globalIndex[0], globalIndex[1], globalIndex[2]);

      const tarch::la::Vector<3, double> momentum(mass * velocity_raw(globalIndex[0], globalIndex[1], globalIndex[2], 0),
                                                  mass * velocity_raw(globalIndex[0], globalIndex[1], globalIndex[2], 1),
                                                  mass * velocity_raw(globalIndex[0], globalIndex[1], globalIndex[2], 2));

      cell->setMicroscopicMass(mass);
      cell->setMicroscopicMomentum(momentum);
    }
  }

  py::array_t<double>* loadRecvVelocity() {
    auto numcells = I09::numberCellsInDomain - tarch::la::Vector<3, unsigned int>(2 * _outerRegion);
    std::vector<unsigned int> shape = {numcells[0], numcells[1], numcells[2], 3};
    py::array_t<double>* res = new py::array_t<double>(shape);

    if (numcells[0] * numcells[1] * numcells[2] != _md2macroBuffer.size())
      throw std::runtime_error("Unexpected _md2macroBuffer.size()");

    auto res_raw = res->mutable_unchecked<4>();

    I01 idx_global;
    coupling::datastructures::CouplingCell<3>* cell;
    for (auto pair : _md2macroBuffer) {
      std::tie(cell, idx_global) = pair;
      tarch::la::Vector<3, double> vel(cell->getMacroscopicMomentum());
      if (cell->getMacroscopicMass() != 0.0) {
        vel = (1.0 / cell->getMacroscopicMass()) * vel;
      }
      const tarch::la::Vector<3, unsigned int> idx{I05{idx_global}.get()};
      res_raw(idx[0], idx[1], idx[2], 0) = vel[0];
      res_raw(idx[0], idx[1], idx[2], 1) = vel[1];
      res_raw(idx[0], idx[1], idx[2], 2) = vel[2];
    }

    return res;
  }

  py::array_t<double>* loadRecvDensity(double cellmass) {
    auto numcells = I09::numberCellsInDomain - tarch::la::Vector<3, unsigned int>(2 * _outerRegion);
    std::vector<unsigned int> shape = {numcells[0], numcells[1], numcells[2]};
    py::array_t<double>* res = new py::array_t<double>(shape);

    const unsigned int numCellsRecv = _md2macroBuffer.size();
    if (numcells[0] * numcells[1] * numcells[2] != numCellsRecv)
      throw std::runtime_error("Unexpected _md2macroBuffer.size()");

    auto res_raw = res->mutable_unchecked<3>();

    I01 idx_global;
    coupling::datastructures::CouplingCell<3>* cell;
    for (auto pair : _md2macroBuffer) {
      std::tie(cell, idx_global) = pair;
      const tarch::la::Vector<3, unsigned int> idx{I05{idx_global}.get()};
      res_raw(idx[0], idx[1], idx[2]) = cell->getMacroscopicMass() / cellmass;
    }

    return res;
  }

  coupling::datastructures::FlexibleCellContainer<3> _macro2MDBuffer;
  coupling::datastructures::FlexibleCellContainer<3> _md2macroBuffer;

private:
  const unsigned int _outerRegion; // defines an offset of cells which is
                                   // considered to be the outer region

  /**
   *  @brief allocates the send buffer (with values for all coupling cells).
   *  @param couetteSolverInterface interface for the continuum solver */
  void allocateMacro2MDBuffer(coupling::interface::MacroscopicSolverInterface<3>& msi, int rank) {
    deleteBuffer(_macro2MDBuffer);
    for (auto idx : I08()) {
      if (!I12::contains(idx)) {
        if (tarch::utils::contains(msi.getSourceRanks(idx), (unsigned int)rank)) {
          coupling::datastructures::CouplingCell<3>* couplingCell = new coupling::datastructures::CouplingCell<3>();
          _macro2MDBuffer << std::make_pair(couplingCell, idx);
          if (couplingCell == nullptr)
            throw std::runtime_error(std::string("ERROR mamico.cpp::allocateMacro2MDBuffer: couplingCells==NULL!"));
        }
      }
    }
  }

  /** allocates the recv-buffer. This buffer contains all global inner coupling cells, but only on rank 0. On all other ranks, no cells are stored and a NULL
   * ptr is returned */
  void allocateMD2macroBuffer(coupling::interface::MacroscopicSolverInterface<3>& msi, int rank) {
    deleteBuffer(_md2macroBuffer);
    for (auto idx : I12()) {
      if (tarch::utils::contains(msi.getTargetRanks(idx), (unsigned int)rank)) {
        coupling::datastructures::CouplingCell<3>* couplingCell = new coupling::datastructures::CouplingCell<3>();
        _md2macroBuffer << std::make_pair(couplingCell, idx);
        if (couplingCell == nullptr)
          throw std::runtime_error(std::string("ERROR CouetteScenario::allocateMacro2mdBuffer: couplingCell==NULL!"));
      }
    }
  }

  /** deletes the buffer */
  void deleteBuffer(coupling::datastructures::FlexibleCellContainer<3>& buffer) {
    // delete all potential entries of buffer
    for (auto pair : buffer)
      if (pair.first != NULL)
        delete pair.first;
  }
};

PYBIND11_MODULE(mamico, mamico) {

  ////////////////////////////////////////////////////////////////////////////////
  //  Module and submodule hierarchy definitions
  //  ////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////

  mamico.doc() = "The macro micro coupling tool";

  py::module coupling = mamico.def_submodule("coupling", "Contains the coupling tool and interface implementations");
  py::module tarch = mamico.def_submodule("tarch", "Contains helper classes and files, e.g. for XML parsing");
  py::module simplemd = mamico.def_submodule("simplemd", "Contains the MD simulation SimpleMD");

  py::module configuration = tarch.def_submodule("configuration", "XML configuration utils");
  py::module utils = tarch.def_submodule("utils", "");

  py::module services = coupling.def_submodule("services", "");
  py::module solvers = coupling.def_submodule("solvers", "");
  py::module interface = coupling.def_submodule("interface", "");
  py::module paralleltopology = coupling.def_submodule("paralleltopology", "");

  utils.def("initMPI", &initMPI, "Calls MPI_Init and returns rank of this process");
  utils.def("finalizeMPI", &MPI_Finalize, "Calls MPI_Finalize");
  utils.def(
      "initIndexing",
      [](tarch::la::Vector<3, double> globalMDDomainSize, tarch::la::Vector<3, unsigned int> mdNumberProcesses, tarch::la::Vector<3, double> couplingCellSize,
         coupling::paralleltopology::ParallelTopologyType parallelTopologyType, unsigned int outerRegion, unsigned int rank) {
        const unsigned int dim = 3;
        return IDXS.initWithMDSize(globalMDDomainSize, tarch::la::Vector<3, double>{0, 0, 0}, mdNumberProcesses, couplingCellSize, parallelTopologyType,
                                   outerRegion, rank);
      },
      "Calls init of the IndexingService singleton object");
  ///////////////////////////////////////////////////////////////////////////////////
  // Bindings for all simplemd configuration classes and getter functions
  // //////////
  ///////////////////////////////////////////////////////////////////////////////////

  py::class_<simplemd::configurations::MolecularDynamicsConfiguration>(configuration, "MolecularDynamicsConfiguration")
      .def("__repr__",
           [](const simplemd::configurations::MolecularDynamicsConfiguration& c) {
             return "<MolecularDynamicsConfiguration read from XML tag \"" + c.getTag() + "\">";
           })
      .def("isValid", &simplemd::configurations::MolecularDynamicsConfiguration::isValid)
      .def("getDomainConfiguration", &simplemd::configurations::MolecularDynamicsConfiguration::getDomainConfiguration, py::return_value_policy::reference)
      .def("getMoleculeConfiguration", &simplemd::configurations::MolecularDynamicsConfiguration::getMoleculeConfiguration, py::return_value_policy::reference)
      .def("getMPIConfiguration", &simplemd::configurations::MolecularDynamicsConfiguration::getMPIConfiguration, py::return_value_policy::reference)
      .def("getVTKConfiguration", &simplemd::configurations::MolecularDynamicsConfiguration::getVTKConfiguration, py::return_value_policy::reference)
      .def("getADIOS2Configuration", &simplemd::configurations::MolecularDynamicsConfiguration::getAdios2Configuration, py::return_value_policy::reference)
      .def("getSimulationConfiguration", &simplemd::configurations::MolecularDynamicsConfiguration::getSimulationConfiguration,
           py::return_value_policy::reference)
      .def("getRDFConfiguration", &simplemd::configurations::MolecularDynamicsConfiguration::getRDFConfiguration, py::return_value_policy::reference)
      .def("getCheckpointConfiguration", &simplemd::configurations::MolecularDynamicsConfiguration::getCheckpointConfiguration,
           py::return_value_policy::reference)
      .def("getProfilePlotterConfigurations", &simplemd::configurations::MolecularDynamicsConfiguration::getProfilePlotterConfigurations,
           py::return_value_policy::reference)
      .def("getExternalForceConfigurations", &simplemd::configurations::MolecularDynamicsConfiguration::getExternalForceConfigurations,
           py::return_value_policy::reference);

  py::class_<simplemd::configurations::DomainConfiguration>(configuration, "DomainConfiguration")
      .def("getMoleculesPerDirection", &simplemd::configurations::DomainConfiguration::getMoleculesPerDirection)
      .def("getGlobalDomainSize", &simplemd::configurations::DomainConfiguration::getGlobalDomainSize)
      .def("getGlobalDomainOffset", &simplemd::configurations::DomainConfiguration::getGlobalDomainOffset)
      .def("getCutoffRadius", &simplemd::configurations::DomainConfiguration::getCutoffRadius)
      .def("getMeshWidth", &simplemd::configurations::DomainConfiguration::getMeshWidth)
      .def("getKB", &simplemd::configurations::DomainConfiguration::getKB)
      .def("getBlockSize", &simplemd::configurations::DomainConfiguration::getBlockSize)
      .def("getBoundary", &simplemd::configurations::DomainConfiguration::getBoundary)
      .def("getCheckpointFilestem", &simplemd::configurations::DomainConfiguration::getCheckpointFilestem)
      .def("initFromCheckpoint", &simplemd::configurations::DomainConfiguration::initFromCheckpoint)
      .def("initFromSequentialCheckpoint", &simplemd::configurations::DomainConfiguration::initFromSequentialCheckpoint);

  py::class_<simplemd::configurations::MoleculeConfiguration>(configuration, "MoleculeConfiguration")
      .def("getMeanVelocity", &simplemd::configurations::MoleculeConfiguration::getMeanVelocity)
      .def("getTemperature", &simplemd::configurations::MoleculeConfiguration::getTemperature)
      .def("getMass", &simplemd::configurations::MoleculeConfiguration::getMass)
      .def("getEpsilon", &simplemd::configurations::MoleculeConfiguration::getEpsilon)
      .def("getSigma", &simplemd::configurations::MoleculeConfiguration::getSigma);

  py::class_<simplemd::configurations::MPIConfiguration>(configuration, "MPIConfiguration")
      .def("getNumberOfProcesses", &simplemd::configurations::MPIConfiguration::getNumberOfProcesses);

  py::class_<simplemd::configurations::VTKConfiguration>(configuration, "VTKConfiguration")
      .def("getFilename", &simplemd::configurations::VTKConfiguration::getFilename)
      .def("getWriteEveryTimestep", &simplemd::configurations::VTKConfiguration::getWriteEveryTimestep);

  py::class_<simplemd::configurations::Adios2Configuration>(configuration, "Adios2Configuration")
      .def("getFilename", &simplemd::configurations::Adios2Configuration::getFilename)
      .def("getWriteEveryTimestep", &simplemd::configurations::Adios2Configuration::getWriteEveryTimestep);

  py::class_<simplemd::configurations::SimulationConfiguration>(configuration, "SimulationConfiguration")
      .def("getDt", &simplemd::configurations::SimulationConfiguration::getDt)
      .def("getNumberOfTimesteps", &simplemd::configurations::SimulationConfiguration::getNumberOfTimesteps)
      .def("getReorganiseMemoryEveryTimestep", &simplemd::configurations::SimulationConfiguration::getReorganiseMemoryEveryTimestep)
      .def("computeMacroscopicQuantitiesEveryTimestep", &simplemd::configurations::SimulationConfiguration::computeMacroscopicQuantitiesEveryTimestep)
      .def("fixSeed", &simplemd::configurations::SimulationConfiguration::fixSeed)
      .def("useOverlappingCommunicationWithForceComputation",
           &simplemd::configurations::SimulationConfiguration::useOverlappingCommunicationWithForceComputation);

  py::class_<simplemd::configurations::RDFConfiguration>(configuration, "RDFConfiguration")
      .def("isDefined", &simplemd::configurations::RDFConfiguration::isDefined)
      .def("getStartAtTimestep", &simplemd::configurations::RDFConfiguration::getStartAtTimestep)
      .def("getEvaluateEveryTimestep", &simplemd::configurations::RDFConfiguration::getEvaluateEveryTimestep)
      .def("getWriteEveryTimestep", &simplemd::configurations::RDFConfiguration::getWriteEveryTimestep)
      .def("getNumberOfPoints", &simplemd::configurations::RDFConfiguration::getNumberOfPoints);

  py::class_<simplemd::configurations::CheckpointConfiguration>(configuration, "CheckpointConfiguration");

  py::class_<simplemd::configurations::ProfilePlotterConfiguration>(configuration, "ProfilePlotterConfiguration")
      .def("getStartCell", &simplemd::configurations::ProfilePlotterConfiguration::getStartCell)
      .def("getRange", &simplemd::configurations::ProfilePlotterConfiguration::getRange)
      .def("getWriteEveryTimestep", &simplemd::configurations::ProfilePlotterConfiguration::getWriteEveryTimestep)
      .def("getSampleEveryTimestep", &simplemd::configurations::ProfilePlotterConfiguration::getSampleEveryTimestep)
      .def("getStartAtTimestep", &simplemd::configurations::ProfilePlotterConfiguration::getStartAtTimestep);

  py::class_<simplemd::configurations::ExternalForceConfiguration>(configuration, "ExternalForceConfiguration")
      .def("getExternalForce", &simplemd::configurations::ExternalForceConfiguration::getExternalForce);

  ///////////////////////////////////////////////////////////////////////////////////
  // Bindings for all coupling configuration classes and getter functions
  // //////////
  ///////////////////////////////////////////////////////////////////////////////////

  py::class_<coupling::configurations::MaMiCoConfiguration<3>>(configuration, "MaMiCoConfiguration")
      .def("__repr__",
           [](const coupling::configurations::MaMiCoConfiguration<3>& c) { return "<MaMiCoConfiguration read from XML tag \"" + c.getTag() + "\">"; })
      .def("isValid", &coupling::configurations::MaMiCoConfiguration<3>::isValid)
      .def("getCouplingCellConfiguration", &coupling::configurations::MaMiCoConfiguration<3>::getCouplingCellConfiguration, py::return_value_policy::reference)
      .def("getParticleInsertionConfiguration", &coupling::configurations::MaMiCoConfiguration<3>::getParticleInsertionConfiguration,
           py::return_value_policy::reference)
      .def("getMomentumInsertionConfiguration", &coupling::configurations::MaMiCoConfiguration<3>::getMomentumInsertionConfiguration,
           py::return_value_policy::reference)
      .def("getBoundaryForceConfiguration", &coupling::configurations::MaMiCoConfiguration<3>::getBoundaryForceConfiguration,
           py::return_value_policy::reference)
      .def("getTransferStrategyConfiguration", &coupling::configurations::MaMiCoConfiguration<3>::getTransferStrategyConfiguration,
           py::return_value_policy::reference)
      .def("getParallelTopologyConfiguration", &coupling::configurations::MaMiCoConfiguration<3>::getParallelTopologyConfiguration,
           py::return_value_policy::reference)
      .def("getThermostatConfiguration", &coupling::configurations::MaMiCoConfiguration<3>::getThermostatConfiguration, py::return_value_policy::reference);

  py::class_<coupling::configurations::CouplingCellConfiguration<3>>(configuration, "CouplingCellConfiguration")
      .def("getCouplingCellSize", &coupling::configurations::CouplingCellConfiguration<3>::getCouplingCellSize)
      .def("getNumberLinkedCellsPerCouplingCell", &coupling::configurations::CouplingCellConfiguration<3>::getNumberLinkedCellsPerCouplingCell)
      .def("getWriteEveryMicroscopicTimestep", &coupling::configurations::CouplingCellConfiguration<3>::getWriteEveryMicroscopicTimestep)
      .def("getMicroscopicFilename", &coupling::configurations::CouplingCellConfiguration<3>::getMicroscopicFilename)
      .def("getWriteEveryMacroscopicTimestep", &coupling::configurations::CouplingCellConfiguration<3>::getWriteEveryMacroscopicTimestep)
      .def("getMacroscopicFilename", &coupling::configurations::CouplingCellConfiguration<3>::getMacroscopicFilename);

  py::class_<coupling::configurations::ParticleInsertionConfiguration>(configuration, "ParticleInsertionConfiguration")
      .def("getParticleInsertionType", &coupling::configurations::ParticleInsertionConfiguration::getParticleInsertionType);

  py::class_<coupling::configurations::MomentumInsertionConfiguration<3>>(configuration, "MomentumInsertionConfiguration")
      .def("getMomentumInsertionType", &coupling::configurations::MomentumInsertionConfiguration<3>::getMomentumInsertionType)
      .def("getInnerOverlap", &coupling::configurations::MomentumInsertionConfiguration<3>::getInnerOverlap);

  py::class_<coupling::configurations::BoundaryForceConfiguration<3>>(configuration, "BoundaryForceConfiguration")
      .def("getBoundaryForceType", &coupling::configurations::BoundaryForceConfiguration<3>::getBoundaryForceType);

  py::class_<coupling::configurations::TransferStrategyConfiguration<3>>(configuration, "TransferStrategyConfiguration")
      .def("getStrategyType", &coupling::configurations::TransferStrategyConfiguration<3>::getStrategyType);

  py::class_<coupling::configurations::ParallelTopologyConfiguration>(configuration, "ParallelTopologyConfiguration")
      .def("getParallelTopologyType", &coupling::configurations::ParallelTopologyConfiguration::getParallelTopologyType);

  py::class_<coupling::configurations::ThermostatConfiguration>(configuration, "ThermostatConfiguration")
      .def("getThermostatRegionType", &coupling::configurations::ThermostatConfiguration::getThermostatRegionType)
      .def("getCells2Use", &coupling::configurations::ThermostatConfiguration::getCells2Use);

  configuration.def("parseMolecularDynamicsConfiguration", &makeConfiguration<simplemd::configurations::MolecularDynamicsConfiguration>,
                    py::return_value_policy::take_ownership);
  configuration.def("parseMaMiCoConfiguration", &makeConfiguration<coupling::configurations::MaMiCoConfiguration<3>>, py::return_value_policy::take_ownership);

  ///////////////////////////////////////////////////////////////////////////////////
  // Bindings for coupling service and interface classes and functions
  // //////////////
  ///////////////////////////////////////////////////////////////////////////////////

  py::class_<tarch::utils::MultiMDService<3>>(utils, "MultiMDService")
      .def(py::init<const Vec3ui&, const unsigned int&>(), "numberProcesses"_a = Vec3ui(1, 1, 1), "totalNumberMDSimulations"_a = 1)
      .def("getLocalNumberOfMDSimulations", &tarch::utils::MultiMDService<3>::getLocalNumberOfMDSimulations)
      .def("getGlobalNumberOfLocalMDSimulation", &tarch::utils::MultiMDService<3>::getGlobalNumberOfLocalMDSimulation)
      .def("getLocalCommunicator", [](const tarch::utils::MultiMDService<3>& s) { return (void*)(s.getLocalCommunicator()); });

  py::class_<coupling::interface::MDSimulation>(coupling, "MDSimulation")
      .def("init", (void(coupling::interface::MDSimulation::*)(const tarch::utils::MultiMDService<3>&, unsigned int)) & coupling::interface::MDSimulation::init)
      .def("switchOffCoupling", &coupling::interface::MDSimulation::switchOffCoupling)
      .def("switchOnCoupling", &coupling::interface::MDSimulation::switchOnCoupling)
      .def("simulateTimesteps", &coupling::interface::MDSimulation::simulateTimesteps, "numberTimesteps"_a = 1, "firstTimestep"_a)
      .def("setCouplingCellService", &coupling::interface::MDSimulation::setCouplingCellService)
      .def("shutdown", &coupling::interface::MDSimulation::shutdown);

  py::class_<coupling::services::CouplingCellService<3>>(services, "CouplingCellService")
      .def("computeAndStoreTemperature", &coupling::services::CouplingCellService<3>::computeAndStoreTemperature)
      .def("plotEveryMacroscopicTimestep", &coupling::services::CouplingCellService<3>::plotEveryMacroscopicTimestep)
      .def(
          "addFilterToSequence",
          [](coupling::services::CouplingCellService<3>* service, const char* filter_sequence, int filter_index,
             // use of std::optional is neccessary because pybind11 doesnt
             // support implicit None conversion for STL datatypes
             std::optional<std::function<py::array_t<double>(py::array_t<double>)>> scalar_filter_func,
             std::optional<std::function<py::array_t<double>(py::array_t<double>)>> vector_filter_func) {
            std::function<py::array_t<double>(py::array_t<double>)>* sf_ptr = nullptr;
            if (scalar_filter_func.has_value())
              sf_ptr = &(scalar_filter_func.value());

            std::function<py::array_t<double>(py::array_t<double>)>* vf_ptr = nullptr;
            if (vector_filter_func.has_value())
              vf_ptr = &(vector_filter_func.value());

            // each coupling::conversion:functionWrapper checks for nullptrs,
            // i.e "None" args
            service->getFilterPipeline()
                ->getSequence(filter_sequence)
                ->addFilter(coupling::conversion::functionWrapper_Scalar(sf_ptr), coupling::conversion::functionWrapper_Vector(vf_ptr), filter_index);
          },
          "filter_sequence"_a, "filter_index"_a, "scalar_filter_func"_a = py::none(), "vector_filter_func"_a = py::none());

  coupling.def(
      "getMDSimulation",
      [](const simplemd::configurations::MolecularDynamicsConfiguration& c1, const coupling::configurations::MaMiCoConfiguration<3>& c2, void* comm) {
#if defined(LS1_MARDYN)
        auto offset = c1.getDomainConfiguration().getGlobalDomainOffset();
        coupling::interface::LS1StaticCommData::getInstance().setConfigFilename("ls1config.xml");
        coupling::interface::LS1StaticCommData::getInstance().setBoxOffsetAtDim(0, offset[0]); // temporary till ls1 offset is natively supported
        coupling::interface::LS1StaticCommData::getInstance().setBoxOffsetAtDim(1, offset[1]);
        coupling::interface::LS1StaticCommData::getInstance().setBoxOffsetAtDim(2, offset[2]);

#endif
        return coupling::interface::SimulationAndInterfaceFactory::getInstance().getMDSimulation(c1, c2, (MPI_Comm)comm);
      },
      "simpleMDConfig"_a, "mamicoConfig"_a, "localComm"_a, py::return_value_policy::take_ownership);

  py::class_<coupling::interface::MDSolverInterface<MY_LINKEDCELL, 3>>(interface, "MDSolverInterface");

  coupling.def(
      "getMDSolverInterface",
      [](const simplemd::configurations::MolecularDynamicsConfiguration& c1, const coupling::configurations::MaMiCoConfiguration<3>& c2,
         coupling::interface::MDSimulation* sim) {
        return coupling::interface::SimulationAndInterfaceFactory::getInstance().getMDSolverInterface(c1, c2, sim);
      },
      py::return_value_policy::take_ownership);

  py::class_<coupling::interface::MacroscopicSolverInterface<3>> maSoIf(interface, "MacroscopicSolverInterface");

  py::enum_<coupling::paralleltopology::ParallelTopologyType>(paralleltopology, "ParallelTopologyType")
      .value("UNDEFINED", coupling::paralleltopology::UNDEFINED)
      .value("XYZ", coupling::paralleltopology::XYZ)
      .value("ZYX", coupling::paralleltopology::ZYX)
      .export_values();

  // we indicate that CouetteSolverInterface is a subclass of
  // MacroscopicSolverInterface by passing maSoIf here
  py::class_<coupling::solvers::CouetteSolverInterface<3>>(solvers, "CouetteSolverInterface", maSoIf)
      .def(py::init<tarch::la::Vector<3, unsigned int>, unsigned int>());

  // expose MamicoInterfaceProvider singleton to python side
  coupling.def("setMacroscopicSolverInterface", [](coupling::interface::MacroscopicSolverInterface<3>* iface) {
    coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL, 3>::getInstance().setMacroscopicSolverInterface(iface);
  });
  coupling.def(
      "getMacroscopicSolverInterface",
      []() { return coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL, 3>::getInstance().getMacroscopicSolverInterface(); },
      py::return_value_policy::reference);
  coupling.def("setMDSolverInterface", [](coupling::interface::MDSolverInterface<MY_LINKEDCELL, 3>* iface) {
    coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL, 3>::getInstance().setMDSolverInterface(iface);
  });
  // coupling.def("getMDSolverInterface", [](){
  //		return
  // coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL,3>::getInstance().getMDSolverInterface();
  //	}, py::return_value_policy::reference );
  coupling.def("setCouplingCellService", [](coupling::services::CouplingCellService<3>* service) {
    coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL, 3>::getInstance().setCouplingCellService(service);
  });

  py::class_<coupling::services::MultiMDCellService<MY_LINKEDCELL, 3>>(services, "MultiMDCellService")
      .def(py::init([](std::vector<coupling::interface::MDSolverInterface<MY_LINKEDCELL, 3>*> mdSolverInterfaces,
                       coupling::interface::MacroscopicSolverInterface<3>* macroscopicSolverInterface,
                       simplemd::configurations::MolecularDynamicsConfiguration& simpleMDConfig, unsigned int rank, unsigned int totalNumberMDSimulations,
                       coupling::configurations::MaMiCoConfiguration<3>& mamicoConfig, const std::string filterPipelineConfiguration,
                       tarch::utils::MultiMDService<3>& multiMDService) {
             return new coupling::services::MultiMDCellService<MY_LINKEDCELL, 3>(mdSolverInterfaces, macroscopicSolverInterface, simpleMDConfig, mamicoConfig,
                                                                                 filterPipelineConfiguration.c_str(), multiMDService);
           }),
           "mdSolverInterfaces"_a, "macroscopicSolverInterface"_a, "simpleMDConfig"_a, "rank"_a = 0, "totalNumberMDSimulations"_a = 1, "mamicoConfig"_a,
           "xmlConfigFilename"_a, "multiMDService"_a)
      .def("getCouplingCellService", &coupling::services::MultiMDCellService<MY_LINKEDCELL, 3>::getCouplingCellService, py::return_value_policy::reference)
      .def("sendFromMacro2MD",
           [](coupling::services::MultiMDCellService<MY_LINKEDCELL, 3>& service, CouplingBuffer& buf) { service.sendFromMacro2MD(buf._macro2MDBuffer); })
      .def("sendFromMD2Macro",
           [](coupling::services::MultiMDCellService<MY_LINKEDCELL, 3>& service, CouplingBuffer& buf) { service.sendFromMD2Macro(buf._md2macroBuffer); })
      .def("constructFilterPipelines", &coupling::services::MultiMDCellService<MY_LINKEDCELL, 3>::constructFilterPipelines);

  py::class_<CouplingBuffer>(coupling, "Buffer")
      .def(py::init<coupling::interface::MacroscopicSolverInterface<3>&, unsigned int, unsigned int>())
      .def("recv2CSV", &CouplingBuffer::recv2CSV)
      .def("store2send", &CouplingBuffer::store2send)
      .def("loadRecvVelocity", &CouplingBuffer::loadRecvVelocity, py::return_value_policy::take_ownership)
      .def("loadRecvDensity", &CouplingBuffer::loadRecvDensity, py::return_value_policy::take_ownership);
}
