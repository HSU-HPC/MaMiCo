// This file is part of the MaMiCo project. For conditions of distribution
// and use, please see the copyright notice in MaMiCo's main folder

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <mpi.h>
#include <string>
#include <limits>
#include <stdexcept>

#include "tarch/utils/MultiMDService.h"
#include "coupling/CouplingMDDefinitions.h"
#include "tarch/configuration/ParseConfiguration.h"
#include "simplemd/configurations/MolecularDynamicsConfiguration.h"
#include "coupling/configurations/MaMiCoConfiguration.h"
#include "coupling/interface/MDSimulationFactory.h"
#include "coupling/services/MultiMDCellService.h"
#include "coupling/solvers/CouetteSolverInterface.h"

// for debugging purposes only
#include <iostream>

/** 
 *  Python bindings for MaMiCo framework
 *
 *  Features:
 *   - Submodule hierarchy structures MaMiCo on the python side 
 *   - Interactive runtime access to all XML configuration values
 *   - Automatic bidirectional conversion between python tuple/list and tarch:la::Vector
 *   - MPI parallel (MPI communication on C++ side, rank is returned to python side)
 *   - Interactive runtime management of coupling services and interface classes
 *       (using python wrappers around C++ objects)
 *   - Automatic garbage collection is responsible for deleting objects that are
 *       created / managed on the python side. (Manual memory management on C++ side.)
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
template<class T> T* makeConfiguration(const std::string filename, const std::string topleveltag){
	T* cfg = new T{};  // corresponding delete will be called by python, because of return_value_policy::take_ownership
	tarch::configuration::ParseConfiguration::parseConfiguration<T>(filename,topleveltag,*cfg);
	return cfg;
}

///////////////////////////////////////////////////////////////////////////////////
// 
// Type casters for automatic runtime conversion between 
// PySequence, e.g. tuple or list, on the python side
// tarch:la::Vector on the MaMiCo C++ side
//
///////////////////////////////////////////////////////////////////////////////////


using Vec3ui = tarch::la::Vector<3,unsigned int>;
namespace pybind11 { namespace detail {
    template <> struct type_caster<Vec3ui> {
    public:
        PYBIND11_TYPE_CASTER(Vec3ui, _("tarch::la::Vector<3,unsigned int>"));

        /** Conversion part 1 (Python->C++): */
        bool load(handle src, bool implicit) {
            PyObject *source = src.ptr();
            if(!PySequence_Check(source)) return false;
            if(PySequence_Length(source) != 3) return false;

            for(unsigned int index = 0; index < 3; index++){
            	PyObject *itm = PySequence_GetItem(source, index);
            	if (!itm) return false;

            	PyObject *tmp = PyNumber_Long(itm);
            	Py_DECREF(itm);
            	if (!tmp)	return false;

            	unsigned long val = PyLong_AsUnsignedLong(tmp);
            	Py_DECREF(tmp);
            	if(PyErr_Occurred()) return false;

            	if(val > std::numeric_limits<unsigned int>::max()) return false;
            	value[index] = (unsigned int)val;
            }
            return true;
        }

        /*** Conversion part 2 (C++ -> Python)  */
        static handle cast(Vec3ui src, return_value_policy /* policy */, handle /* parent */) {
        	return PyTuple_Pack(3, PyLong_FromUnsignedLong(src[0]), PyLong_FromUnsignedLong(src[1]), PyLong_FromUnsignedLong(src[2]));
        }
    };
}} // namespace pybind11::detail

using Vec3d = tarch::la::Vector<3,double>;
namespace pybind11 { namespace detail {
    template <> struct type_caster<Vec3d> {
    public:
        PYBIND11_TYPE_CASTER(Vec3d, _("tarch::la::Vector<3,double>"));

        /** Conversion part 1 (Python->C++): */
        bool load(handle src, bool implicit) {
            PyObject *source = src.ptr();
            if(!PySequence_Check(source)) return false;
            if(PySequence_Length(source) != 3) return false;

            for(unsigned int index = 0; index < 3; index++){
            	PyObject *itm = PySequence_GetItem(source, index);
            	if (!itm) return false;
            	value[index] = PyFloat_AsDouble(itm);
            	Py_DECREF(itm);
            	if(PyErr_Occurred()) return false;
            }
            return true;
        }

        /*** Conversion part 2 (C++ -> Python)  */
        static handle cast(Vec3d src, return_value_policy /* policy */, handle /* parent */) {
        	return PyTuple_Pack(3, PyFloat_FromDouble(src[0]), PyFloat_FromDouble(src[1]), PyFloat_FromDouble(src[2]));
        }
    };
}} // namespace pybind11::detail

///////////////////////////////////////////////////////////////////////////////////
//   End of type casters //////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

// Helper function, calls MPI_Init and returns this rank
int initMPI(){
	// Get argc and argv
	py::list argv_list = py::module::import("sys").attr("argv");
	char** argv = new char * [argv_list.size()];
	int argc = 0;
	for (auto arg : argv_list)
		// (todo mini memleak here)
		argv[argc++] = (new std::string(PyUnicode_AsUTF8(arg.ptr())))->data();  

	int rank = 0;
	MPI_Init(&argc,&argv);
	delete[] argv;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
return rank;
}


////////////////////////////////////////////////////////////////
// Coupling Buffer management //////////////////////////////////
////////////////////////////////////////////////////////////////

class CouplingBuffer{
public:
CouplingBuffer(const coupling::IndexConversion<3>& indexConversion,
	coupling::interface::MacroscopicSolverInterface<3> &macroscopicSolverInterface,
	unsigned int rank, unsigned int outerRegion):
  _sendBuffer(), _globalCellIndices4SendBuffer(), _recvBuffer(), _globalCellIndices4RecvBuffer(),
  _idcv(indexConversion), _outerRegion(outerRegion)
{
	allocateSendBuffer(indexConversion, macroscopicSolverInterface, rank);
	allocateRecvBuffer(indexConversion, macroscopicSolverInterface, rank);
}

~CouplingBuffer(){
	// free buffers/arrays
	deleteBuffer(_sendBuffer);
	if (_globalCellIndices4SendBuffer != NULL){
	  delete [] _globalCellIndices4SendBuffer; 
	  _globalCellIndices4SendBuffer = NULL; 
	}
	deleteBuffer(_recvBuffer);
	if (_globalCellIndices4RecvBuffer != NULL){
	  delete [] _globalCellIndices4RecvBuffer; 
	  _globalCellIndices4RecvBuffer = NULL; 
	}
}

void recv2CSV( const std::string filename) const {
	std::ofstream file(filename.c_str());
	if (!file.is_open())
		throw std::runtime_error("Could not open file " + filename);
	// loop over received cells; read macroscopic mass+momentum buffers and write cell index, mass and velocity to one line in the csv-file
	const unsigned int numCellsRecv = _recvBuffer.size();
	for (unsigned int i = 0; i < numCellsRecv; i++){
	  tarch::la::Vector<3,double> vel(_recvBuffer[i]->getMacroscopicMomentum());
	  if (_recvBuffer[i]->getMacroscopicMass()!=0.0){ vel = (1.0/_recvBuffer[i]->getMacroscopicMass())*vel; }
	  const tarch::la::Vector<3,unsigned int> counter(_idcv.getGlobalVectorCellIndex(_globalCellIndices4RecvBuffer[i]));
	  file << counter[0] << " ; " << counter[1] << " ; " << counter[2] << " ; " << vel[0] << " ; " << vel[1] << " ; " << vel[2] << " ; " << _recvBuffer[i]->getMacroscopicMass() << ";";
	  file << std::endl;
	}
	file.close();
}

void store2send(const double cellmass,
py::array_t<double> velocity, 
py::array_t<double> density) {
	const unsigned int size = _sendBuffer.size();

	// Some sanity checks to see if dimension and shape of
	// velocity and density numpy arrays actually match sendBuffer 
	if(velocity.ndim() != 4)
		throw std::runtime_error("velocity.ndim() != 4");
	if(density.ndim() != 3)
		throw std::runtime_error("density.ndim() != 3");
	auto numcells = _idcv.getGlobalNumberMacroscopicCells();
	for(unsigned int d=0;d<3; d++){
		if(velocity.shape(d) != numcells[d])
			throw std::runtime_error("velocity.shape(" + 
				std::to_string(d) + ") should be " + std::to_string(numcells[d]) +
				" but is " + std::to_string(velocity.shape(d)));
		if(density.shape(d) != numcells[d])
			throw std::runtime_error("density.shape(" + 
				std::to_string(d) + ") should be " + std::to_string(numcells[d]) +
				" but is " + std::to_string(density.shape(d)));
	}
	if(velocity.shape(3) != 3)
		throw std::runtime_error("velocity.shape(3) != 3");

	// Create unchecked proxy objects for direct element access
	// without internal checking of dimensions and bounds
	auto velocity_raw = velocity.unchecked<4>();
	auto density_raw = density.unchecked<3>();
	
	for (unsigned int i = 0; i < size; i++){
		const tarch::la::Vector<3,unsigned int> globalIndex(_idcv.getGlobalVectorCellIndex(_globalCellIndices4SendBuffer[i]) - tarch::la::Vector<3,unsigned int>(1));
		double mass = cellmass
			* density_raw(globalIndex[0],globalIndex[1],globalIndex[2]);

		const tarch::la::Vector<3,double> momentum(
		  mass * velocity_raw(globalIndex[0],globalIndex[1],globalIndex[2],0), 
		  mass * velocity_raw(globalIndex[0],globalIndex[1],globalIndex[2],1),
		  mass * velocity_raw(globalIndex[0],globalIndex[1],globalIndex[2],2) 
		);

		_sendBuffer[i]->setMicroscopicMass(mass);
		_sendBuffer[i]->setMicroscopicMomentum(momentum);
	}
}

py::array_t<double>* loadRecvVelocity(){
  auto numcells = _idcv.getGlobalNumberMacroscopicCells() - tarch::la::Vector<3,unsigned int>(2*_outerRegion);
  std::vector<unsigned int> shape = {numcells[0],numcells[1],numcells[2], 3};
  py::array_t<double>* res = new py::array_t<double>(shape);

  const unsigned int numCellsRecv = _recvBuffer.size();
  if(numcells[0]*numcells[1]*numcells[2] != numCellsRecv)
	throw std::runtime_error("Unexpected _recvBuffer.size()");

  auto res_raw = res->mutable_unchecked<4>();

  for (unsigned int i = 0; i < numCellsRecv; i++){
	tarch::la::Vector<3,double> vel(_recvBuffer[i]->getMacroscopicMomentum());
	if (_recvBuffer[i]->getMacroscopicMass()!=0.0){ vel = (1.0/_recvBuffer[i]->getMacroscopicMass())*vel; }
	const tarch::la::Vector<3,unsigned int> idx(
	  _idcv.getGlobalVectorCellIndex(_globalCellIndices4RecvBuffer[i])
	  - tarch::la::Vector<3,unsigned int>(_outerRegion + 1));
	res_raw(idx[0], idx[1], idx[2], 0) = vel[0];
	res_raw(idx[0], idx[1], idx[2], 1) = vel[1];
	res_raw(idx[0], idx[1], idx[2], 2) = vel[2];
  }

  return res;
}

py::array_t<double>* loadRecvDensity(double cellmass){
  auto numcells = _idcv.getGlobalNumberMacroscopicCells() - tarch::la::Vector<3,unsigned int>(2*_outerRegion);
  std::vector<unsigned int> shape = {numcells[0],numcells[1],numcells[2]};
  py::array_t<double>* res = new py::array_t<double>(shape);

  const unsigned int numCellsRecv = _recvBuffer.size();
  if(numcells[0]*numcells[1]*numcells[2] != numCellsRecv)
	throw std::runtime_error("Unexpected _recvBuffer.size()");

  auto res_raw = res->mutable_unchecked<3>();

  for (unsigned int i = 0; i < numCellsRecv; i++){
	const tarch::la::Vector<3,unsigned int> idx(
	  _idcv.getGlobalVectorCellIndex(_globalCellIndices4RecvBuffer[i])
	  - tarch::la::Vector<3,unsigned int>(_outerRegion + 1));
	res_raw(idx[0], idx[1], idx[2]) = _recvBuffer[i]->getMacroscopicMass() / cellmass;
  }

  return res;
}

std::vector<coupling::datastructures::MacroscopicCell<3>* > _sendBuffer;
unsigned int *_globalCellIndices4SendBuffer;
std::vector<coupling::datastructures::MacroscopicCell<3>* > _recvBuffer;
unsigned int *_globalCellIndices4RecvBuffer;

private:
const coupling::IndexConversion<3>& _idcv;
const unsigned int _outerRegion; // defines an offset of cells which is considered to be the outer region

void allocateSendBuffer(const coupling::IndexConversion<3>& indexConversion,
coupling::interface::MacroscopicSolverInterface<3> &macroscopicSolverInterface,
unsigned int rank) {
  // determine global number of cells
  const tarch::la::Vector<3,unsigned int> cells(indexConversion.getGlobalNumberMacroscopicCells()+tarch::la::Vector<3,unsigned int>(2));
  const unsigned int num = cells[0]*cells[1]*cells[2];

  // delete all potential entries of sendBuffer
  deleteBuffer(_sendBuffer);

  // count number of cells to be sent from this process; therefore, loop over all global macroscopic cells...
  unsigned int numCellsSent=0;
  for (unsigned int i =0; i < num; i++){
	 // ... and find out, if the current cell should be send to MD from this solver process
	 if(macroscopicSolverInterface.sendMacroscopicQuantityToMDSolver(indexConversion.getGlobalVectorCellIndex(i))){
	   std::vector<unsigned int> ranks = macroscopicSolverInterface.getSourceRanks(indexConversion.getGlobalVectorCellIndex(i));
	   bool containsThisRank=false;
	   for (unsigned int k = 0; k < ranks.size(); k++){
		 containsThisRank = containsThisRank || (ranks[k]==rank);
	   }
	   if (containsThisRank){ numCellsSent++; }
	}
  }

  // allocate array for cell indices
  unsigned int* indices = new unsigned int [numCellsSent];

  // allocate sendBuffer and initialise all entries, incl. indices
  for (unsigned int i = 0; i < num; i++){
	if (macroscopicSolverInterface.sendMacroscopicQuantityToMDSolver(indexConversion.getGlobalVectorCellIndex(i))){
	   std::vector<unsigned int> ranks = macroscopicSolverInterface.getSourceRanks(indexConversion.getGlobalVectorCellIndex(i));
	   bool containsThisRank=false;
	   for (unsigned int k = 0; k < ranks.size(); k++){
		 containsThisRank = containsThisRank || (ranks[k]==rank);
	   }
	   if (containsThisRank){
		 _sendBuffer.push_back(new coupling::datastructures::MacroscopicCell<3>());
		 indices[_sendBuffer.size()-1] = i;
	   }
	}
  }
  _globalCellIndices4SendBuffer = indices;
}

void allocateRecvBuffer(const coupling::IndexConversion<3>& indexConversion,
coupling::interface::MacroscopicSolverInterface<3> &macroscopicSolverInterface,
unsigned int rank) {

  // determine global number of cells
  const tarch::la::Vector<3,unsigned int> cells(indexConversion.getGlobalNumberMacroscopicCells()+tarch::la::Vector<3,unsigned int>(2));
  const unsigned int num = cells[0]*cells[1]*cells[2];

  // delete all potential entries of sendBuffer
  deleteBuffer(_recvBuffer);

  // determine number of cells that should be received
  unsigned int numCellsRecv = 0;
  for (unsigned int i = 0; i < num; i++){
	if(macroscopicSolverInterface.receiveMacroscopicQuantityFromMDSolver(indexConversion.getGlobalVectorCellIndex(i))){
	  std::vector<unsigned int> ranks = macroscopicSolverInterface.getTargetRanks(indexConversion.getGlobalVectorCellIndex(i));
	  bool containsThisRank=false;
	  for (unsigned int k = 0; k < ranks.size(); k++){
		containsThisRank = containsThisRank || (ranks[k]==rank);
	  }
	  if (containsThisRank){ numCellsRecv++; }
	}
  }
  // allocate array for cell indices
  unsigned int* indices = new unsigned int [numCellsRecv];

  // allocate recvBuffer and initialise all entries, incl. indices
  for (unsigned int i = 0; i < num; i++){
	if (macroscopicSolverInterface.receiveMacroscopicQuantityFromMDSolver(indexConversion.getGlobalVectorCellIndex(i))){
	  std::vector<unsigned int> ranks = macroscopicSolverInterface.getTargetRanks(indexConversion.getGlobalVectorCellIndex(i));
	  bool containsThisRank=false;
	  for (unsigned int k = 0; k < ranks.size(); k++){
		containsThisRank = containsThisRank || (ranks[k]==rank);
	  }
	  if (containsThisRank){
		_recvBuffer.push_back(new coupling::datastructures::MacroscopicCell<3>());
		// set linearized index
		indices[_recvBuffer.size()-1] = i;
	  }
	}
  }
  _globalCellIndices4RecvBuffer = indices;
}

/** deletes the buffer */
void deleteBuffer(std::vector<coupling::datastructures::MacroscopicCell<3>* >& buffer) {
	// delete all potential entries of buffer
	for (unsigned int i = 0; i < buffer.size(); i++){if (buffer[i]!=NULL){ delete buffer[i]; buffer[i]=NULL;}}
	buffer.clear();
}
};

PYBIND11_MODULE(mamico, mamico) {

////////////////////////////////////////////////////////////////////////////////
//  Module and submodule hierarchy definitions  ////////////////////////////////
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

utils.def("initMPI", &initMPI, "Calls MPI_Init and returns rank of this process");
utils.def("finalizeMPI", &MPI_Finalize, "Calls MPI_Finalize");

///////////////////////////////////////////////////////////////////////////////////
// Bindings for all simplemd configuration classes and getter functions  //////////
///////////////////////////////////////////////////////////////////////////////////

py::class_<simplemd::configurations::MolecularDynamicsConfiguration>(configuration, "MolecularDynamicsConfiguration")
	.def("__repr__",
		[](const simplemd::configurations::MolecularDynamicsConfiguration &c) {
			return "<MolecularDynamicsConfiguration read from XML tag \"" + c.getTag() + "\">";
		}
	)
	.def("isValid", &simplemd::configurations::MolecularDynamicsConfiguration::isValid)
	.def("getDomainConfiguration", &simplemd::configurations::MolecularDynamicsConfiguration::getDomainConfiguration, py::return_value_policy::reference)
	.def("getMoleculeConfiguration", &simplemd::configurations::MolecularDynamicsConfiguration::getMoleculeConfiguration, py::return_value_policy::reference)
	.def("getMPIConfiguration", &simplemd::configurations::MolecularDynamicsConfiguration::getMPIConfiguration, py::return_value_policy::reference)
	.def("getVTKConfiguration", &simplemd::configurations::MolecularDynamicsConfiguration::getVTKConfiguration, py::return_value_policy::reference)
	.def("getSimulationConfiguration", &simplemd::configurations::MolecularDynamicsConfiguration::getSimulationConfiguration, py::return_value_policy::reference)
	.def("getRDFConfiguration", &simplemd::configurations::MolecularDynamicsConfiguration::getRDFConfiguration, py::return_value_policy::reference)
	.def("getCheckpointConfiguration", &simplemd::configurations::MolecularDynamicsConfiguration::getCheckpointConfiguration, py::return_value_policy::reference)
	.def("getProfilePlotterConfigurations", &simplemd::configurations::MolecularDynamicsConfiguration::getProfilePlotterConfigurations, py::return_value_policy::reference)
	.def("getExternalForceConfigurations", &simplemd::configurations::MolecularDynamicsConfiguration::getExternalForceConfigurations, py::return_value_policy::reference);

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

py::class_<simplemd::configurations::SimulationConfiguration>(configuration, "SimulationConfiguration")	
	.def("getDt", &simplemd::configurations::SimulationConfiguration::getDt)
	.def("getNumberOfTimesteps", &simplemd::configurations::SimulationConfiguration::getNumberOfTimesteps)
	.def("getReorganiseMemoryEveryTimestep", &simplemd::configurations::SimulationConfiguration::getReorganiseMemoryEveryTimestep)
	.def("computeMacroscopicQuantitiesEveryTimestep", &simplemd::configurations::SimulationConfiguration::computeMacroscopicQuantitiesEveryTimestep)
	.def("fixSeed", &simplemd::configurations::SimulationConfiguration::fixSeed)
	.def("useOverlappingCommunicationWithForceComputation", &simplemd::configurations::SimulationConfiguration::useOverlappingCommunicationWithForceComputation);

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
// Bindings for all coupling configuration classes and getter functions  //////////
///////////////////////////////////////////////////////////////////////////////////

py::class_<coupling::configurations::MaMiCoConfiguration<3>>(configuration, "MaMiCoConfiguration")
	.def("__repr__",
		[](const coupling::configurations::MaMiCoConfiguration<3> &c) {
			return "<MaMiCoConfiguration read from XML tag \"" + c.getTag() + "\">";
		}
	)
	.def("isValid", &coupling::configurations::MaMiCoConfiguration<3>::isValid)
	.def("getMacroscopicCellConfiguration", &coupling::configurations::MaMiCoConfiguration<3>::getMacroscopicCellConfiguration, py::return_value_policy::reference)
	.def("getParticleInsertionConfiguration", &coupling::configurations::MaMiCoConfiguration<3>::getParticleInsertionConfiguration, py::return_value_policy::reference)
	.def("getMomentumInsertionConfiguration", &coupling::configurations::MaMiCoConfiguration<3>::getMomentumInsertionConfiguration, py::return_value_policy::reference)
	.def("getBoundaryForceConfiguration", &coupling::configurations::MaMiCoConfiguration<3>::getBoundaryForceConfiguration, py::return_value_policy::reference)
	.def("getTransferStrategyConfiguration", &coupling::configurations::MaMiCoConfiguration<3>::getTransferStrategyConfiguration, py::return_value_policy::reference)
	.def("getParallelTopologyConfiguration", &coupling::configurations::MaMiCoConfiguration<3>::getParallelTopologyConfiguration, py::return_value_policy::reference);

py::class_<coupling::configurations::MacroscopicCellConfiguration<3>>(configuration, "MacroscopicCellConfiguration")
	.def("getMacroscopicCellSize", &coupling::configurations::MacroscopicCellConfiguration<3>::getMacroscopicCellSize)
	.def("getNumberLinkedCellsPerMacroscopicCell", &coupling::configurations::MacroscopicCellConfiguration<3>::getNumberLinkedCellsPerMacroscopicCell)
	.def("getWriteEveryMicroscopicTimestep", &coupling::configurations::MacroscopicCellConfiguration<3>::getWriteEveryMicroscopicTimestep)
	.def("getMicroscopicFilename", &coupling::configurations::MacroscopicCellConfiguration<3>::getMicroscopicFilename)
	.def("getWriteEveryMacroscopicTimestep", &coupling::configurations::MacroscopicCellConfiguration<3>::getWriteEveryMacroscopicTimestep)
	.def("getMacroscopicFilename", &coupling::configurations::MacroscopicCellConfiguration<3>::getMacroscopicFilename);

py::class_<coupling::configurations::ParticleInsertionConfiguration>(configuration, "ParticleInsertionConfiguration")
	.def("getParticleInsertionType", &coupling::configurations::ParticleInsertionConfiguration::getParticleInsertionType);

py::class_<coupling::configurations::MomentumInsertionConfiguration>(configuration, "MomentumInsertionConfiguration")
	.def("getMomentumInsertionType", &coupling::configurations::MomentumInsertionConfiguration::getMomentumInsertionType)
	.def("getInnerOverlap", &coupling::configurations::MomentumInsertionConfiguration::getInnerOverlap);

py::class_<coupling::configurations::BoundaryForceConfiguration<3>>(configuration, "BoundaryForceConfiguration")
	.def("getBoundaryForceType", &coupling::configurations::BoundaryForceConfiguration<3>::getBoundaryForceType);

py::class_<coupling::configurations::TransferStrategyConfiguration<3>>(configuration, "TransferStrategyConfiguration")
	.def("getStrategyType", &coupling::configurations::TransferStrategyConfiguration<3>::getStrategyType);

py::class_<coupling::configurations::ParallelTopologyConfiguration>(configuration, "ParallelTopologyConfiguration")
	.def("getParallelTopologyType", &coupling::configurations::ParallelTopologyConfiguration::getParallelTopologyType);

configuration.def("parseMolecularDynamicsConfiguration", &makeConfiguration<simplemd::configurations::MolecularDynamicsConfiguration>, py::return_value_policy::take_ownership);
configuration.def("parseMaMiCoConfiguration", &makeConfiguration<coupling::configurations::MaMiCoConfiguration<3>>, py::return_value_policy::take_ownership);

///////////////////////////////////////////////////////////////////////////////////
// Bindings for coupling service and interface classes and functions //////////////
///////////////////////////////////////////////////////////////////////////////////

py::class_<tarch::utils::MultiMDService<3>>(utils, "MultiMDService")
	.def(py::init<const Vec3ui &, const unsigned int &>(),
		"numberProcesses"_a=Vec3ui(1,1,1), 
		"totalNumberMDSimulations"_a=1)
	.def("getLocalNumberOfMDSimulations", &tarch::utils::MultiMDService<3>::getLocalNumberOfMDSimulations)
	.def("getGlobalNumberOfLocalMDSimulation", &tarch::utils::MultiMDService<3>::getGlobalNumberOfLocalMDSimulation)
	.def("getLocalCommunicator", [](const tarch::utils::MultiMDService<3>& s){return (void*)(s.getLocalCommunicator());} );

py::class_<coupling::interface::MDSimulation>(coupling, "MDSimulation")
	.def("init", (void (coupling::interface::MDSimulation::*)(const tarch::utils::MultiMDService<3>&,unsigned int))
		&coupling::interface::MDSimulation::init)
	.def("switchOffCoupling", &coupling::interface::MDSimulation::switchOffCoupling)
	.def("switchOnCoupling", &coupling::interface::MDSimulation::switchOnCoupling)
	.def("simulateTimesteps", &coupling::interface::MDSimulation::simulateTimesteps, 
		"numberTimesteps"_a=1, "firstTimestep"_a)
	.def("setMacroscopicCellService", &coupling::interface::MDSimulation::setMacroscopicCellService)
	.def("shutdown", &coupling::interface::MDSimulation::shutdown); 

py::class_<coupling::services::MacroscopicCellService<3>>(services, "MacroscopicCellService")
	.def("computeAndStoreTemperature", &coupling::services::MacroscopicCellService<3>::computeAndStoreTemperature)
	.def("getIndexConversion", &coupling::services::MacroscopicCellService<3>::getIndexConversion, py::return_value_policy::reference)
	.def("plotEveryMacroscopicTimestep", &coupling::services::MacroscopicCellService<3>::plotEveryMacroscopicTimestep)
	.def("addFilterToSequence", []
			(
				 coupling::services::MacroscopicCellService<3>* service,
				 const char* filter_sequence,
				 int filter_index,
				 //use of std::optional is neccessary because pybind11 doesnt support implicit None conversion for STL datatypes
				 std::optional<std::function<py::array_t<double> (py::array_t<double>)>> scalar_filter_func,
				 std::optional<std::function<py::array_t<double> (py::array_t<double>)>> vector_filter_func
			)
			{
				std::function<py::array_t<double> (py::array_t<double>)> * sf_ptr = nullptr;
				if(scalar_filter_func.has_value()) sf_ptr = &(scalar_filter_func.value());

				std::function<py::array_t<double> (py::array_t<double>)> * vf_ptr = nullptr;
				if(vector_filter_func.has_value()) vf_ptr = &(vector_filter_func.value());

				//each coupling::conversion:functionWrapper checks for nullptrs, i.e "None" args
				service->addFilterToSequence( 
						filter_sequence, 
						coupling::conversion::functionWrapper_Scalar(sf_ptr),
						coupling::conversion::functionWrapper_Vector(vf_ptr),
						filter_index);
			}, 
			"filter_sequence"_a,
			"filter_index"_a,
			"scalar_filter_func"_a = py::none(),
			"vector_filter_func"_a = py::none()
			);


    coupling.def("getMDSimulation", []
    	(const simplemd::configurations::MolecularDynamicsConfiguration& c1,
    	 const coupling::configurations::MaMiCoConfiguration<3>& c2,
    	 void* comm){
    		return coupling::interface::SimulationAndInterfaceFactory::getInstance().getMDSimulation(c1,c2,(MPI_Comm) comm);
    	}, "simpleMDConfig"_a, "mamicoConfig"_a, "localComm"_a, py::return_value_policy::take_ownership);

    py::class_<coupling::interface::MDSolverInterface<MY_LINKEDCELL,3>>(interface, "MDSolverInterface");

    coupling.def("getMDSolverInterface", []
    	(const simplemd::configurations::MolecularDynamicsConfiguration& c1,
    	 const coupling::configurations::MaMiCoConfiguration<3>& c2,
    	 coupling::interface::MDSimulation* sim){
    		return coupling::interface::SimulationAndInterfaceFactory::getInstance().getMDSolverInterface(c1,c2,sim);
    	},  py::return_value_policy::take_ownership);

    py::class_<coupling::interface::MacroscopicSolverInterface<3>> maSoIf(interface, "MacroscopicSolverInterface");

    // we indicate that CouetteSolverInterface is a subclass of MacroscopicSolverInterface by passing maSoIf here
    py::class_<coupling::solvers::CouetteSolverInterface<3>>(solvers, "CouetteSolverInterface", maSoIf)
    	.def(py::init<tarch::la::Vector<3,unsigned int>,unsigned int>());

    // expose MamicoInterfaceProvider singleton to python side
    coupling.def("setMacroscopicSolverInterface", []
    	(coupling::interface::MacroscopicSolverInterface<3> * iface){
    		coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL,3>::getInstance().setMacroscopicSolverInterface(iface);
    	} );
    coupling.def("getMacroscopicSolverInterface", [](){
    		return coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL,3>::getInstance().getMacroscopicSolverInterface();
    	}, py::return_value_policy::reference );
    coupling.def("setMDSolverInterface", []
    	(coupling::interface::MDSolverInterface<MY_LINKEDCELL,3> * iface){
    		coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL,3>::getInstance().setMDSolverInterface(iface);
    	});
    //coupling.def("getMDSolverInterface", [](){
    //		return coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL,3>::getInstance().getMDSolverInterface();
    //	}, py::return_value_policy::reference );
    coupling.def("setMacroscopicCellService", []
        (coupling::services::MacroscopicCellService<3>* service){
            coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL,3>::getInstance().setMacroscopicCellService(service);
        });

    py::class_<coupling::services::MultiMDCellService<MY_LINKEDCELL,3>>(services, "MultiMDCellService")
    	.def(py::init([](
	    		std::vector<coupling::interface::MDSolverInterface<MY_LINKEDCELL,3>*> mdSolverInterfaces,
			    coupling::interface::MacroscopicSolverInterface<3>* macroscopicSolverInterface,
			    const simplemd::configurations::MolecularDynamicsConfiguration& simpleMDConfig,
			    unsigned int rank,
			    unsigned int totalNumberMDSimulations,
			    const coupling::configurations::MaMiCoConfiguration<3>& mamicoConfig,
          const std::string filename,
			    const tarch::utils::MultiMDService<3>& multiMDService
    			){
    			return new coupling::services::MultiMDCellService<MY_LINKEDCELL,3>(
			        mdSolverInterfaces,
			        macroscopicSolverInterface, 
			        simpleMDConfig.getMPIConfiguration().getNumberOfProcesses(), 
			        rank, 
			        totalNumberMDSimulations,
			        mamicoConfig.getParticleInsertionConfiguration(), 
			        mamicoConfig.getMomentumInsertionConfiguration(), 
			        mamicoConfig.getBoundaryForceConfiguration(),
			        mamicoConfig.getTransferStrategyConfiguration(), 
			        mamicoConfig.getParallelTopologyConfiguration(), 
			        simpleMDConfig.getSimulationConfiguration().getNumberOfTimesteps(),
			        mamicoConfig.getMacroscopicCellConfiguration(), 
              filename.c_str(),
			        multiMDService
			    ); }
    		),
    		"mdSolverInterfaces"_a, "macroscopicSolverInterface"_a, "simpleMDConfig"_a, 
    		"rank"_a=0, "totalNumberMDSimulations"_a=1, "mamicoConfig"_a, "xmlConfigFilename"_a,
    		"multiMDService"_a)
    	.def("getMacroscopicCellService", &coupling::services::MultiMDCellService<MY_LINKEDCELL,3>::getMacroscopicCellService , py::return_value_policy::reference)
    	.def("sendFromMacro2MD", [](coupling::services::MultiMDCellService<MY_LINKEDCELL,3>& service, CouplingBuffer& buf){
            service.sendFromMacro2MD(buf._sendBuffer, buf._globalCellIndices4SendBuffer);
        })
    	.def("sendFromMD2Macro", [](coupling::services::MultiMDCellService<MY_LINKEDCELL,3>& service, CouplingBuffer& buf){
            service.sendFromMD2Macro(buf._recvBuffer, buf._globalCellIndices4RecvBuffer);
        } )
    ;

    py::class_<coupling::IndexConversion<3>>(coupling, "IndexConversion");

    py::class_<CouplingBuffer>(coupling, "Buffer")
        .def(py::init<const coupling::IndexConversion<3>&, 
                coupling::interface::MacroscopicSolverInterface<3> &,
                unsigned int,
                unsigned int
            >()
        )
        .def("recv2CSV", &CouplingBuffer::recv2CSV)
        .def("store2send", &CouplingBuffer::store2send)
        .def("loadRecvVelocity", &CouplingBuffer::loadRecvVelocity, py::return_value_policy::take_ownership)
        .def("loadRecvDensity", &CouplingBuffer::loadRecvDensity, py::return_value_policy::take_ownership)
    ;
}
