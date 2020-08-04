#include <pybind11/pybind11.h>
#include <mpi.h>
#include <string>
#include <limits>

#include "tarch/utils/MultiMDService.h"
#include "coupling/CouplingMDDefinitions.h"
#include "tarch/configuration/ParseConfiguration.h"
#include "simplemd/configurations/MolecularDynamicsConfiguration.h"
#include "coupling/configurations/MaMiCoConfiguration.h"
#include "coupling/interface/MDSimulationFactory.h"

// for debugging purposes only
#include <iostream>

namespace py = pybind11;
using namespace pybind11::literals;

int initMPI(){
	py::list argv_list = py::module::import("sys").attr("argv");

	char** argv = new char * [argv_list.size()];
	int argc = 0;
	for (auto arg : argv_list)
		argv[argc++] = (new std::string(PyUnicode_AsUTF8(arg.ptr())))->data();  //todo mini memleak here

	// DEBUG
	// std::cout << "argc = " << argc << std::endl;
	// for(int i = 0; i < argc; i++)
	// std::cout << "argv[" << i << "] = " << argv[i] << std::endl;

	int rank = 0;
	MPI_Init(&argc,&argv);
	delete[] argv;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	return rank;
}

template<class T> T* makeConfiguration(const std::string filename, const std::string topleveltag){
	T* cfg = new T{};
	tarch::configuration::ParseConfiguration::parseConfiguration<T>(filename,topleveltag,*cfg);
	return cfg;
}

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

PYBIND11_MODULE(mamico, mamico) {
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

	py::class_<coupling::configurations::MaMiCoConfiguration<3>>(configuration, "MaMiCoConfiguration")
		.def("__repr__",
	        [](const coupling::configurations::MaMiCoConfiguration<3> &c) {
	            return "<MaMiCoConfiguration read from XML tag \"" + c.getTag() + "\">";
	        }
	    )
	    .def("isValid", &coupling::configurations::MaMiCoConfiguration<3>::isValid);

	configuration.def("parseMolecularDynamicsConfiguration", &makeConfiguration<simplemd::configurations::MolecularDynamicsConfiguration>, py::return_value_policy::take_ownership);
	configuration.def("parseMaMiCoConfiguration", &makeConfiguration<coupling::configurations::MaMiCoConfiguration<3>>, py::return_value_policy::take_ownership);

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
    	; 

    coupling.def("getMDSimulation", []
    	(const simplemd::configurations::MolecularDynamicsConfiguration& c1,
    	 const coupling::configurations::MaMiCoConfiguration<3>& c2,
    	 void* comm){
    		return coupling::interface::SimulationAndInterfaceFactory::getInstance().getMDSimulation(c1,c2,(MPI_Comm) comm);
    	}, "simpleMDConfig"_a, "mamicoConfig"_a, "localComm"_a, py::return_value_policy::take_ownership);
}
