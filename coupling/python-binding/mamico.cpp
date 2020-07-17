#include <pybind11/pybind11.h>
#include <mpi.h>
#include <string>

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
}