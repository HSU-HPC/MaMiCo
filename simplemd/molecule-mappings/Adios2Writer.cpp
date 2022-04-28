// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#if BUILD_WITH_ADIOS2

#include "simplemd/molecule-mappings/Adios2Writer.h"

constexpr char const* rx_name = "rx";
constexpr char const* ry_name = "ry";
constexpr char const* rz_name = "rz";
constexpr char const* vx_name = "vx";
constexpr char const* vy_name = "vy";
constexpr char const* vz_name = "vz";
constexpr char const* simtime_name = "simulationtime";
long unsigned int local;
long unsigned int global;
uint64_t offset;

int counter = 0;

simplemd::moleculemappings::Adios2Writer::Adios2Writer(
const simplemd::services::ParallelTopologyService& parallelTopologyService,
const simplemd::services::MoleculeService& moleculeService,
const std::string& filename,
tarch::la::Vector<MD_DIM,double> offset,
const simplemd::configurations::MolecularDynamicsConfiguration &configuration
#if (MD_PARALLEL==MD_YES)
,MPI_Comm communicator
#endif
):
_parallelTopologyService(parallelTopologyService),_moleculeService(moleculeService),
_filename(filename), _timestep(0), _offset(offset),_configuration(configuration)
{
  //std::cout<<"constructor"<<std::endl;
  #if (MD_PARALLEL==MD_YES)
  initAdios2(communicator);
  #else
  initAdios2();
  #endif
}


void simplemd::moleculemappings::Adios2Writer::initAdios2(
  #if (MD_PARALLEL==MD_YES)
  MPI_Comm communicator
  #endif
  ) {
  #if (MD_PARALLEL==MD_YES)
  _inst = std::make_shared<adios2::ADIOS>(communicator);
  #else
  _inst = std::make_shared<adios2::ADIOS>();
  #endif

  _numberProcesses = _configuration.getMPIConfiguration().getNumberOfProcesses()[0]*_configuration.getMPIConfiguration().getNumberOfProcesses()[1]*_configuration.getMPIConfiguration().getNumberOfProcesses()[2];

  _io = std::make_shared<adios2::IO>(_inst->DeclareIO("Output"));
  _io->SetEngine("BP4");
  _engine = std::make_shared<adios2::Engine>(_io->Open("adios2output.bp", adios2::Mode::Write));

  offset = 0;


  _io->DefineAttribute<int>("num_processes", _numberProcesses);
  _io->DefineAttribute<int>("num_components", 1);

  std::string component_id = "component_0";
  std::vector<std::array<double, 3>> adios_centers;
  std::array<double, 3> _domainCenter;
  std::vector<double> adios_sigmas;
  std::vector<double> adios_mass;
  std::vector<double> adios_epsilon;
  
    _domainCenter[0] = (_configuration.getDomainConfiguration().getGlobalDomainSize().operator[](0)/2);
    _domainCenter[1] = (_configuration.getDomainConfiguration().getGlobalDomainSize().operator[](1)/2);
    _domainCenter[2] = (_configuration.getDomainConfiguration().getGlobalDomainSize().operator[](2)/2);
    adios_centers.emplace_back(_domainCenter);
    adios_sigmas.emplace_back(_configuration.getMoleculeConfiguration().getSigma());
    adios_mass.emplace_back(_configuration.getMoleculeConfiguration().getMass());
    adios_epsilon.emplace_back(_configuration.getMoleculeConfiguration().getEpsilon());
  

  _io->DefineAttribute<double>(component_id + "_centers", adios_centers[0].data(),adios_centers.size() * 3);
  _io->DefineAttribute<double>(component_id + "_sigma", adios_sigmas.data(), adios_sigmas.size());
  _io->DefineAttribute<double>(component_id + "_mass", adios_mass.data(), adios_mass.size());
  _io->DefineAttribute<double>(component_id + "_epsilon", adios_epsilon.data(), adios_epsilon.size());
  _io->DefineAttribute<std::string>(component_id + "_name", "component_1");
  _io->DefineAttribute<std::string>(component_id + "_element_names", "SimpleMD_Particle");

  if(rx_var)
    std::cout<<"valid"<<std::endl;

  //int _writefrequency = 50000;
  //_append_mode = "OFF";
  //_compression = "none";
  //_compression_accuracy = "0.00001";
  //_compression_rate = "8";
  //_num_files = -1;*/

}

simplemd::moleculemappings::Adios2Writer::~Adios2Writer(){
 _engine->Close();
}



void simplemd::moleculemappings::Adios2Writer::setTimestep(const unsigned int &timestep){
  _timestep = timestep;
}


void simplemd::moleculemappings::Adios2Writer::beginMoleculeIteration(){


}


void simplemd::moleculemappings::Adios2Writer::handleMolecule(Molecule &molecule){
  
    _positionsx.push_back((float)molecule.getConstPosition()[0]);
    _velocitiesx.push_back((float)molecule.getConstVelocity()[0]);
    _positionsy.push_back((float)molecule.getConstPosition()[1]);
    _velocitiesy.push_back((float)molecule.getConstVelocity()[1]);
    _positionsz.push_back((float)molecule.getConstPosition()[2]);
    _velocitiesz.push_back((float)molecule.getConstVelocity()[2]);
    _component_id.push_back((uint64_t)0);
    counter++;

}



void simplemd::moleculemappings::Adios2Writer::endMoleculeIteration(){
  _engine->BeginStep();
  _io->RemoveAllVariables();
  
  local = _moleculeService.getNumberMolecules();
  global = local;
  #if (MD_PARALLEL == MD_YES)
  MPI_Exscan(&local, &offset, 1, MPI_UINT64_T, MPI_SUM, communicator);
  MPI_Allreduce(&local, &global, 1, MPI_UINT64_T, MPI_SUM, communicator);
  #endif

  rx_var = _io->DefineVariable<float>(rx_name, {global}, {offset}, {local}, adios2::ConstantDims);
  ry_var = _io->DefineVariable<float>(ry_name, {global}, {offset}, {local}, adios2::ConstantDims);
  rz_var = _io->DefineVariable<float>(rz_name, {global}, {offset}, {local}, adios2::ConstantDims);
  vx_var = _io->DefineVariable<float>(vx_name, {global}, {offset}, {local}, adios2::ConstantDims);
  vy_var = _io->DefineVariable<float>(vy_name, {global}, {offset}, {local}, adios2::ConstantDims);
  vz_var = _io->DefineVariable<float>(vz_name, {global}, {offset}, {local}, adios2::ConstantDims);
  time_var = _io->DefineVariable<double>(simtime_name);
  global_box_var = _io->DefineVariable<float>("global_box", {6}, {0}, {6}, adios2::ConstantDims);
  component_id_var = _io->DefineVariable<uint64_t>("component_id", {global}, {offset}, {local},adios2::ConstantDims);
  //component_var = io->DefineVariable<uint64_t>

  std::array<double, 6> tmp_global_box = {_offset[0], _offset[1], _offset[2], _domainCenter[0]*2, _domainCenter[1]*2, _domainCenter[2]*2};
	std::array<float, 6> global_box;
	std::copy(tmp_global_box.begin(), tmp_global_box.end(), global_box.begin());

  /*const auto rxVar = _io->InquireVariable<float>(rx_name);
  const auto ryVar = _io->InquireVariable<float>(ry_name);
  const auto rzVar = _io->InquireVariable<float>(rz_name);
  const auto vxVar = _io->InquireVariable<float>(vx_name);
  const auto vyVar = _io->InquireVariable<float>(vy_name);
  const auto vzVar = _io->InquireVariable<float>(vz_name);
  const auto timeVar = _io->InquireVariable<double>(simtime_name);*/
  //std::cout<<_positionsx.size()<<std::endl;

  _engine->Put<float>(rx_var, _positionsx.data());
  _engine->Put<float>(ry_var, _positionsy.data());
  _engine->Put<float>(rz_var, _positionsz.data());
  _engine->Put<float>(vx_var, _velocitiesx.data());
  _engine->Put<float>(vy_var, _velocitiesy.data());
  _engine->Put<float>(vz_var, _velocitiesz.data());
  _engine->Put<float>(global_box_var, global_box.data());
  _engine->Put<uint64_t>(component_id_var, _component_id.data());
  _engine->Put<double>(time_var, _timestep);

  _engine->EndStep();

  _positionsx.clear();
  _positionsy.clear();
  _positionsz.clear();
  _velocitiesx.clear();
  _velocitiesy.clear();
  _velocitiesz.clear();
  _component_id.clear();
}

#endif