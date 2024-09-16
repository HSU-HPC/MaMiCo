// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#include "simplemd/molecule-mappings/Adios2Writer.h"
#include "math.h"
#include <cmath>

constexpr char const* rx_name = "rx";
constexpr char const* ry_name = "ry";
constexpr char const* rz_name = "rz";
constexpr char const* vx_name = "vx";
constexpr char const* vy_name = "vy";
constexpr char const* vz_name = "vz";
constexpr char const* v_abs_name = "v_abs";
constexpr char const* simtime_name = "simulationtime";
long unsigned int local;
long unsigned int global;
uint64_t offset;

int counter = 0;

simplemd::moleculemappings::Adios2Writer::Adios2Writer(const simplemd::services::ParallelTopologyService& parallelTopologyService,
                                                       const simplemd::services::MoleculeService& moleculeService,
                                                       const simplemd::configurations::MolecularDynamicsConfiguration& configuration
#if (MD_PARALLEL == MD_YES)
                                                       ,
                                                       MPI_Comm communicator
#endif
                                                       )
    : _parallelTopologyService(parallelTopologyService), _moleculeService(moleculeService), _timestep(0), _configuration(configuration)
#if (MD_PARALLEL == MD_YES)
      ,
      _communicator(communicator)
#endif
{
  initAdios2();
}

void simplemd::moleculemappings::Adios2Writer::initAdios2() {
#if (MD_PARALLEL == MD_YES)
  _inst = std::make_shared<adios2::ADIOS>(_communicator);
#else
  _inst = std::make_shared<adios2::ADIOS>();
#endif

  _numberProcesses = _configuration.getMPIConfiguration().getNumberOfProcesses()[0] * _configuration.getMPIConfiguration().getNumberOfProcesses()[1] *
                     _configuration.getMPIConfiguration().getNumberOfProcesses()[2];

  // Setup Adios2 IO
  _io = std::make_shared<adios2::IO>(_inst->DeclareIO("Output"));
  _io->SetEngine("BP4");
  _engine = std::make_shared<adios2::Engine>(_io->Open(_configuration.getAdios2Configuration().getFilename(), adios2::Mode::Write));
  offset = 0;

  _io->DefineAttribute<int>("num_processes", _numberProcesses);
  _io->DefineAttribute<int>("num_components", 1);

  // simpleMD only has 1 component
  std::string component_id = "component_0";
  std::vector<std::array<double, 3>> adios_centers;
  std::vector<double> adios_sigmas;
  std::vector<double> adios_mass;
  std::vector<double> adios_epsilon;
  std::array<double, 3> _domainCenter;

  for (int i = 0; i < MD_DIM; i++) {
    _domainSize[i] = (float)(_configuration.getDomainConfiguration().getGlobalDomainSize().operator[](i));
    _domainCenter[i] = (_configuration.getDomainConfiguration().getGlobalDomainSize().operator[](i) / 2);
    _offset[i] = (float)(_configuration.getDomainConfiguration().getGlobalDomainOffset().operator[](i));
  }

  adios_centers.emplace_back(_domainCenter);
  adios_sigmas.emplace_back(_configuration.getMoleculeConfiguration().getSigma());
  adios_mass.emplace_back(_configuration.getMoleculeConfiguration().getMass());
  adios_epsilon.emplace_back(_configuration.getMoleculeConfiguration().getEpsilon());

  _io->DefineAttribute<double>(component_id + "_centers", adios_centers[0].data(), adios_centers.size() * 3);
  _io->DefineAttribute<double>(component_id + "_sigma", adios_sigmas.data(), adios_sigmas.size());
  _io->DefineAttribute<double>(component_id + "_mass", adios_mass.data(), adios_mass.size());
  _io->DefineAttribute<double>(component_id + "_epsilon", adios_epsilon.data(), adios_epsilon.size());
  _io->DefineAttribute<std::string>(component_id + "_name", "component_1");
  _io->DefineAttribute<std::string>(component_id + "_element_names", "SimpleMD_Particle");
}

simplemd::moleculemappings::Adios2Writer::~Adios2Writer() { _engine->Close(); }

void simplemd::moleculemappings::Adios2Writer::setTimestep(const unsigned int& timestep) { _timestep = timestep; }

void simplemd::moleculemappings::Adios2Writer::beginMoleculeIteration() {

  _engine->BeginStep();
  _io->RemoveAllVariables(); // empty buffer

  local = _moleculeService.getNumberMolecules();
  global = local; // valid for sequential execution, overwritten in parallel case

#if (MD_PARALLEL == MD_YES)
  MPI_Exscan(&local, &offset, 1, MPI_UINT64_T, MPI_SUM, _communicator);
  MPI_Allreduce(&local, &global, 1, MPI_UINT64_T, MPI_SUM, _communicator);
#endif

  // define velocity- and position-variables in Adios2 IO
  rx_var = _io->DefineVariable<float>(rx_name, {global}, {offset}, {local}, adios2::ConstantDims);
  ry_var = _io->DefineVariable<float>(ry_name, {global}, {offset}, {local}, adios2::ConstantDims);
  rz_var = _io->DefineVariable<float>(rz_name, {global}, {offset}, {local}, adios2::ConstantDims);
  vx_var = _io->DefineVariable<float>(vx_name, {global}, {offset}, {local}, adios2::ConstantDims);
  vy_var = _io->DefineVariable<float>(vy_name, {global}, {offset}, {local}, adios2::ConstantDims);
  vz_var = _io->DefineVariable<float>(vz_name, {global}, {offset}, {local}, adios2::ConstantDims);
  v_abs_var = _io->DefineVariable<float>(v_abs_name, {global}, {offset}, {local}, adios2::ConstantDims);
  time_var = _io->DefineVariable<double>(simtime_name);
  global_box_var = _io->DefineVariable<float>("global_box", {6}, {0}, {6}, adios2::ConstantDims);
  component_id_var = _io->DefineVariable<uint64_t>("component_id", {global}, {offset}, {local}, adios2::ConstantDims);
}

void simplemd::moleculemappings::Adios2Writer::handleMolecule(Molecule& molecule) {

  _positionsx.push_back((float)molecule.getConstPosition()[0]);
  _velocitiesx.push_back((float)molecule.getConstVelocity()[0]);
  _positionsy.push_back((float)molecule.getConstPosition()[1]);
  _velocitiesy.push_back((float)molecule.getConstVelocity()[1]);
  _positionsz.push_back((float)molecule.getConstPosition()[2]);
  _velocitiesz.push_back((float)molecule.getConstVelocity()[2]);
  _velocities_abs.push_back(
      (float)sqrt(pow(molecule.getConstVelocity()[0], 2) + pow(molecule.getConstVelocity()[1], 2) + pow(molecule.getConstVelocity()[2], 2)));
  _component_id.push_back((uint64_t)0);
  counter++;
}

void simplemd::moleculemappings::Adios2Writer::endMoleculeIteration() {

  std::array<float, 6> global_box = {_offset[0], _offset[1], _offset[2], _domainSize[0] + _offset[0], _domainSize[1] + _offset[1], _domainSize[2] + _offset[2]};

  // write data
  _engine->Put<float>(rx_var, _positionsx.data());
  _engine->Put<float>(ry_var, _positionsy.data());
  _engine->Put<float>(rz_var, _positionsz.data());
  _engine->Put<float>(vx_var, _velocitiesx.data());
  _engine->Put<float>(vy_var, _velocitiesy.data());
  _engine->Put<float>(vz_var, _velocitiesz.data());
  _engine->Put<float>(v_abs_var, _velocities_abs.data());
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
  _velocities_abs.clear();
  _component_id.clear();
}