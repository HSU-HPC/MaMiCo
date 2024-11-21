#pragma once
#include "coupling/interface/MacroscopicSolverInterface.h"
#include "coupling/indexing/IndexTypes.h"

namespace coupling {
namespace preciceadapter {
  enum class DataType { scalar, vector };

  struct DataDescription {
    std::string name;
    DataType type;
  };

  template <unsigned int dim> class PreciceInterface;
} // namespace preciceadapter
} // namespace coupling

template <unsigned int dim> class coupling::preciceadapter::PreciceInterface : public coupling::interface::MacroscopicSolverInterface<dim> {
public:
  /**
   * Constructor
   */
  PreciceInterface(const bool twoWayCoupling) : _twoWayCoupling(twoWayCoupling) {}

  /**
   * Return true if the coupling cell with index idx is a macro to MD coupling cell
   */
  virtual bool isMacro2MD(I01 idx) const {
    return (I08::contains(idx) && !I12::contains(idx));
  }

  /**
   * Return true if the coupling cell with index idx is a MD to macro coupling cell
   */
  virtual bool isMD2Macro(I01 idx) const {
    return (I12::contains(idx));
  }

  /**
   * get the mesh name associated to this cell
   * this mesh is used to couple the macroscopic solver to the MD solver
   * @return mesh name
   */
  virtual std::string getMacro2MDSolverMeshName(I01 idx) const = 0;

  /**
   * get the mesh name associated to this cell
   * this mesh is used to couple the MD solver to the macroscopic solver
   * @return mesh name
   */
  virtual std::string getMD2MacroSolverMeshName(I01 idx) const = 0;

  /**
   * get the mesh off set associated to this cell
   * @return mesh offset
   */
  virtual tarch::la::Vector<dim, double> getMacro2MDSolverMeshOffset(I01 idx) const {
    return tarch::la::Vector<3, double>(0.0);
  }

  /**
   * get the mesh off set associated to this cell
   * @return mesh offset
   */
  virtual tarch::la::Vector<dim, double> getMD2MacroSolverMeshOffset(I01 idx) const {
    return tarch::la::Vector<3, double>(0.0);
  }

  /**
   * return true if two way coupling, i.e. macro to md and md to macro, is activated
   */
  bool twoWayCoupling() const {return _twoWayCoupling;} 

  /**
   * Return true if the coupling cell with index idx belongs to the mesh named meshName
   */
  virtual bool contains(std::string meshName, I01 idx) const = 0;

  /**
   * Add a data description description to this interface
   */
  void addDataDescription(const std::string& meshName, const DataDescription& dataDescription) {
    _dataDescriptions[meshName][dataDescription.name] = dataDescription;
  }

  /**
   * get a vector of data descriptions associated to this mesh
   * @return std::vector<DataDescription>
   */
  std::vector<DataDescription> getDataDescriptions(const std::string meshName) const {
    if (_dataDescriptions.find(meshName) == _dataDescriptions.end()) throw std::runtime_error("PreciceAdapter::getDataDescriptions: no mesh named " + meshName);
    std::vector<DataDescription> dataDescriptions;
    for (auto const& [dataName, dataDescription] : _dataDescriptions.at(meshName)) {
      dataDescriptions.push_back(dataDescription);
    }
    return dataDescriptions;
  }

  /**
   * get a data description associated to this mesh and data names
   * @return data
   */
  DataDescription getDataDescription(const std::string& meshName, const std::string& dataName) const {
    if (_dataDescriptions.find(meshName) == _dataDescriptions.end()) throw std::runtime_error("PreciceAdapter::getDataDescription: no mesh named " + meshName);
    if (_dataDescriptions.at(meshName).find(dataName) == _dataDescriptions.at(meshName).end()) throw std::runtime_error("PreciceAdapter::getDataDescription: no data " + dataName + " for mesh " + meshName);
    return _dataDescriptions.at(meshName).at(dataName);
  }

  /**
   * read the vector data sent from the macroscopic solver to precice and update the macroscopic cell
   */
  virtual void readVectorData(const std::string& meshName, const std::string& dataName, coupling::datastructures::CouplingCell<dim>* const cell, 
    const I01& idx, const double& vx, const double& vy, const double& vz) const {}

  /**
   * read the scalar data sent from the macroscopic solver to precice and update the macroscopic cell
   */
  virtual void readScalarData(const std::string& meshName, const std::string& dataName, coupling::datastructures::CouplingCell<dim>* const cell,
  const I01& idx, const double v) const {}

  /**
   * write the vector data sent from mamico and update the double values vx, vy and vz
   */
  virtual void writeVectorData(const std::string& meshName, const std::string& dataName, const coupling::datastructures::CouplingCell<dim>* const cell, 
  const I01& idx, double& vx, double& vy, double& vz) const {}

  /**
   * write the scalar data sent from mamico and update the double values vx, vy and vz
   */
  virtual void writeScalarData(const std::string& meshName, const std::string& dataName, const coupling::datastructures::CouplingCell<dim>* const cell,
  const I01& idx, double& v) const {}

private:
  std::map<std::string, std::map<std::string, DataDescription>> _dataDescriptions;
  const bool _twoWayCoupling;
};