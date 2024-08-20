#pragma once
#include "coupling/interface/MacroscopicSolverInterface.h"

namespace coupling {
namespace preciceadapter {
  enum class DataType { scalar, vector };

  struct Data {
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
  PreciceInterface(bool twoWayCoupling) : _twoWayCoupling(twoWayCoupling) {}

  /**
   * get the mesh name associated to this cell
   * this mesh is used to couple the macroscopic solver to the MD solver
   * @return mesh name
   */
  virtual std::string getMacroscopicToMDSolverMeshName(tarch::la::Vector<dim, unsigned int> globalCellIndex) = 0;

  /**
   * get the mesh name associated to this cell
   * this mesh is used to couple the MD solver to the macroscopic solver
   * @return mesh name
   */
  virtual std::string getMDToMacroscopicSolverMeshName(tarch::la::Vector<dim, unsigned int> globalCellIndex) = 0;

  /**
   * get the mesh off set associated to this cell
   * @return mesh offset
   */
  virtual tarch::la::Vector<dim, double> getMacroscopicToMDSolverMeshOffset(tarch::la::Vector<dim, unsigned int> globalCellIndex) {
    return tarch::la::Vector<3, double>(0.0);
  }

  /**
   * get the mesh off set associated to this cell
   * @return mesh offset
   */
  virtual tarch::la::Vector<dim, double> getMDToMacroscopicSolverMeshOffset(tarch::la::Vector<dim, unsigned int> globalCellIndex) {
    return tarch::la::Vector<3, double>(0.0);
  }

  /**
   * return true if two way coupling, i.e. macro to md and md to macro, is activated
   */
  bool twoWayCoupling() {return _twoWayCoupling;} 

  /**
   * Add data to this interface
   */
  void addData(std::string meshName, Data data) {
    _data[meshName][data.name] = data;
  }

  /**
   * get a vector of data associated to this mesh
   * @return std::vector<Data>
   */
  std::vector<Data> getData(std::string meshName) {
    if (_data.find(meshName) == _data.end()) throw std::runtime_error("PreciceAdapter::getData: no mesh named " + meshName);
    std::vector<Data> data;
    for (std::map<std::string, Data>::iterator it = _data[meshName].begin(); it != _data[meshName].end(); ++it) {
      data.push_back(it->second);
    }
    return data;
  }

  /**
   * get a data object associated to this mesh and data names
   * @return data
   */
  Data getData(std::string meshName, std::string dataName) {
    if (_data.find(meshName) == _data.end()) throw std::runtime_error("PreciceAdapter::getData: no mesh named " + meshName);
    if (_data[meshName].find(dataName) == _data[meshName].end()) throw std::runtime_error("PreciceAdapter::getData: no data " + dataName + " for mesh " + meshName);
    return _data[meshName][dataName];
  }

  /**
   * read the vector data sent from the macroscopic solver to precice and update the macroscopic cell
   */
  virtual void readVectorData(std::string meshName, std::string dataName, coupling::datastructures::MacroscopicCell<dim>* const cell, const double vx, const double vy, const double vz) {}

  /**
   * read the scalar data sent from the macroscopic solver to precice and update the macroscopic cell
   */
  virtual void readScalarData(std::string meshName, std::string dataName, coupling::datastructures::MacroscopicCell<dim>* const cell, const double v) {}

  /**
   * write the vector data sent from mamico and update the double values vx, vy and vz
   */
  virtual void writeVectorData(std::string meshName, std::string dataName, const coupling::datastructures::MacroscopicCell<dim>* const cell, double& vx, double& vy, double& vz) {}

  /**
   * write the scalar data sent from mamico and update the double values vx, vy and vz
   */
  virtual void writeScalarData(std::string meshName, std::string dataName, const coupling::datastructures::MacroscopicCell<dim>* const cell, double& v) {}

private:
  std::map<std::string, std::map<std::string, Data>> _data;
  bool _twoWayCoupling;
};