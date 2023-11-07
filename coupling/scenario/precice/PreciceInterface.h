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
   * get the mesh name associated to this cell
   * @return mesh name
   */
  virtual std::string getMeshName(tarch::la::Vector<dim, unsigned int> globalCellIndex) = 0;

  /**
   * get the mesh off set associated to this cell
   * @return mesh offset
   */
  virtual tarch::la::Vector<dim, double> getMeshOffset(tarch::la::Vector<dim, unsigned int> globalCellIndex) = 0;

  /**
   * get a vector of data associated to this mesh
   * @return std::vector<Data>
   */
  virtual std::vector<Data> getData(std::string meshName) = 0;

  /**
   * get a data object associated to this mesh and data names
   * @return data
   */
  virtual Data getData(std::string meshName, std::string dataName) = 0;

  /**
   * read the vector data sent from the macroscopic solver to precice and update the macroscopic cell
   */
  virtual void readVectorData(std::string meshName, std::string dataName, coupling::datastructures::MacroscopicCell<dim>* const cell, const double vx, const double vy, const double vz) = 0;

  /**
   * read the scalar data sent from the macroscopic solver to precice and update the macroscopic cell
   */
  virtual void readScalarData(std::string meshName, std::string dataName, coupling::datastructures::MacroscopicCell<dim>* const cell, const double v) = 0;

  /**
   * write the vector data sent from mamico and update the double values vx, vy and vz
   */
  virtual void writeVectorData(std::string meshName, std::string dataName, const coupling::datastructures::MacroscopicCell<dim>* const cell, double& vx, double& vy, double& vz) = 0;

  /**
   * write the scalar data sent from mamico and update the double values vx, vy and vz
   */
  virtual void writeScalarData(std::string meshName, std::string dataName, const coupling::datastructures::MacroscopicCell<dim>* const cell, double& v) = 0;
};