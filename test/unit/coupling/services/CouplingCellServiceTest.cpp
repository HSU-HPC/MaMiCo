#include "coupling/CouplingMDDefinitions.h"
#include "coupling/services/CouplingCellService.h"
#include "coupling/configurations/BoundaryForceConfiguration.h"
#include "coupling/configurations/CouplingCellConfiguration.h"
#include "coupling/configurations/MomentumInsertionConfiguration.h"
#include "coupling/configurations/ParallelTopologyConfiguration.h"
#include "coupling/configurations/ParticleInsertionConfiguration.h"
#include "coupling/configurations/ThermostatConfiguration.h"
#include "coupling/configurations/TransferStrategyConfiguration.h"
#include "coupling/interface/MDSolverInterface.h"
#include "coupling/interface/MacroscopicSolverInterface.h"
#include "coupling/interface/impl/SimpleMD/SimpleMDLinkedCellWrapper.h"
#include "tarch/utils/MultiMDService.h"
#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

using namespace coupling::indexing;

/**
 *  @author Louis Viot
 */
class CouplingCellServiceTest : public CppUnit::TestFixture {

  CPPUNIT_TEST_SUITE(CouplingCellServiceTest);
  CPPUNIT_TEST(test<3>);
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp() {}

  void tearDown() {}

  template <unsigned int dim> void test() {
    int rank = 0;
    const tarch::la::Vector<dim, unsigned int> numberOfProcesses{1, 1, 1};

    IndexingService<3>::getInstance().initWithCells({1}, numberOfProcesses, coupling::paralleltopology::XYZ, 1, rank);

    coupling::interface::MDSolverInterface<coupling::interface::SimpleMDLinkedCellWrapper, dim>* mdSolverInterface = new TestMDSolverInterface<dim>();
    coupling::interface::MacroscopicSolverInterface<dim>* macroscopicSolverInterface = new TestMacroscopicSolverInterface<dim>();

    TestParticleInsertionConfiguration particleInsertionConfiguration;
    TestMomentumInsertionConfiguration<dim> momentumInsertionConfiguration;
    TestBoundaryForceConfiguration<dim> boundaryForceConfiguration;
    TestTransferStrategyConfiguration<dim> transferStrategyConfiguration;
    TestParallelTopologyConfiguration parallelTopologyConfiguration;
    TestThermostatConfiguration thermostatConfiguration;
    const unsigned int numberOfTimeSteps = 1;
    TestCouplingCellConfiguration<dim> couplingCellConfiguration;
    const char* filterPipelineConfigurationFile = nullptr;
    tarch::utils::MultiMDService<dim> multiMDService{tarch::la::Vector<dim, unsigned int>{1}, 1};
    const unsigned int topologyOffset = 1;
    new coupling::services::CouplingCellServiceImpl<coupling::interface::SimpleMDLinkedCellWrapper, dim>(
        1, mdSolverInterface, macroscopicSolverInterface, numberOfProcesses, rank, particleInsertionConfiguration, momentumInsertionConfiguration,
        boundaryForceConfiguration, transferStrategyConfiguration, parallelTopologyConfiguration, thermostatConfiguration, numberOfTimeSteps,
        couplingCellConfiguration, filterPipelineConfigurationFile, multiMDService, topologyOffset);
  }

private:
  template <unsigned int dim> class TestMDSolverInterface : public coupling::interface::MDSolverInterface<coupling::interface::SimpleMDLinkedCellWrapper, dim> {
  public:
    TestMDSolverInterface() : coupling::interface::MDSolverInterface<coupling::interface::SimpleMDLinkedCellWrapper, dim>() {}
    virtual ~TestMDSolverInterface() {}

    coupling::interface::SimpleMDLinkedCellWrapper&
    getLinkedCell(const typename coupling::interface::MDSolverInterface<coupling::interface::SimpleMDLinkedCellWrapper, dim>::CellIndex_T& couplingCellIndex,
                  const tarch::la::Vector<dim, unsigned int>& linkedCellInCouplingCell,
                  const tarch::la::Vector<dim, unsigned int>& linkedCellsPerCouplingCell) override {
      return _linkedCellWrapper;
    }

    virtual tarch::la::Vector<dim, double> getGlobalMDDomainSize() const override { return tarch::la::Vector<dim, double>(1.0); }

    virtual tarch::la::Vector<dim, double> getGlobalMDDomainOffset() const override { return tarch::la::Vector<dim, double>(0.0); }

    virtual double getMoleculeMass() const override { return 1.0; }
    virtual double getKB() const override { return 1.0; }
    virtual double getMoleculeSigma() const override { return 1.0; }
    virtual double getMoleculeEpsilon() const override { return 1.0; }
    virtual void getInitialVelocity(const tarch::la::Vector<dim, double>& meanVelocity, const double& kB, const double& temperature,
                                    tarch::la::Vector<dim, double>& initialVelocity) const override {}
    virtual void deleteMoleculeFromMDSimulation(const coupling::interface::Molecule<dim>& molecule,
                                                coupling::interface::SimpleMDLinkedCellWrapper& cell) override {}
    virtual void addMoleculeToMDSimulation(const coupling::interface::Molecule<dim>& molecule) override {}
    virtual void setupPotentialEnergyLandscape(const tarch::la::Vector<dim, unsigned int>& indexOfFirstCouplingCell,
                                               const tarch::la::Vector<dim, unsigned int>& rangeCouplingCells,
                                               const tarch::la::Vector<dim, unsigned int>& linkedCellsPerCouplingCell) override {}
    virtual tarch::la::Vector<dim, unsigned int> getLinkedCellIndexForMoleculePosition(const tarch::la::Vector<dim, double>& position) override {
      return tarch::la::Vector<dim, unsigned int>(0);
    }
    virtual void calculateForceAndEnergy(coupling::interface::Molecule<dim>& molecule) override {}
    virtual void synchronizeMoleculesAfterMassModification() override {}
    virtual void synchronizeMoleculesAfterMomentumModification() override {}
    virtual double getDt() override { return 1.0; }
    virtual coupling::interface::MoleculeIterator<coupling::interface::SimpleMDLinkedCellWrapper, dim>*
    getMoleculeIterator(coupling::interface::SimpleMDLinkedCellWrapper& cell) override {
      return NULL;
    }
    coupling::interface::SimpleMDLinkedCellWrapper _linkedCellWrapper;
  };

  template <unsigned int dim> class TestMacroscopicSolverInterface : public coupling::interface::MacroscopicSolverInterface<dim> {
  public:
    ~TestMacroscopicSolverInterface() {}
    std::vector<unsigned int> getRanks(I01 globalCellIndex) override { return {1}; };
    unsigned int getOuterRegion() override { return 1; }
  };

  class TestParticleInsertionConfiguration : public coupling::configurations::ParticleInsertionConfiguration {
  public:
    TestParticleInsertionConfiguration()
        : coupling::configurations::ParticleInsertionConfiguration(
              1, 1, 1, 1, 1, 1, 1, 1, 1, coupling::configurations::ParticleInsertionConfiguration::ParticleInsertionType::NO_INSERTION) {}

    ~TestParticleInsertionConfiguration() {}

    void parseSubtag(tinyxml2::XMLElement* node) override {};

    std::string getTag() const override { return ""; };

    bool isValid() const override { return true; };
  };

  template <unsigned int dim> class TestMomentumInsertionConfiguration : public coupling::configurations::MomentumInsertionConfiguration<dim> {
  public:
    TestMomentumInsertionConfiguration()
        : coupling::configurations::MomentumInsertionConfiguration<dim>(
              coupling::configurations::MomentumInsertionConfiguration<dim>::MomentumInsertionType::NO_INSERTION) {}

    ~TestMomentumInsertionConfiguration() {}

    void parseSubtag(tinyxml2::XMLElement* node) override {};

    std::string getTag() const override { return ""; };

    bool isValid() const override { return true; };
  };

  template <unsigned int dim> class TestBoundaryForceConfiguration : public coupling::configurations::BoundaryForceConfiguration<dim> {
  public:
    TestBoundaryForceConfiguration()
        : coupling::configurations::BoundaryForceConfiguration<dim>(
              coupling::configurations::BoundaryForceConfiguration<dim>::BoundaryForceType::NO_BOUNDARYFORCE, 1, 1, {true}) {}

    ~TestBoundaryForceConfiguration() {}

    void parseSubtag(tinyxml2::XMLElement* node) override {};

    std::string getTag() const override { return ""; };

    bool isValid() const override { return true; };
  };

  template <unsigned int dim> class TestTransferStrategyConfiguration : public coupling::configurations::TransferStrategyConfiguration<dim> {
  public:
    TestTransferStrategyConfiguration()
        : coupling::configurations::TransferStrategyConfiguration<dim>(
              coupling::configurations::TransferStrategyConfiguration<dim>::StrategyType::DirectTransferStrategy, {true}, 1) {}

    ~TestTransferStrategyConfiguration() {}

    void parseSubtag(tinyxml2::XMLElement* node) override {};

    std::string getTag() const override { return ""; };

    bool isValid() const override { return true; };
  };

  class TestParallelTopologyConfiguration : public coupling::configurations::ParallelTopologyConfiguration {
  public:
    TestParallelTopologyConfiguration() : coupling::configurations::ParallelTopologyConfiguration(coupling::paralleltopology::XYZ) {}

    ~TestParallelTopologyConfiguration() {}

    void parseSubtag(tinyxml2::XMLElement* node) override {};

    std::string getTag() const override { return ""; };

    bool isValid() const override { return true; };
  };

  class TestThermostatConfiguration : public coupling::configurations::ThermostatConfiguration {
  public:
    TestThermostatConfiguration()
        : coupling::configurations::ThermostatConfiguration(coupling::configurations::ThermostatConfiguration::ThermostatRegion::none) {}

    ~TestThermostatConfiguration() {}

    void parseSubtag(tinyxml2::XMLElement* node) override {};

    std::string getTag() const override { return ""; };

    bool isValid() const override { return true; };
  };

  template <unsigned int dim> class TestCouplingCellConfiguration : public coupling::configurations::CouplingCellConfiguration<dim> {
  public:
    TestCouplingCellConfiguration() : coupling::configurations::CouplingCellConfiguration<dim>({1}, {1}, 1, "", 1, "") {}

    ~TestCouplingCellConfiguration() {}

    void parseSubtag(tinyxml2::XMLElement* node) override {};

    std::string getTag() const override { return ""; };

    bool isValid() const override { return true; };
  };
};

CPPUNIT_TEST_SUITE_REGISTRATION(CouplingCellServiceTest);