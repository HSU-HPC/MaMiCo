#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

#include "coupling/interface/impl/ls1/LS1RegionWrapper.h"

#include "coupling/interface/impl/ls1/LS1Molecule.h"

#define MY_LINKEDCELL ls1::LS1RegionWrapper

class LS1RegionWrapperTest : public CppUnit::TestFixture
{
	CPPUNIT_TEST_SUITE(LS1RegionWrapperTest);
	CPPUNIT_TEST(testPotentialEnergyAtPoint);
	CPPUNIT_TEST_SUITE_END();

public:
	void setUp()
	{
		global_log = new Log::Logger(Log::Error);
		coupling::interface::LS1StaticCommData::getInstance().setBoxOffsetAtDim(0,0);
		coupling::interface::LS1StaticCommData::getInstance().setBoxOffsetAtDim(1,0);
		coupling::interface::LS1StaticCommData::getInstance().setBoxOffsetAtDim(2,0);
		coupling::interface::LS1StaticCommData::getInstance().setLocalCommunicator(MPI_COMM_WORLD);
		_testSimulation = new Simulation();
		global_simulation = _testSimulation;
		_testSimulation->disableFinalCheckpoint();
		_testSimulation->readConfigFile("../test/unit/coupling/interface/impl/ls1/ls1gridconfig.xml");
		//_testSimulation->getDomain()->thermostatOff();
		//_testSimulation->getDomain()->setExplosionHeuristics(false);
		// after this point the mamico plugin exists and is accessible
		_testSimulation->prepare_start();
		_testSimulation->preSimLoopSteps();
	}
	void tearDown()
	{
		if(_testSimulation != nullptr)
		{
			delete _testSimulation;
			_testSimulation = nullptr;
		}
	}
	void testPotentialEnergyAtPoint()
	{
		std::cout << "Particle density: " << _testSimulation->getTotalNumberOfMolecules() << 
			" Resultant density: " << static_cast<double>(_testSimulation->getTotalNumberOfMolecules())/(10*10*10) << std::endl;
		for(int i = 4; i < 8; i++)
		{
			for(int j = 4; j < 8; j++)
			{
				for(int k = 4; k < 8; k++)
				{
					tarch::la::Vector<3,double> position{i+0.5,j+0.5,k+0.5};
					ls1::LS1RegionWrapper fullDomainWrapper(_testSimulation->getEnsemble()->domain()->rmin(), _testSimulation->getEnsemble()->domain()->rmax(), _testSimulation);
					tarch::la::Vector<3,double> force (0.0);
					double pseudoPotential = 0.0;
					std::tie(force, pseudoPotential) = fullDomainWrapper.calculateForceAndPotentialAtPoint(position, true);

					::Molecule ls1Molecule;
					coupling::interface::LS1Molecule mamicoMolecule(&ls1Molecule);
					mamicoMolecule.setPosition(position);
					double realPotential = mamicoMolecule.getPotentialEnergy();
					CPPUNIT_ASSERT_MESSAGE( "position assersion",mamicoMolecule.getPosition() == position );
					CPPUNIT_ASSERT_MESSAGE( "potential assertion",pseudoPotential == realPotential );
				}
			}
			
		}
	}
private:
	Simulation* _testSimulation;
};

CPPUNIT_TEST_SUITE_REGISTRATION(LS1RegionWrapperTest);