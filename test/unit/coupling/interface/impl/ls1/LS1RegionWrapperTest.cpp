#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

#include "coupling/interface/impl/ls1/LS1RegionWrapper.h"

#include "coupling/interface/impl/ls1/LS1Molecule.h"

#define MY_LINKEDCELL ls1::LS1RegionWrapper

class LS1RegionWrapperTest : public CppUnit::TestFixture
{
	CPPUNIT_TEST_SUITE(LS1RegionWrapperTest);
	CPPUNIT_TEST(testPointIsInRegion);
	CPPUNIT_TEST(testSetRegion);
	CPPUNIT_TEST(testPotentialEnergyAtPoint);
	CPPUNIT_TEST(testAddAndDeleteParticle);
	CPPUNIT_TEST(testResetIterator);
	CPPUNIT_TEST(testIteration);
	CPPUNIT_TEST_SUITE_END();

public:
	void setUp()
	{
		global_log = new Log::Logger(Log::None);
		global_log->set_mpi_output_root(0);
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
	void testAddAndDeleteParticle()
	{
		//choose a position
		tarch::la::Vector<3,double> position = {2,2,2};

		//init
		double bBoxMin[3];
		double bBoxMax[3];
		global_simulation->domainDecomposition().getBoundingBoxMinMax(global_simulation->getDomain(), bBoxMin, bBoxMax);
		ls1::LS1RegionWrapper wrapper(bBoxMin, bBoxMax, _testSimulation);

		//verify that position is empty
		bool found = false;
		while(wrapper.iteratorValid())
		{
			::Molecule* temp = wrapper.getParticleAtIterator();
			if(temp->r(0) == position[0] && temp->r(1) == position[1] && temp->r(2) == position[2])
			{
				found = true;
				break;
			}
			wrapper.iteratorNext();
		}
		CPPUNIT_ASSERT_MESSAGE("Found before insertion", !found);

		//insert particle at position
		wrapper.setupIDcounterForParticleAddition();
		::Molecule blank;
		coupling::interface::LS1Molecule tempParticle(&blank);
		tempParticle.setPosition(position);
		wrapper.addMolecule(tempParticle);

		//verify that position is filled
		found = false; //sanity
		wrapper.iteratorReset();
		while(wrapper.iteratorValid())
		{
			::Molecule* temp = wrapper.getParticleAtIterator();
			if(temp->r(0) == position[0] && temp->r(1) == position[1] && temp->r(2) == position[2])
			{
				found = true;
				break;
			}
			wrapper.iteratorNext();
		}
		CPPUNIT_ASSERT_MESSAGE("Found after insertion", !(found ^ wrapper.isInRegion(position)));

		//deletion
		wrapper.deleteMolecule(tempParticle);

		//verify
		found = false; //sanity
		wrapper.iteratorReset();
		while(wrapper.iteratorValid())
		{
			::Molecule* temp = wrapper.getParticleAtIterator();
			if(temp->r(0) == position[0] && temp->r(1) == position[1] && temp->r(2) == position[2])
			{
				found = true;
				break;
			}
			wrapper.iteratorNext();
		}
		CPPUNIT_ASSERT_MESSAGE("Found after deletion", !found);
	}
	void testPointIsInRegion()
	{
		double start[3] = {0,0,0}, end[3] = {5,5,5};
		ls1::LS1RegionWrapper wrapper(start, end, _testSimulation);
		CPPUNIT_ASSERT(wrapper.isInRegion({1,1,1}));
		CPPUNIT_ASSERT(!wrapper.isInRegion({5,6,5}));
	}
	void testSetRegion()
	{
		double start[3] = {0,0,0}, end[3] = {5,5,5};
		ls1::LS1RegionWrapper wrapper(start, end, _testSimulation);
		CPPUNIT_ASSERT(wrapper.isInRegion({4,4,4}));
		end[0] = end[1] = end[2] = 3;
		wrapper.setRegion(start, end);
		CPPUNIT_ASSERT(!wrapper.isInRegion({4,4,4}));
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
	void testResetIterator()
	{
		ls1::LS1RegionWrapper fullRegion(_testSimulation->getEnsemble()->domain()->rmin(), _testSimulation->getEnsemble()->domain()->rmax(), _testSimulation);
		fullRegion.iteratorReset();
		if(!fullRegion.iteratorValid())
			return; //no particles in subdomain
		
		//grab first particle
		::Molecule* firstParticle = fullRegion.getParticleAtIterator();

		//iterate through and reset
		while(fullRegion.iteratorValid())
		{
			fullRegion.iteratorNext();
		}

		//reset and check first particle
		fullRegion.iteratorReset();
		CPPUNIT_ASSERT( firstParticle->r(0) == fullRegion.getParticleAtIterator()->r(0)
					 && firstParticle->r(1) == fullRegion.getParticleAtIterator()->r(1)
					 && firstParticle->r(2) == fullRegion.getParticleAtIterator()->r(2) );
	}
	void testIteration()
	{
		ls1::LS1RegionWrapper fullRegion(_testSimulation->getEnsemble()->domain()->rmin(), _testSimulation->getEnsemble()->domain()->rmax(), _testSimulation);
		long unsigned int particleCount = 0;
		fullRegion.iteratorReset();
		while(fullRegion.iteratorValid())
		{
			particleCount++;
			fullRegion.iteratorNext();
		}

#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
		long unsigned int globalCount;
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Allreduce(&particleCount, &globalCount, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
		particleCount = globalCount;
#endif
		CPPUNIT_ASSERT( particleCount == _testSimulation->getTotalNumberOfMolecules() );

		//reset
		particleCount = 0;
		fullRegion.iteratorReset();
		while(fullRegion.iteratorValid())
		{
			particleCount++;
			fullRegion.iteratorNext();
		}

#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
		globalCount = 0;
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Allreduce(&particleCount, &globalCount, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
		particleCount = globalCount;
#endif
		CPPUNIT_ASSERT( particleCount == _testSimulation->getTotalNumberOfMolecules() );
	}
private:
	Simulation* _testSimulation;
};

CPPUNIT_TEST_SUITE_REGISTRATION(LS1RegionWrapperTest);