#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

#include "coupling/interface/impl/ls1/LS1MDSolverInterface.h"

#include "tarch/tinyxml2/tinyxml2.h"

#include <sstream>

#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
#include <mpi.h>
#endif

#define MY_LINKEDCELL ls1::LS1RegionWrapper

class LS1MDSolverInterfaceTest : public CppUnit::TestFixture
{
	CPPUNIT_TEST_SUITE(LS1MDSolverInterfaceTest);
	CPPUNIT_TEST(testConstants);
	CPPUNIT_TEST(testAddAndDeleteParticle);
	CPPUNIT_TEST(testGetCellAndIterator);
	CPPUNIT_TEST_SUITE_END();

public:
	void setUp()
	{
		global_log = new Log::Logger(Log::None);
		global_log->set_mpi_output_root(0);
		coupling::interface::LS1StaticCommData::getInstance().setBoxOffsetAtDim(0,0);
		coupling::interface::LS1StaticCommData::getInstance().setBoxOffsetAtDim(1,0);
		coupling::interface::LS1StaticCommData::getInstance().setBoxOffsetAtDim(2,0);
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
		coupling::interface::LS1StaticCommData::getInstance().setLocalCommunicator(MPI_COMM_WORLD);
#endif
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
	void testConstants()
	{
		//create interface
		coupling::interface::LS1MDSolverInterface interface({5,5,5},{1,1,1});
		//load config file to read consts
		tinyxml2::XMLDocument conffile;
		tinyxml2::XMLElement* siteInfo = nullptr;
		conffile.LoadFile("../test/unit/coupling/interface/impl/ls1/ls1gridconfig.xml");
		tinyxml2::XMLHandle fileHandle(conffile);

		siteInfo = fileHandle.FirstChildElement("mardyn").FirstChildElement("simulation").FirstChildElement("ensemble").FirstChildElement("components").FirstChildElement("moleculetype").FirstChildElement("site").ToElement();
		if(siteInfo)
		{
			double mass,sigma,epsilon;
			siteInfo->FirstChildElement("mass")->QueryDoubleText(&mass);
			siteInfo->FirstChildElement("sigma")->QueryDoubleText(&sigma);
			siteInfo->FirstChildElement("epsilon")->QueryDoubleText(&epsilon);
			CPPUNIT_ASSERT( mass == interface.getMoleculeMass() );
			CPPUNIT_ASSERT( sigma == interface.getMoleculeSigma() );
			CPPUNIT_ASSERT( epsilon == interface.getMoleculeEpsilon() );
			CPPUNIT_ASSERT( interface.getKB() == 1.0);
		}
		else
		{
			CPPUNIT_ASSERT_MESSAGE( "XML node could not be found", false );
		}
		
		siteInfo = fileHandle.FirstChildElement("mardyn").FirstChildElement("simulation").FirstChildElement("integrator").FirstChildElement("timestep").ToElement();
		if(siteInfo)
		{
			double timestep;
			siteInfo->QueryDoubleText(&timestep);
			CPPUNIT_ASSERT( timestep == interface.getDt() );
		}
		else
		{
			CPPUNIT_ASSERT_MESSAGE( "XML node could not be found", false );
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
		//create interface
		coupling::interface::LS1MDSolverInterface interface({5,5,5},{1,1,1});

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
		::Molecule blank;
		coupling::interface::LS1Molecule tempParticle(&blank);
		tempParticle.setPosition(position);
		interface.addMoleculeToMDSimulation(tempParticle);

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
		CPPUNIT_ASSERT_MESSAGE("Found after insertion", found == wrapper.isInRegion(position));

		//deletion
		interface.deleteMoleculeFromMDSimulation(tempParticle, wrapper);

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
	void testGetCellAndIterator()
	{
		//create interface with macroscopic cell size 10 and 2 linked cells per mac.cell, hence linked cells are size 5,5,5
		coupling::interface::LS1MDSolverInterface interface({10,10,10},{2,2,2});

	}
private:
	Simulation* _testSimulation;
};

CPPUNIT_TEST_SUITE_REGISTRATION(LS1MDSolverInterfaceTest);