#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

#include "coupling/CouplingMDDefinitions.h"
#include "coupling/ErrorEstimation.h"

/** tests several properties of the ErrorEstimation. 
 *  @author Vahid Jafarimann
 */
class ErrorEstimationTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(ErrorEstimationTest);
  CPPUNIT_TEST(testVelocityError);
  CPPUNIT_TEST(testDensityError);
  CPPUNIT_TEST(testTemperatureError);
  CPPUNIT_TEST(testPressureError);
  CPPUNIT_TEST_SUITE_END();

public:
  void testVelocityError() {
    double velocity = 5;
    double temperature = 18;
    double numberOfParticle = 16;
    double particleMass = 2.0;
    double soundSpeed = 1.0;
    double numberOfSamples = 9;
    double cellVolume = 1.0;
    std::cout << "Run velocityError test..." << std::endl;
    _testVelocityError(velocity, temperature, numberOfParticle, particleMass, soundSpeed, numberOfSamples, cellVolume);
  }
  
  void testDensityError() {
    double velocity = 5;
    double temperature = 16;
    double numberOfParticle = 25;
    double particleMass = 1.667;
    double soundSpeed = 1.0;
    double numberOfSamples = 64;
    double cellVolume = 1.0;
    std::cout << "Run densityError test..." << std::endl;
    _testDensityError(velocity, temperature, numberOfParticle, particleMass, soundSpeed, numberOfSamples, cellVolume);
  }
  
  void testTemperatureError() {
    double velocity = 5;
    double temperature = 16;
    double numberOfParticle = 66.78688306;
    double particleMass = 1.667;
    double soundSpeed = 1.0;
    double numberOfSamples = 100;
    double cellVolume = 1.0;
    std::cout << "Run temperatureError test..." << std::endl;
    _testTemperatureError(velocity, temperature, numberOfParticle, particleMass, soundSpeed, numberOfSamples, cellVolume);
  }
  
  void testPressureError() {
    double velocity = 5;
    double temperature = 16;
    double numberOfParticle = 25;
    double particleMass = 1.667;
    double soundSpeed = 1.0;
    double numberOfSamples = 1.667;
    double cellVolume = 1.0;
    std::cout << "Run pressureError test..." << std::endl;
    _testPressureError(velocity, temperature, numberOfParticle, particleMass, soundSpeed, numberOfSamples, cellVolume);
  }


private:
  void _testVelocityError(double velocity, double temperature, double numberOfParticle, double particleMass, double soundSpeed, double numberOfSamples, double cellVolume) {
      
    coupling::error::ErrorEstimation errorControl(velocity, temperature, numberOfParticle, particleMass, soundSpeed, numberOfSamples, cellVolume);
    double errorVelo = errorControl.getError(coupling::error::ErrorEstimation::Velocity, coupling::error::ErrorEstimation::Absolute);
    CPPUNIT_ASSERT(errorVelo == 0.25);
    
    errorVelo = errorControl.getError(coupling::error::ErrorEstimation::Velocity, coupling::error::ErrorEstimation::Relative);
    CPPUNIT_ASSERT(errorVelo == 0.05);
    
    errorControl.setAbsVelocityError(0.05);
    double NoMD = errorControl.getCorrectorNumberOfSamples(coupling::error::ErrorEstimation::Velocity, coupling::error::ErrorEstimation::Absolute);
    CPPUNIT_ASSERT(NoMD == 225.0);
    
    errorControl.setRelVelocityError(0.01);
    NoMD = errorControl.getCorrectorNumberOfSamples(coupling::error::ErrorEstimation::Velocity, coupling::error::ErrorEstimation::Relative);
    CPPUNIT_ASSERT(NoMD == 225.0);
  }

  void _testDensityError(double velocity, double temperature, double numberOfParticle, double particleMass, double soundSpeed, double numberOfSamples, double cellVolume) {
      
    coupling::error::ErrorEstimation errorControl(velocity, temperature, numberOfParticle, particleMass, soundSpeed, numberOfSamples, cellVolume);
    double errorD = errorControl.getError(coupling::error::ErrorEstimation::Density, coupling::error::ErrorEstimation::Relative);
    CPPUNIT_ASSERT(errorD == 0.1);
    
    errorD = errorControl.getError(coupling::error::ErrorEstimation::Density, coupling::error::ErrorEstimation::Absolute);
    CPPUNIT_ASSERT(errorD == 2.5);
    
    errorControl.setRelDensityError(0.2);
    double NoMD = errorControl.getCorrectorNumberOfSamples(coupling::error::ErrorEstimation::Density, coupling::error::ErrorEstimation::Relative);
    CPPUNIT_ASSERT(std::abs(NoMD- 16)<0.000001);
    
  }
  
  void _testTemperatureError(double velocity, double temperature, double numberOfParticle, double particleMass, double soundSpeed, double numberOfSamples, double cellVolume) {
      
    coupling::error::ErrorEstimation errorControl(velocity, temperature, numberOfParticle, particleMass, soundSpeed, numberOfSamples, cellVolume);
    double errorT = errorControl.getError(coupling::error::ErrorEstimation::Temperature, coupling::error::ErrorEstimation::Relative);
    CPPUNIT_ASSERT(std::abs(errorT - 0.01)<0.000001);
    
    errorT = errorControl.getError(coupling::error::ErrorEstimation::Temperature, coupling::error::ErrorEstimation::Absolute);
    CPPUNIT_ASSERT(std::abs(errorT - 0.16)<0.000001);
    
    errorControl.setRelDensityError(0.2);
    double NoMD = errorControl.getCorrectorNumberOfSamples(coupling::error::ErrorEstimation::Temperature, coupling::error::ErrorEstimation::Relative);
    CPPUNIT_ASSERT(std::abs(NoMD - 3)<0.000001);
    
  }
  
  void _testPressureError(double velocity, double temperature, double numberOfParticle, double particleMass, double soundSpeed, double numberOfSamples, double cellVolume) {
      
    coupling::error::ErrorEstimation errorControl(velocity, temperature, numberOfParticle, particleMass, soundSpeed, numberOfSamples, cellVolume);
    double errorP = errorControl.getError(coupling::error::ErrorEstimation::Pressure, coupling::error::ErrorEstimation::Relative);
    CPPUNIT_ASSERT(std::abs(errorP - 1.0)<0.000001);
    
    errorControl.setRelDensityError(0.2);
    double NoMD = errorControl.getCorrectorNumberOfSamples(coupling::error::ErrorEstimation::Pressure, coupling::error::ErrorEstimation::Relative);
    CPPUNIT_ASSERT(std::abs(NoMD - 4)<0.000001);
    
  }

 
};

CPPUNIT_TEST_SUITE_REGISTRATION(ErrorEstimationTest);
