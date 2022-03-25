#ifndef _Error_ESTIMATION_H_
#define _Error_ESTIMATION_H_

namespace coupling {
namespace error {

class ErrorEstimation;

}
} // namespace coupling

/** This class provides the basic eqations for analytical error estimation based in the statistical mechanic.
 * 	There are two typs of error estimation in this class bases on velocity ans density fluctuation.
 *	Also it can be used to estimate the required number of MD-simulations as a function of required maximum error.
 *	@brief This class is used to analyse the error and predict the required number of instances.
 *  @author Vahid Jafarittmer
 */
class coupling::error::ErrorEstimation {
public:
  /** Constructor: Initializes the private member of the class.
   * 	@param velocity
   * 	@param temperature
   * 	@param numberOfParticle
   * 	@param particleMass
   * 	@param soundSpeed
   * 	@param numberOfSamples
   * 	@param cellVolume
   */
  ErrorEstimation(double velocity, double temperature, double numberOfParticle, double particleMass, double soundSpeed, int numberOfSamples, double cellVolume)
      : _velocity(velocity), _temperature(temperature), _numberOfParticle(numberOfParticle), _particleMass(particleMass), _soundSpeed(soundSpeed),
        _numberOfSamples(numberOfSamples), _cellVolume(cellVolume), _desiredAbsErrorVelocity(0.08), _desiredRelErrorVelocity(0.1),
        _desiredAbsErrorDensity(0.05), _desiredRelErrorDensity(0.05), _desiredAbsErrorTemperature(0.05), _desiredRelErrorTemperature(0.05),
        _desiredAbsErrorPressure(0.05), _desiredRelErrorPressure(0.05), _gamma(1.667), _C_v(1.4973), _k(1) {}

  /** Destructor:
   */
  ~ErrorEstimation() {
#ifdef DEBUG_Error
    std::cout << " Error estimation deconstructed" << std::endl;
#endif
  }

  /** The quantiy, its error we want to analyse
   *	@enum errorBaseQuantity
   */
  enum errorBaseQuantity {
    Velocity = 0 /**< error in Velocity*/,
    Density = 1 /**< error in Density*/,
    Temperature = 2 /**< error in Temperature*/,
    Pressure = 3 /**< error in Pressure*/
  };
  /** relative/absolute error
   *	@enum errorType
   */
  enum errorType { Relative = 0 /**< Relative error */, Absolute = 1 /**< Absolute error */ };

  /** This function predict the relative or absolute error of the quantites veocty, density, tempetature or pressure based on its arguments.
   *	 @param errorBaseQuantity
   *	 @param errorType
   */
  double getError(errorBaseQuantity BaseQuantity, errorType Type) {

    double error;

    if (BaseQuantity == Velocity) {

      error = getErrorVelocity(_numberOfSamples, Type == Relative ? _velocity : 1, _temperature, _numberOfParticle, _particleMass);

    } else if (BaseQuantity == Density) {

      error = getErrorDensity(_numberOfSamples, _soundSpeed, _temperature, Type == Relative ? _numberOfParticle : 1 / _numberOfParticle, _particleMass);

    } else if (BaseQuantity == Temperature) {

      error = getErrorTemperature(_numberOfSamples, _numberOfParticle, Type == Relative ? 1 : _temperature);

    } else if (BaseQuantity == Pressure) {

      error = getErrorPressure(_numberOfSamples, _numberOfParticle, _temperature, _soundSpeed, _particleMass, _cellVolume, 1);

    } else {

      std::cout << "ERROR coupling::ErrorEstimation(): Base quantity invalid! " << std::endl;
      exit(EXIT_FAILURE);
    }

    return error;
  }

  /** This function predict the required number of MD instances to keep the relative or absolute error of the quantites veocty, density, tempetature or pressure
   *based on its arguments under a certain value. The desired error value has to be set before, this function is called.
   *	 @param errorBaseQuantity
   *	 @param errorType
   */
  double getCorrectorNumberOfSamples(errorBaseQuantity BaseQuantity, errorType Type) {

    if (BaseQuantity == Velocity) {

      return requiredSamplesV(Type == Relative ? _desiredRelErrorVelocity : _desiredAbsErrorVelocity, _temperature, _soundSpeed,
                              Type == Relative ? _velocity : 1, _numberOfParticle, _particleMass);

    } else if (BaseQuantity == Density) {

      return requiredSamplesD(_desiredAbsErrorDensity, _temperature, _soundSpeed, _velocity, _particleMass, _numberOfParticle);

    } else if (BaseQuantity == Temperature) {

      return 3.0;

    } else if (BaseQuantity == Pressure) {

      return 4.0;
    }

    std::cout << "Error Estimation::getCorrectorNumberOfSamples not a Valid quantity" << std::endl;
    std::exit(EXIT_FAILURE);
    return 0.0;
  }

  /** This function estimates the velocity error
   *	 @param numberOfSamples
   *	 @param velocity
   *	 @param temperature
   *	 @param numberOfParticle
   *	 @param particleMass
   *	 @return error veocty error
   */
  double getErrorVelocity(int numberOfSamples, double velocity, double temperature, double numberOfParticle, double particleMass) {
    double error = 1 / (velocitySNR(velocity, temperature, numberOfParticle, particleMass) * std::sqrt(numberOfSamples));
    return error;
  }

  /** This function estimates the density error
   *	 @param numberOfSamples
   *	 @param soundSpeed
   *	 @param temperature
   *	 @param numberOfParticle
   *	 @param particleMass
   *	 @return error density error
   */
  double getErrorDensity(int numberOfSamples, double soundSpeed, double temperature, double numberOfParticle, double particleMass) {
    double refSP = referenceSoundSpeed(temperature, particleMass);
    double Ac = acousticNumber(soundSpeed, refSP);
    double error = 1 / std::sqrt(numberOfParticle * numberOfSamples) / Ac;
    //	std::cout << "numberOfParticle "  << numberOfParticle << "numberOfSamples  " << numberOfSamples << "Ac  "  << Ac <<  std::endl;
    return error;
  }

  /** This function estimates the temperature error
   *	 @param numberOfSamples
   *	 @param temperature
   *	 @param numberOfParticle
   *	 @return error temperature error
   */
  double getErrorTemperature(int numberOfSamples, double numberOfParticle, double temperature) {

    //	    double deltaT = k*temperature*temperature/_C_v/numberOfParticle;

    double error = std::sqrt(_k / (_C_v * numberOfParticle * numberOfSamples)) * temperature;
    return error;
  }

  /** This function estimates the pressure error
   *	 @param numberOfSamples
   *	 @param numberOfParticle
   *	 @param temperature
   *	 @param soundSpeed
   *	 @param particleMass
   *	 @param cellVolume
   *	 @param pressure
   *	 @return error pressure error
   */
  double getErrorPressure(int numberOfSamples, double numberOfParticle, double temperature, double soundSpeed, double particleMass, double cellVolume,
                          double pressure) {

    double Ac = acousticNumber(soundSpeed, referenceSoundSpeed(temperature, particleMass));
    double referenceP = referencePressure(temperature, numberOfParticle, cellVolume);

    double error = referenceP / pressure * Ac * std::sqrt(_gamma / (numberOfParticle * numberOfSamples)) * temperature;
    return error;
  }

  /** This function estimates the number of MD instances required to keep the error of the veloctiy under _desiredAbsErrorVelocity or _desiredrelErrorVelocity
   *	 @param desiredError
   *	 @param temperature
   *	 @param soundSpeed
   *	 @param velocity
   *	 @param numberOfParticle
   *	 @param particleMass
   *	 @return desiredNumber required number of MD instances
   */
  double requiredSamplesV(double desiredError, double temperature, double soundSpeed, double velocity, double numberOfParticle, double particleMass) {
    //	double refSP= referenceSoundSpeed( temperature, particleMass);
    //	double Ac = acousticNumber(soundSpeed, refSP);
    double SNR = velocitySNR(velocity, temperature, numberOfParticle, particleMass);
    double desiredNumber = 1 / ((SNR * desiredError) * (SNR * desiredError));

    return desiredNumber;
  }

  /** This function estimates the number of MD instances required to keep the error of the density under _desiredAbsErrorDensity or _desiredrelErrorDensity
   *	 @param desiredError
   *	 @param temperature
   *	 @param soundSpeed
   *	 @param velocity
   *	 @param numberOfParticle
   *	 @param particleMass
   *	 @return desiredNumber required number of MD instances
   */
  double requiredSamplesD(double desiredError, double temperature, double soundSpeed, double velocity, double particleMass, double numberOfParticle) {
    double refSP = referenceSoundSpeed(temperature, particleMass);
    double Ac = acousticNumber(soundSpeed, refSP);
    double desiredNumber = (numberOfParticle * desiredError) * (numberOfParticle * desiredError) * Ac;

    return desiredNumber;
  }

  /** This function estimates the number of MD instances required to keep the error of the temperature under _desiredAbsErrorTemperature or
   *_desiredrelErrorTemperature
   *	 @param desiredError
   *	 @param temperature
   *	 @param soundSpeed
   *	 @param velocity
   *	 @param numberOfParticle
   *	 @param particleMass
   *	 @return desiredNumber required number of MD instances
   */
  double reqiredSamplesT(double desiredError) {

    double desiredNumber = (_k / (_C_v * _numberOfParticle * desiredError)) * (_k / (_C_v * _numberOfParticle * desiredError));

    return desiredNumber;
  }

  /** signal-to-noise ratio as the average fluid velocity over the standard deviation (square root of the veloyity fluctuation=
   *	 @param velocity
   *	 @param temperature
   *	 @param numberOfParticle
   *	 @param particleMass
   *	 @return desiredNumber required number of MD instances
   *	 @sa veloyityFluctuation
   */
  double velocitySNR(double velocity, double temperature, double numberOfParticle, double particleMass) {

    return velocity / veloyityFluctuation(temperature, numberOfParticle, particleMass);
  }

  /** veloyity fluctuation from equilibrium statistical mechanics
   *	 @param temperature
   *	 @param numberOfParticle
   *	 @param particleMass
   *	 @return veloyity fluctuation
   *	 @sa velocitySNR
   */
  double veloyityFluctuation(double temperature, double numberOfParticle, double particleMass) {

    return std::sqrt(_k * temperature / particleMass / numberOfParticle);
  }

  /** sets the required absloute error of the velocity @param error */
  void setAbsVelocityError(double error) { _desiredAbsErrorVelocity = error; }

  /** sets the required relative error of the velocity @param error */
  void setRelVelocityError(double error) { _desiredRelErrorVelocity = error; }

  /** sets the required absloute error of the density @param error */
  void setAbsDensityError(double error) { _desiredAbsErrorDensity = error; }

  /** sets the required relative error of the density @param error */
  void setRelDensityError(double error) { _desiredRelErrorDensity = error; }

  /** sets the required absloute error of the temperature @param error */
  void setAbsTemperatureError(double error) { _desiredAbsErrorTemperature = error; }

  /** sets the required relative error of the temperature @param error */
  void setRelTemperatureError(double error) { _desiredRelErrorTemperature = error; }

  /** sets the required absloute error of the pressure @param error */
  void setAbsPressureError(double error) { _desiredAbsErrorPressure = error; }

  /** sets the required relative error of the pressure @param error */
  void setRelPressureError(double error) { _desiredRelErrorPressure = error; }

private:
  /** This function estimates the velocity error
   *	 @return velocity error
   */
  double getErrorVelocity() {
    double er = getErrorVelocity(_numberOfSamples, _velocity, _temperature, _numberOfParticle, _particleMass);
    return er;
  }

  /** This function estimates the density error
   *	 @return density error
   */
  double getErrorDensity() {
    double er = getErrorDensity(_numberOfSamples, _soundSpeed, _temperature, _numberOfParticle, _particleMass);
    return er;
  }

  /** This function estimates the temperature error
   *	 @return temperature error
   */
  double getErrorTemperature() {
    double er = getErrorTemperature(_numberOfSamples, _numberOfParticle, _temperature);
    return er;
  }

  double acousticNumber(double soundSpeed, double soundSpeed_ref) { return (soundSpeed / soundSpeed_ref); }

  //------------------------- Reference Properties ---------------------------------------------------
  /** calculates the sound speed of a reference ideal gas at the same temperature
   * 	@param temperature
   * 	@param particleMass
   * 	@return referenc sound speed
   */
  double referenceSoundSpeed(double temperature, double particleMass) { return std::sqrt(_gamma * _k * temperature / particleMass); }

  /** calculates thepressure of an ideal gas under the same conditions
   * 	@param temperature
   * 	@param numberOfParticle
   * 	@param cellVolume
   * 	@return reference pressure
   */
  double referencePressure(double temperature, double numberOfParticle, double cellVolume) {
    return std::sqrt(numberOfParticle * _k * temperature / cellVolume);
  }

  double _velocity;
  double _temperature;
  double _numberOfParticle;
  double _particleMass;
  double _soundSpeed;
  int _numberOfSamples;
  double _cellVolume;

  double _desiredAbsErrorVelocity; //=0.08;
  double _desiredRelErrorVelocity; // =0.1;

  double _desiredAbsErrorDensity; //=0.05;
  double _desiredRelErrorDensity; //=0.05;

  double _desiredAbsErrorTemperature; //=0.05;
  double _desiredRelErrorTemperature; //=0.05;

  double _desiredAbsErrorPressure; //=0.05;
  double _desiredRelErrorPressure; //=0.05;

  double _gamma;
  double _C_v;
  /** Boltzman Constant unit: [J/K] */
  double _k;
};

#endif
