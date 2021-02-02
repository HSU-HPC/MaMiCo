#ifndef _Error_ESTIMATION_H_
#define _Error_ESTIMATION_H_



namespace coupling { namespace error {

//  template <unsigned int dim>
  class ErrorEstimation;

  }
}



/** This class provides the basic eqations for analytical error estimation based in the statistical mechanic.
  * There are two typs of error estimation in this class bases on velocity ans density fluctuation.
  * Also it can be used to calculate the required number of MD-simulations as a function of desired max error. 
 */


//template <unsigned int dim>
class coupling::error::ErrorEstimation {
  public:
    ErrorEstimation(const double averageVelocity, const double averageTemperature,
			const double averageNumberOfParticle, const double averageParticleMass,
			double soundSpeed, const int numberOfSamples):
			_averageVelocity(averageVelocity), _averageTemperature(averageTemperature),
			_numberOfParticle(averageNumberOfParticle), _particleMass(averageParticleMass),
			_soundSpeed(soundSpeed), _numberOfSamples(numberOfSamples),
			_desiredErrorDensity(0.05), _desiredErrorVelocity(0.1), _desiredAbsoluteErrorVelocity(0.08){}
      
      ~ErrorEstimation(){
       	#ifdef DEBUG_Error
        std::cout << " Error estimation deconstructed" << std::endl;
        #endif
	}

    enum errorBaseQuantity{Velocity=0, Density=1, Temperature=2};
    enum errorType{Relative=0, Absolute=1};

      //************************************************************************************
      //************************************************************************************

     double getError(errorBaseQuantity BaseQuantity, errorType Type){

        double error;

	if (BaseQuantity==Velocity){

        error = getErrorVelocity(_numberOfSamples, Type==Relative ? _averageVelocity : 1, _averageTemperature, _numberOfParticle, _particleMass);

	}else if (BaseQuantity==Density){

        error = getErrorDensity(_numberOfSamples, _soundSpeed, _averageTemperature, _numberOfParticle, _particleMass);

        }else if (BaseQuantity==Temperature){

        }else{

        std::cout << "ERROR coupling::ErrorEstimation(): Base quantity invalid! " << std::endl;
        exit(EXIT_FAILURE);

	}	

        return error;
      }


      //************************************************************************************
      //************************************************************************************

      double getErrorVelocity(int numberOfSamples, double velocity, double temperature, const double numberOfParticle,const double particleMass ) {
        double error= 1/(velocitySNR( velocity, temperature, numberOfParticle, particleMass)*std::sqrt(numberOfSamples));
        return error;
      }

      double getErrorDensity(int numberOfSamples, double soundSpeed, double temperature, double numberOfParticle,const double particleMass) {
        double refSP = referenceSoundSpeed( temperature, particleMass);
        double Ac = acousticNumber(soundSpeed, refSP);
        double error= 1/std::sqrt(numberOfParticle*numberOfSamples)/Ac;
	std::cout << "numberOfParticle "  << numberOfParticle << "numberOfSamples  " << numberOfSamples << "Ac  "  << Ac <<  std::endl;
        return error;
      }

      double getErrorTemperature(int numberOfSamples, const double numberOfParticle) {

        double error = std::sqrt(_k/(_C_v*numberOfParticle*numberOfSamples));
        return error;
      }


      //************************************************************************************
      //************************************************************************************


      double getCorrectorNumberOfSamples(errorBaseQuantity BaseQuantity, errorType Type){

//        double er = getErrorVelocity();

        if (BaseQuantity==Velocity){

        return requiredSamplesV(Type==Relative ? _desiredErrorVelocity : _desiredAbsoluteErrorVelocity, _averageTemperature, _soundSpeed, Type==Relative ? _averageVelocity : 1, _numberOfParticle, _particleMass);

        }else if (BaseQuantity==Density){

        return requiredSamplesD(_desiredErrorDensity, _averageTemperature, _soundSpeed, _averageVelocity,  _particleMass, _numberOfParticle);

        }else if (BaseQuantity==Temperature){
	
	return 3.0;

        }

	return 4.0;

      }

      //************************************************************************************
      //************************************************************************************

      double requiredSamplesV(double desiredError, double temperature, double soundSpeed, double velocity, double numberOfParticle, double particleMass){ 
//	double refSP= referenceSoundSpeed( temperature, particleMass);
//	double Ac = acousticNumber(soundSpeed, refSP);
	double SNR = velocitySNR( velocity, temperature,  numberOfParticle, particleMass);
	double desiredNumber = 1/((SNR*desiredError)*(SNR*desiredError));

	return desiredNumber;
	}


      double requiredSamplesD(double desiredError, double temperature, double soundSpeed, double velocity,  double particleMass, double numberOfParticle){
        double refSP= referenceSoundSpeed( temperature, particleMass);
        double Ac = acousticNumber(soundSpeed, refSP);
        double desiredNumber = (numberOfParticle*desiredError)*(numberOfParticle*desiredError)*Ac;

        return desiredNumber;
        }

      double reqiredSamplesT(double desiredError){

        double desiredNumber = (_k/(_C_v*_numberOfParticle*desiredError))*(_k/(_C_v*_numberOfParticle*desiredError));

        return desiredNumber;
        }
      //*******************************************************************************
      //*******************************************************************************
      const double velocitySNR(double velocity, double temperature, const double numberOfParticle,const double particleMass){
	return velocity/veloyityFluctuation( temperature, numberOfParticle, particleMass); }  // velocity deviation   

      const double veloyityFluctuation(double temperature, const double numberOfParticle,const double particleMass){
	return std::sqrt(_k*temperature/particleMass/numberOfParticle);}  // velocity deviation


      //*************************************************************************************
      void setAbsVelocityError(double error){ _desiredAbsoluteErrorVelocity=error;}

      void setAbsDensityError(double error){ _desiredErrorDensity=error;}

      void setRelVelocityError(double error){ _desiredErrorVelocity=error;}

 private:


      double getErrorVelocity(){
        double er = getErrorVelocity(_numberOfSamples, _averageVelocity, _averageTemperature, _numberOfParticle, _particleMass);
        return er;
      }

      double getErrorTemperature(){
        double er = getErrorTemperature(_numberOfSamples, _numberOfParticle);
        return er;
      }

      double getErrorDensity(){
        double er = getErrorDensity(_numberOfSamples, _soundSpeed, _averageTemperature, _numberOfParticle, _particleMass);
        return er;
      }

        double requiredSamplesV(){
        return requiredSamplesV(_desiredErrorVelocity, _averageTemperature, _soundSpeed, _averageVelocity, _numberOfParticle, _particleMass);
        }

//**********************************************************
//***************************************************************
   	
      double _averageVelocity;
      double _averageTemperature;
      double _numberOfParticle;
      double _particleMass;
      double _gamma=1.667;
      double _soundSpeed;
      int _numberOfSamples;
      
    double _desiredErrorDensity; //=0.05;
    double _desiredErrorVelocity; // =0.1;
    double _desiredAbsoluteErrorVelocity; //=0.08;

      double _C_v=1.4973;
      double _k=1;            // Boltzman Constant unit: [J/K]
    ////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////

    double acousticNumber(double soundSpeed, double soundSpeed_ref){return (soundSpeed>0 ? soundSpeed/soundSpeed_ref : 1);}
    double referenceSoundSpeed(double temperature, double particleMass){return std::sqrt(_gamma*_k*temperature/particleMass);}
   

};

#endif
