// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#pragma once
#include <string>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#define DEBUG_POD
#include "coupling/filtering/FilterInterface.h"
#include "tarch/utils/MultiMDService.h"

namespace coupling {
    template<unsigned int dim>
    class POD;
}

/** Noise reduction algorithm using Proper orthogonal decomposition
 *
 *  @author Piet Jarmatz, Felix Maurer
 */
template<unsigned int dim>
class coupling::POD : public coupling::FilterInterface<dim>{
    public:
        POD(  	const std::vector<coupling::datastructures::MacroscopicCell<dim> *>& inputCellVector,
				const std::vector<coupling::datastructures::MacroscopicCell<dim> *>& outputCellVector,
				const std::vector<tarch::la::Vector<dim, unsigned int>> cellIndices,
				bool filteredValues[7],
				const tarch::utils::MultiMDService<dim>& multiMDService,
				int tws,
				int kmax
				):
				coupling::FilterInterface<dim>(inputCellVector, outputCellVector, cellIndices, filteredValues),
				_multiMDService(multiMDService),
				_timeWindowSize(tws),
				_kMax(kmax),
				_cycleCounter(0),
				_spatialIndex(0),
				_t(0),
				_data(NULL),
				_C(NULL),
				_A(NULL),
				_A_T(NULL)
		{
			int spatialSize = cellIndices.size();
			_data = new Eigen::MatrixXd[dim+1]; // separate data matrices for: mass, momentum0, momentum1, momentum2
			_C = new Eigen::MatrixXd[dim+1];
 			_A = new Eigen::MatrixXd[dim+1];
  			_A_T = new Eigen::MatrixXd[dim+1];
  			for(unsigned int i=0;i<dim+1;i++){
    			_data[i] = Eigen::MatrixXd::Constant(_timeWindowSize, spatialSize, (i==0)?1:0);
 			}

        	#ifdef DEBUG_POD
			std::cout << "		POD: Created Proper Orthogonal Decomposition instance." << std::endl;
			//TODO selection of filtered properties
			std::cout << "			WARNING: Regardless of configuration, POD always filters macroscopic mass and momentum." << std::endl;
       		#endif
        }

        ~POD(){
			if (_data!=NULL){ delete [] _data; _data=NULL;}
			if (_C!=NULL){ delete [] _C; _C=NULL;}
			if (_A!=NULL){ delete [] _A; _A=NULL;}
 			if (_A_T!=NULL){ delete [] _A_T; _A_T=NULL;}

        	#ifdef DEBUG_POD
			std::cout << "		POD: Deleted Proper Orthogonal Decomposition instance." << std::endl;
        	#endif
        }

     
	    void operator()();
	private:
		const tarch::utils::MultiMDService<dim>& _multiMDService;
		unsigned int _timeWindowSize; // number of snapshots / coupling cycles taken into consideration for noise reduction
    	const unsigned int _kMax; // number of dominant eigenvalues
    	unsigned int _cycleCounter; // coupling cycle counter, indicates how many data snapshots are available already
    	unsigned int _spatialIndex; // cell counter, should run from zero to getLocalNumberMacroscopicCells()-1 within an iteration of ProcessInnerMacroscopicCell
    	unsigned int _t; // active temporal index, iterates cyclic between zero and _timeWindowSize
    	Eigen::MatrixXd *_data;  // set of snapshots (sampled by transferStrategy)
    	Eigen::MatrixXd *_C;  // temporal auto-correlation covariance matrix of _data
    	Eigen::MatrixXd *_A;  // POD temporal modes / eigenvectors of C
    	Eigen::MatrixXd *_A_T;  // Transpose of A


};

//include implementation of header
#include "POD.cpph"


