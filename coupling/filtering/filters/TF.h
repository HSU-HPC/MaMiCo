#pragma once

#include "coupling/filtering/interfaces/FilterInterface.h"


template<unsigned int dim>
class coupling::TF : public coupling::FilterInterface<dim>{
    public:
        Gauss(){}

     
	void operator()();

	private:

};

