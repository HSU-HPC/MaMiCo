#ifndef _MOLECULARDYNAMICS_COUPLING_NOISEREDUCTION_IDENTITY_H_
#define _MOLECULARDYNAMICS_COUPLING_NOISEREDUCTION_IDENTITY_H_

#include "coupling/noisereduction/NoiseReduction.h"

namespace coupling {
  namespace noisereduction {
    template<class LinkedCell, unsigned int dim>
    class IdentityTransform;
  }
}

/** Dummy noise reduction class for not performing any filtering at all. Data from MD is kept unmodified
 *
 *  @author Piet Jarmatz
 */
template<class LinkedCell,unsigned int dim>
class coupling::noisereduction::IdentityTransform:
public coupling::noisereduction::NoiseReduction<LinkedCell,dim> {
  public:
    IdentityTransform(const coupling::IndexConversion<dim> &indexConversion):
      coupling::noisereduction::NoiseReduction<LinkedCell,dim>(indexConversion){}
      
    virtual ~IdentityTransform(){}
};
#endif // _MOLECULARDYNAMICS_COUPLING_NOISEREDUCTION_IDENTITY_H_
