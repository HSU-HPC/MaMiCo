#ifndef _MOLECULARDYNAMICS_COUPLING_NOISEREDUCTION_IDENTITY_H_
#define _MOLECULARDYNAMICS_COUPLING_NOISEREDUCTION_IDENTITY_H_

#include "coupling/noisereduction/NoiseReduction.h"

namespace coupling {
  namespace noisereduction {
    template<unsigned int dim>
    class IdentityTransform;
  }
}

/** Dummy noise reduction class for not performing any filtering at all. Data from MD is kept unmodified
 *
 *  @author Piet Jarmatz
 */
template<unsigned int dim>
class coupling::noisereduction::IdentityTransform:
public coupling::noisereduction::NoiseReduction<dim> {
  public:
    IdentityTransform(const coupling::IndexConversion<dim> &indexConversion, const tarch::utils::MultiMDService<dim>& multiMDService):
      coupling::noisereduction::NoiseReduction<dim>(indexConversion, multiMDService){}
      
    virtual ~IdentityTransform(){}
};
#endif // _MOLECULARDYNAMICS_COUPLING_NOISEREDUCTION_IDENTITY_H_
