#ifndef _MOLECULARDYNAMICS_COUPLING_INTERFACE_SIMPLEMDLINKEDCELLWRAPPER_H_
#define _MOLECULARDYNAMICS_COUPLING_INTERFACE_SIMPLEMDLINKEDCELLWRAPPER_H_

#include "simplemd/MolecularDynamicsDefinitions.h"

namespace coupling {
namespace interface {
class SimpleMDLinkedCellWrapper;
}
} // namespace coupling

/** interface to access molecule of SimpleMD by substituting the simplemd::LinedCell for a plain integer cell index.
 */
class coupling::interface::SimpleMDLinkedCellWrapper {
public:
    SimpleMDLinkedCellWrapper(const tarch::la::Vector<MD_DIM, unsigned int> cellIndex) : _cellIndex(cellIndex) {}

    inline const tarch::la::Vector<MD_DIM, unsigned int> getCellIndex() const { return _cellIndex; } 

private:
    const tarch::la::Vector<MD_DIM, unsigned int> _cellIndex;
};
#endif /* _MOLECULARDYNAMICS_COUPLING_INTERFACE_SIMPLEMDLINKEDCELLWRAPPER_H_ */
