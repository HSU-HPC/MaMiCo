#ifndef LS1_MD_SIMULATION_FACTORY_H_
#define LS1_MD_SIMULATION_FACTORY_H_

#include "coupling/interface/MDSimulationFactory.h"
#include "coupling/interface/impl/ls1/LS1MamicoCouplingSwitch.h"

namespace coupling
{
    namespace interface
    {
        class LS1MDSimulation : public coupling::interface::MDSimulation
        {
        public:
            /** switches coupling on/off*/
            void switchOffCoupling()
            {
                coupling::interface::LS1MamicoCouplingSwitch::getInstance().setCouplingStateOff();
            }
            void switchOnCoupling()
            {
                coupling::interface::LS1MamicoCouplingSwitch::getInstance().setCouplingStateOn();
            }

            /** simulates numberTimesteps time steps and starts at time step no. firstTimestep*/
            virtual void simulateTimesteps(const unsigned int &numberTimesteps, const unsigned int &firstTimestep) = 0;
            /** simulates a single time step*/
            //virtual void simulateTimestep(const unsigned int &thisTimestep ){const unsigned int steps=1; simulateTimesteps(thisTimestep,steps);} TODO BUG
            virtual void sortMoleculesIntoCells() = 0;

            virtual void setMacroscopicCellService(coupling::services::MacroscopicCellService<MDSIMULATIONFACTORY_DIMENSION> *macroscopicCellService) = 0;
            virtual void init() = 0;
            virtual void init(const tarch::utils::MultiMDService<MDSIMULATIONFACTORY_DIMENSION>& multiMDService,unsigned int localMDSimulation) = 0;
            virtual void shutdown() = 0;
        };
    }
}

/** TODO
 * divide domain into ls1regionwrappers as per mamico grid
 * set cutoffs from simulation objects
 */

#endif