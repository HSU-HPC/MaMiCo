# *********************************************************************************************
# README
# *********************************************************************************************
#
# Available configuration Parameters for ExaHyPE
#
# export variable  |  default-value  |  further values         |  description
#-------------------------------------------------------------------------------------------------------
# COMPILER            Intel             GNU                          Used compiler (and linker)
# MODE                Release           Debug,                       Verbosity and Debug level
#                                       Asserts,PeanoProfile,Profile          
# SHAREDMEM           None              OMP, TBB                     Shared-memory parallelisation
# DISTRIBUTEDMEM      MPI               None                         Distributed-memory parallelisation
# BOUNDARYCONDITIONS  None              Periodic                     Type of boundary conditions
# USE_IPO             Off               On, Unsafe                   IPO during compilation (intel only)
# ******************************************************************************************************

# Helper
SHELL = bash
tolower = $(shell echo $(1) | tr '[:upper:]' '[:lower:]')

# Default Paths and name
PEANO_KERNEL_PEANO_PATH=/home/hpc/a2c06/lu23vov/exahype-engine/./Peano/peano
PEANO_KERNEL_TARCH_PATH=/home/hpc/a2c06/lu23vov/exahype-engine/./Peano/tarch
PEANO_TOOLBOX_MULTISCALELINKEDCELL_PATH=/home/hpc/a2c06/lu23vov/exahype-engine/./Peano/multiscalelinkedcell
PEANO_TOOLBOX_SHAREDMEMORY_ORACLES_PATH=/home/hpc/a2c06/lu23vov/exahype-engine/./Peano/sharedmemoryoracles
PEANO_TOOLBOX_MPI_BLANCING_PATH=/home/hpc/a2c06/lu23vov/exahype-engine/./Peano/mpibalancing
EXAHYPE_PATH=/home/hpc/a2c06/lu23vov/exahype-engine/./ExaHyPE
PROJECT_PATH=/home/hpc/a2c06/lu23vov/exahype-engine/../mamico-noisereduction/coupling/solvers/exahype
EXECUTABLE=ExaHyPE-couette

# Default flags from specification file constants
PROJECT_CFLAGS+= -DDim3
ARCHITECTURE=CPU

# Default sharedmem/distributedmem required by the specification file
ifeq ($(SHAREDMEM),)
    SHAREDMEM=None
endif
ifeq ($(DISTRIBUTEDMEM),)
    DISTRIBUTEDMEM=MPI
endif

# If using TBB only, link pthread (comment this if you want to link a specific pthread)
ifeq ($(call tolower,$(SHAREDMEM)),tbb)
    ifeq ($(call tolower,$(DISTRIBUTEDMEM)),none)
        PROJECT_LFLAGS+= -pthread
    endif
endif

# Default options required by the specification file



# Main Makefile include
-include /home/hpc/a2c06/lu23vov/exahype-engine/./ExaHyPE/Makefile


all: 
	@echo couette
	@echo =================
	@echo An ExaHyPE solver