# Copyright (C) 2015 Technische Universitaet Muenchen
# This file is part of the Mamico project. For conditions of distribution
# and use, please see the copyright notice in Mamico's main folder, or at
# www5.in.tum.de/mamico
import os;
import sys;

#
##### Prints out the available options for command line parameter "target"
#
def printTargets():   
   print
   print "       Possible target-names: "
   print "         simplemd       : simple molecular dynamics simulation"
   print "         libsimplemd    : SimpleMD built as library"


#########################################################################
##### MAIN CODE
#########################################################################

##### Initialize build variables
#
cxx = ''
cppdefines = []
cpppath = ['./']
ccflags = []
linkerflags = []
libpath = []
libs = []
sourcesMolecularDynamics = []

libs.append('rt')
ccflags.append('-std=c++0x')

##### Determine dimension for which to build
#
dim = ARGUMENTS.get('dim', 2)   # Read command line parameter
if int(dim) == 2:
   cppdefines.append('MDDim2')
elif int(dim) == 3:
   cppdefines.append('MDDim3')
else:
   print "ERROR: dim must be either '2' or '3'!"
   sys.exit(1)

##### Add build parameter specific build variable settings:
# This section only defines Peano-specific flags. It does not
# set compiler specific stuff.
#
build = ARGUMENTS.get('build', 'debug')   # Read command line parameter
if build == 'debug':
   cppdefines.append('MDDebug')
   cppdefines.append('TarchDebug')
   cppdefines.append('MDError')
   
elif build == 'release':
   pass
else:
   print "ERROR: build must be 'debug' or 'release'!"
   sys.exit(1)
   
##### Determine Machine
#
machine = ARGUMENTS.get('machine', 'generic') # Read command line parameter
mpiLibrary='mpi'
mpiIncludePath='/usr/lib/openmpi/include'
mpiLibraryPath='/usr/lib/openmpi/lib'
pthreadLibrary='pthread'

if machine== 'mac-cluster':
   mpiLibrary=''
   mpiIncludePath=''
   mpiLibraryPath=''

##### Determine MPI-Parallelization
#
parallel = ARGUMENTS.get('parallel', 'parallel_no') # Read command line parameter
if parallel == 'yes' or parallel == 'parallel_yes':
   cppdefines.append('MDParallel')
   cppdefines.append('TarchParallel')
   if machine != 'mac-cluster':
     cpppath.append(mpiIncludePath)
     libpath.append(mpiLibraryPath)
     libs.append('pthread')
   ccflags.append('-DMPICH_IGNORE_CXX_SEEK')
elif parallel == 'no' or parallel == 'parallel_no':
   pass
else:
   print "ERROR: parallel must be = 'yes', 'parallel_yes', 'no' or 'parallel_no'!"
   sys.exit(1)
   
##### Determine GProf usage
# 
gprof = ARGUMENTS.get('gprof', 'gprof_no') # Read command line parameter
if gprof == 'no' or gprof == 'gprof_no':
   pass
elif gprof == 'yes' or gprof == 'gprof_yes':
   ccflags.append('-p')
   ccflags.append('-pg')
   linkerflags.append('-p')
   linkerflags.append('-pg')
else:
   print "ERROR: gprof must be = 'yes', 'gprof_yes', 'no' or 'gprof_no'!"
   sys.exit(1)
   
##### Switch Compiler
#
compiler = ARGUMENTS.get('compiler', 'gcc') # Read command line parameter
if compiler == 'gcc':
   if(parallel=='parallel_no' or parallel=='no'):
     cxx = 'g++'
   else:
     cxx = 'mpicxx'

   ccflags.append('-Wall')
   ccflags.append('-Werror')
   ccflags.append('-pedantic')
   ccflags.append('-pedantic-errors')
   ccflags.append('-Wstrict-aliasing')
   ccflags.append('-fstrict-aliasing')
   ccflags.append('-Wno-long-long')
   ccflags.append('-Wno-unknown-pragmas')
   ccflags.append('-Wconversion')
   ccflags.append('-Wno-non-virtual-dtor')
   if build == 'debug':
      ccflags.append('-g3')
      ccflags.append('-O0')
   elif build == 'release':
      ccflags.append('-O3') 
elif compiler == 'icc':
   if parallel == 'yes':
     cxx = 'mpiicpc'
   else:
     cxx = 'icpc'

   ccflags.append('-fstrict-aliasing')
   ccflags.append('-qpack_semantic=gnu')
   ccflags.append('-ipo')
   linkerflags.append('-ipo')
   if build == 'debug':
      ccflags.append('-O0')
   elif build == 'release':
      ccflags.append('-fast')
      ccflags.append('-w')
      ccflags.append('-Werror-all')
      ccflags.append('-align')
      ccflags.append('-ansi-alias')
else:
   print "ERROR: compiler must be = 'gcc' or 'icc'!"
   sys.exit(1)
      
##### Determine build path
#
target = ARGUMENTS.get('target', '')  # Read command line parameter
build_offset = ARGUMENTS.get('buildoffset', 'build')
buildpath = build_offset  + '/' + str(target) + '/' + str(build) + '/dim' + str(dim) + '/' 
if parallel == 'yes' or parallel == 'parallel_yes':
   buildpath = buildpath + 'parallel_yes/'
else:
   buildpath = buildpath + 'parallel_no/'
if compiler == 'icc':
   buildpath = buildpath + 'icc/'
if compiler == 'gcc':
   buildpath = buildpath + 'gcc/'
if gprof == 'yes' or gprof == 'gprof_yes':
   buildpath = buildpath + 'gprof/'
else:
   buildpath = buildpath + 'gprof_no/'
   
##### Print options used to build
#
print "Target: " + target
print "Options: build = " + str(build) + ", dim = " + str(dim) + ", build-offset = " + str(build_offset) + ", parallel = " + str(parallel) + ", gprof = " + str(gprof) + ", compiler = " + str(compiler)
print "Buildpath: " + buildpath


VariantDir (buildpath, './', duplicate=0)  # Change build directory
libpath.append(buildpath)


#define general sources of molecular dynamics
sourcesMolecularDynamics = [
  Glob('simplemd/BoundaryTreatment.cpp'),
  Glob('simplemd/MolecularDynamicsSimulation.cpp'),
  Glob('simplemd/ProfilePlotter.cpp'),
  Glob('simplemd/cell-mappings/*.cpp'),
  Glob('simplemd/configurations/*.cpp'),
  Glob('simplemd/molecule-mappings/*.cpp'),
  Glob('simplemd/services/*.cpp'),
  Glob('tarch/tinyxml2/tinyxml2.cpp'),
  Glob('tarch/utils/RandomNumberService.cpp')
]



##### Setup construction environment:
#
env = Environment ( 
   CPPDEFINES = cppdefines,
   LIBPATH    = libpath,
   LIBS       = libs, 
   CPPPATH    = cpppath,
   CCFLAGS    = ccflags,
   LINKFLAGS  = linkerflags,
   CXX        = cxx,
   ENV        = os.environ
   )

if target=='simplemd':
  env.Program (
    target = buildpath + 'simplemd',
    source = [
      Glob('simplemd/main.cpp'),
      sourcesMolecularDynamics
    ]
  )
elif target=='multi-simplemd':
  env.Program (
    target = buildpath + 'multi-simplemd',
    source = [
      Glob('simplemd/main_multi.cpp'),
      sourcesMolecularDynamics
    ]
  )
elif target=='libsimplemd':
  env.Library(target= buildpath+'libsimplemd', source = [sourcesMolecularDynamics])
else:
  print "ERROR: Target", target, "is invalid!"
  printTargets ()
  sys.exit(1)
