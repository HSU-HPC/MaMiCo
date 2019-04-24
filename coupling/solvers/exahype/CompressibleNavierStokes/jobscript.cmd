#!/bin/bash
#@ job_type     = parallel
#@ class        = micro
#@ node         = 1
#@ tasks_per_node = 1
#@ island_count = 1
#@ network.MPI = sn_all,not_shared,us 
#@ energy_policy_tag = NONE
#@ hold = user
#@ wall_clock_limit = 48:00:00
#@ job_name = NSBenchmark
#@ node_usage = not_shared
#@ initialdir = /gss/scratch/pr83no/ga24dib2/exahype_bench/logs/working_dir/{dir}
#@ error  =  /gss/scratch/pr83no/ga24dib2/exahype_bench/logs/job.$(schedd_host).$(jobid).err
#@ output =  /gss/scratch/pr83no/ga24dib2/exahype_bench/logs/job.$(schedd_host).$(jobid).out
#@ notification=always
#@ notify_user=lukas.krenz@tum.de
#@ queue

. /etc/profile
. /etc/profile.d/modules.sh
module switch intel/18.0
module switch tbb/2018
module switch gcc/5

export OMP_NUM_THREADS=28
export MP_TASK_AFFINITY=core:28

poe /home/hpc/pr83no/ga24dib2/ExaHyPE-Engine/ApplicationExamples/CompressibleNavierStokes/ExaHyPE-NavierStokes {config_file} 
