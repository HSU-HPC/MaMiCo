#!/bin/bash
#SBATCH --job-name=couette_middle_2_oscil_0_5_u_wall # specifies a user-defined job name

#SBATCH --nodes=1             # number of compute nodes to be used
#SBATCH --ntasks=1           # number of MPI processes
#SBATCH --partition=small #
#SBATCH --cpus-per-task=1    # number of cores per process
#SBATCH --time=12:00:00
#SBATCH --output=/beegfs/project/MaMiCo/mamico-ml/MaMiCo/results/couette_oscil/2_oscil/u_wall_0_5/1_middle/couette_middle_2_oscil_0_5_u_wall%j.log  # log file which will contain all output

# commands to be executed
source /beegfs/project/MaMiCo/setup-env.sh
cd /beegfs/project/MaMiCo/mamico-ml/MaMiCo/results/couette_oscil/2_oscil/u_wall_0_5/1_middle/

./couette
