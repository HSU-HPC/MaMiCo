#!/bin/bash
#SBATCH --job-name=kvs_init_24000_NE # specifies a user-defined job name

#SBATCH --nodes=1             # number of compute nodes to be used
#SBATCH --ntasks=1           # number of MPI processes
#SBATCH --partition=small_gpu #
#SBATCH --gpus=1
#SBATCH --cpus-per-task=1    # number of cores per process
#SBATCH --time=06:00:00
#SBATCH --output=/beegfs/project/MaMiCo/mamico-ml/MaMiCo/results/kvs/2_init_timesteps_24000/0_NE/kvs_init_24000_NE_%j.log  # log file which will contain all output

# commands to be executed
source /beegfs/project/MaMiCo/setup-env.sh
cd /beegfs/project/MaMiCo/mamico-ml/MaMiCo/results/kvs/2_init_timesteps_24000/0_NE/

python kvstest.py
