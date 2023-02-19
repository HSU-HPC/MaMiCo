#!/bin/bash
#SBATCH --job-name=kvs # specifies a user-defined job name

#SBATCH --nodes=1             # number of compute nodes to be used
#SBATCH --ntasks=1           # number of MPI processes
#SBATCH --partition=small_gpu #
#SBATCH --gpus=1
#SBATCH --cpus-per-task=1    # number of cores per process
#SBATCH --time=12:00:00
#SBATCH --output=kvs_%j.log  # log file which will contain all output

# commands to be executed
source /beegfs/project/MaMiCo/setup-env.sh
cd /beegfs/project/MaMiCo/mamico-ml/results/kvs/

python kvstest.py
