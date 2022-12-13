#!/bin/bash
#SBATCH -o /home/hpc/pr63so/lu28caq2/workspace/moleculardynamics/coupling/tests/testmardyn.out
#SBATCH -D /home/hpc/pr63so/lu28caq2/
#SBATCH -J lbmd_30_mardyn
#SBATCH --get-user-env
#SBATCH --partition=snb
#SBATCH --nodelist=mac-snb17
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=2 # one thread per process
#SBATCH --mail-type=end
#SBATCH --mail-user=neumanph@in.tum.de
#SBATCH --export=NONE
#SBATCH --time=24:00:00
source /etc/profile.d/modules.sh

# -ppn = processes per node
# -n   = total number of processes
cd /home/hpc/pr63so/lu28caq2/workspace/moleculardynamics/coupling/tests
date

./test_mardyn

date
