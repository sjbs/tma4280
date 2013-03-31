#!/bin/bash

#PBS -N hybridtest
#PBS -lnodes=2:ppn=2:default
#PBS -lwalltime=00:10:00
#PBS -lwalltime=00:01:00
#PBS -lpmem=2000MB
#PBS -A freecycle
#PBS -q optimist
#PBS -j oe

cd /home/spence/tmacode-checkout/ex6/release

module load intelcomp
module load openmpi/1.4.3-intel
module load cmake
KMP_AFFINITY="granularity=fine, compact"
echo "testrun 2 nodes hybrid 1:"
export OMP_NUM_THREADS=1
mpirun --hostfile $PBS_NODEFILE ./mpi-poisson 4096
echo "testrun 2 nodes hybrid 3:"
export OMP_NUM_THREADS=3
mpirun --hostfile $PBS_NODEFILE ./mpi-poisson 4096
echo "testrun 2 nodes hybrid 6:"
export OMP_NUM_THREADS=6
mpirun --hostfile $PBS_NODEFILE ./mpi-poisson 4096



