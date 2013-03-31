#!/bin/bash

#PBS -N omptest
#PBS -lnodes=1:ppn=1:default
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
echo "testrun omp 1:"
export OMP_NUM_THREADS=1
mpirun -np 1 mpi-poisson 4096
echo "testrun omp 6:"
export OMP_NUM_THREADS=6
mpirun -np 1 mpi-poisson 4096
echo "testrun omp 12:"
export OMP_NUM_THREADS=12
mpirun -np 1 mpi-poisson 4096



