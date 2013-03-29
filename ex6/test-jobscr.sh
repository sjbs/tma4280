#!/bin/bash

#PBS -N test
#PBS -lnodes=1:ppn=12:default
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
echo "testrun mpi 1:"
export OMP_NUM_THREADS=1
mpirun -np 1 4096
echo "testrun mpi 12:"
mpirun -np 12 4096
echo "testrun omp 12:"
export OMP_NUM_THREADS=12
mpirun -np 1 4096
echo "testrun mixed p=2,t=6:"
export OMP_NUM_THREADS=6
mpirun -np 2 4096
echo "testrun mixed p=6,t=2:"
export OMP_NUM_THREADS=2
mpirun -np 6 4096



