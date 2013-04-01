#!/bin/bash

#PBS -N pnntest
#PBS -lnodes=1:ppn=12:default
#PBS -lwalltime=00:05:00
#PBS -lpmem=2000MB
#PBS -A freecycle
#PBS -q optimist
#PBS -j oe

cd /home/spence/tmacode-checkout/ex6/release

module load intelcomp
module load openmpi/1.4.3-intel
module load cmake
KMP_AFFINITY="granularity=fine, compact"
echo "testrun mpi 12:"
export OMP_NUM_THREADS=1
mpirun -hostfile $PBS_NODEFILE ./mpi-poisson 4096

#PBS -lnodes=1:ppn=6:default
echo "testrun mpi 6:"
export OMP_NUM_THREADS=1
mpirun -hostfile $PBS_NODEFILE ./mpi-poisson 4096

#PBS -lnodes=1:ppn=6:default
echo "testrun omp 2:"
export OMP_NUM_THREADS=2
mpirun -hostfile $PBS_NODEFILE ./mpi-poisson 4096
