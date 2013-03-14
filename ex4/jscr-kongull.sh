#!/bin/bash

#PBS -N test
#PBS -lnodes=1:ppn=12:default
#PBS -lwalltime=00:01:00
#PBS -lpmem=2000MB
#PBS -A freecycle
#PBS -q optimist
#PBS -j oe

cd /home/spence/tmacode-checkout/ex4/release

module load intelcomp
module load openmpi/1.4.3-intel
module load cmake
KMP_AFFINITY="granularity=fine, compact"
echo "serial:"
./serial2k
echo "omp 4 :"
OMP_NUM_THREADS=4 ./omp2k
echo "mpi 4:"
mpirun -np 4 mpi-splitv
echo "omp 16 :"
OMP_NUM_THREADS=16 ./omp2k
echo "mpi 16:"
mpirun -np 16 mpi-splitv

