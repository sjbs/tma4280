#!/bin/bash

#PBS -N serial
#PBS -lnodes=1:ppn=1:default
#PBS -lwalltime=00:20:00
#PBS -lpmem=20000MB
#PBS -A freecycle
#PBS -q optimist
#PBS -j oe

cd /home/spence/tmacode-checkout/ex6/release

module load intelcomp
module load openmpi/1.4.3-intel
module load cmake
KMP_AFFINITY="granularity=fine, compact"
echo "testrun serial:"
export OMP_NUM_THREADS=1
for n in 64 256 512 2048 4096 8192 16384 32768
do
	echo "n=$n"
	mpirun -hostfile $PBS_NODEFILE ./mpi-poisson $n
done




