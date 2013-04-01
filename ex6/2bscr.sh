#!/bin/bash

#PBS -N t6_1-2
#PBS -lnodes=1:ppn=2:default
#PBS -lwalltime=00:10:00
#PBS -lpmem=2000MB
#PBS -A freecycle
#PBS -q optimist
#PBS -j oe

cd /home/spence/tmacode-checkout/ex6/release

module load intelcomp
module load openmpi/1.4.3-intel
module load cmake
KMP_AFFINITY="granularity=fine, compact"
echo "nn=1:"
export OMP_NUM_THREADS=6
for n in 64 256 512 2048 4096 8192 16384
do
	echo "n=$n"
	mpirun -hostfile $PBS_NODEFILE ./mpi-poisson $n
done




