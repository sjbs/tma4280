#!/bin/bash

#PBS -N c-t9
#PBS -lnodes=4:ppn=1:default
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
echo "nn=4:"
export OMP_NUM_THREADS=9
for n in 16384
do
	echo "n=$n"
	mpirun -hostfile $PBS_NODEFILE ./mpi-poisson $n
done




