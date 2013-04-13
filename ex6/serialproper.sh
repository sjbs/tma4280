#!/bin/bash

#PBS -N serialproper
#PBS -lnodes=1:ppn=12:default
#PBS -lwalltime=00:40:00
#PBS -lpmem=2000MB
#PBS -A freecycle
#PBS -q optimist
#PBS -j oe

cd /home/spence/tmacode-checkout/ex6/release

module load intelcomp
module load openmpi/1.4.3-intel
module load cmake
KMP_AFFINITY="granularity=fine, compact"
export OMP_NUM_THREADS=1
for n in 64 256 512 2048 4096 8192 16384 
do
	echo "n=$n"
	echo "serial:"
	./poisson $n
done






