#!/bin/bash

#PBS -N partb_nn_1
#PBS -lnodes=1:ppn=12:default
#PBS -lwalltime=03:00:00
#PBS -lpmem=2000MB
#PBS -A freecycle
#PBS -q optimist
#PBS -j oe

cd /home/spence/tmacode-checkout/ex6/release

module load intelcomp
module load openmpi/1.4.3-intel
module load cmake
KMP_AFFINITY="granularity=fine, compact"

cat $PBS_NODEFILE
echo "   "
echo "nn=1"
for npern in 2 3 6 12
do
	for t in 1 2 3 6 12
	do
		let "temp=$npern*$t"
		if [ "$temp" -le "12" ]
		then
			for n in 64 256 512 2048 4096 8192 16384 
			do
				echo "---------"
				echo "n=$n"
				export OMP_NUM_THREADS=$t
				mpirun -npernode $npern ./mpi-poisson $n
			done
		fi
	done
done

echo "f=initfunc:"
export OMP_NUM_THREADS=1
mpirun -npernode 12 ./mpi-poisson-initf 8192






