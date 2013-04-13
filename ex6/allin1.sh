#!/bin/bash

#PBS -N allinonetest
#PBS -lnodes=3:ppn=12:default
#PBS -lwalltime=00:5:00
#PBS -lpmem=2000MB
#PBS -A freecycle
#PBS -q optimist
#PBS -j oe

cd /home/spence/tmacode-checkout/ex6/release

module load intelcomp
module load openmpi/1.4.3-intel
module load cmake
KMP_AFFINITY="granularity=fine, compact"

echo "serial:"
./poisson 1024

for npers in 1 2 3 6
do
	for t in 1 2 3 6
	do
		let "temp=$npers*$t"
		if [ "$temp" -le "6" ]
		then
			for n in 64 
			do
				echo "---------"
				echo "n=$n"
				export OMP_NUM_THREADS=$t
				echo "nn=1"
				mpirun -H 0 -npersocket %npers ./mpi-poisson $n
				echo "nn=2"
				mpirun -H 0,1 -npersocket %npers ./mpi-poisson $n
				echo "nn=3"
				mpirun -H 0,1,2 -npersocket %npers ./mpi-poisson $n
				echo "-------"
			done
		fi
	done
done

echo "f=initfunc:"
export OMP_NUM_THREADS=1
mpirun -H 0,1,2 -npersocket 6 ./mpi-poisson-initf 8192



echo "64 256 512 2048 4096 8192 16384"




