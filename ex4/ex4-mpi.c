#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "common.h"

int main(int argc, char **argv)
{
	
	double mysum, sum, x, pi, t1, t2;
	int n, myid, nproc, i;
	MPI_Init (&argc, &argv);
	MPI_Comm_size (MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank (MPI_COMM_WORLD, &myid);

//	printf("comm_size is %i \n", nproc);
//	printf("comm_rank is %i \n", myid);

	if (argc<=1){
		if (myid==0) printf("Need vector length, n, as input argument. Exiting.\n");
		MPI_Finalize;
		return 1;
	}
	
	n=atoi(argv[1]);
	
	if (n<=0){
		if (myid==0) printf("Need a positive vector length. Shmuck.\n");
		MPI_Finalize;
		return 1;
	}

	
	if (myid == 0) t1 = MPI_Wtime();

//	MPI_Bcast (&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	sum = 0.0;
	for (i = myid; i <= n-1; i += nproc) {
		double temp=1.0/((double)(i+1)*(double)(i+1));
		//v->data[i] = temp;	
		mysum += temp;
	}

	MPI_Reduce (&mysum, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	if (myid == 0) {
		t2 = MPI_Wtime();
		pi = 4.0 * atan(1.0);
//		printf("sum: %f\n", sum);
		printf("difference: %1.16f\n", pi*pi/6.0-sum);
//		printf("Wall Time is %1.4f\n",t2-t1);
	}
	
	MPI_Finalize();
	return(0);
}

