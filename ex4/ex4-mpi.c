#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "common.h"

int main(int argc, char **argv)
{
	if (argc<=1){
	  printf("touch my wiener\n\n");
	}
	int testn=atoi(argv[1]);
	if (testn<=0){
	  printf("testn less than zero.\n");
	  MPI_Finalize;
	  return 1;
	}
	double mysum, sum, x, pi, t1, t2;
	int n, myid, nproc, i;
	MPI_Init (&argc, &argv);
	MPI_Comm_size (MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank (MPI_COMM_WORLD, &myid);

	//printf("comm_size is %i \n", nproc);
	//printf("comm_rank is %i \n", myid);

	//ask about input - arg vs scanf
	if (myid == 0) {
		printf (" Enter vector length:\n");
		scanf ("%i",&n);
		if (n<=0){
			printf("need a positive vector length, idiot.");
			MPI_Finalize();
			return 1;
		}
		t1 = MPI_Wtime();
	}
	MPI_Bcast (&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
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
		printf("sum: %f\n", sum);
		printf("difference: %1.16f\n", pi*pi/6.0-sum);
		printf("Wall Time is %1.4f\n",t2-t1);
	}
	
	MPI_Finalize();
	return(0);
}
