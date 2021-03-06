#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "common.h"

//remember to implement 2^k!

int main(int argc, char** argv)
{
	if (argc < 2) {
		printf("Need n, number of vector elements. Tard.\n");
		return 1;
	}
	
	int N=atoi(argv[1]);
	
	if (N<=0) {
		printf("Need a positive vector length, mormon.\n");
		return 1;
	}
	
	double sum=0;
	
	Vector v = createVector(N);

	double t1=WallTime();
	
	#pragma omp parallel for schedule(static) reduction(+:sum)
	for (int i=0;i<N;++i) {
		double temp=1.0/((double)(i+1)*(double)(i+1));
		v->data[i] = temp;	
		sum += temp;
	}

	double t2=WallTime();

	//exact value of pi
	double pi=(4.0*atan(1.0));
	
	printf("sum: %f\n", sum);
	printf("difference: %1.16f\n", pi*pi/6.0-sum);
	printf("Wall Time is %1.4f\n",t2-t1);
}
	
	

  
	
	
