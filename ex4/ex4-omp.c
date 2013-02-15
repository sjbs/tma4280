#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "common.h"

//current code doesn't generate vector, retard. fix it.

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
	
	Matrix v = createMatrix(1,N);
	
	#pragma omp parallel for schedule(static) reduction(+:v)
	for (long int i=1;i<N+1;++i) {
		double temp=1.0/(i*i);
		v->data[i-1][1] = temp;
		sum += temp;
	}
	
	//exact value of pi
	double pi=(4.0*atan(1.0));
	
	printf("sum: %f\n", sum);
	printf("difference: %1.16f\n", pi*pi/6.0-sum);
	
}
	
	

  
	
	
