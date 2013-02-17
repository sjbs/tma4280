#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "common.h"
#include "ctype.h"


int main(int argc, char** argv)
{
	// if (argc < 2) {
	// 	printf("Need n, number of vector elements. Tard.\n");
	// 	return 1;
	// }
	
	//int N=atoi(argv[1]);
	
	int N = pow(2,14);
	printf("\n %i \n",N);

	double pi=(4.0*atan(1.0));
	// if (N<=0) {
	// 	printf("Need a positive vector length, mormon.\n");
	// 	return 1;
	// }
	
	double sum=0;
	
	Vector v = createVector(N);
	
	for (int i=0;i<N;++i) {
		double temp=1.0/((i+1)*(i+1));
		v->data[i] = temp;	
		sum += temp;
		if (i<=18&&i>=16){
			printf("difference at i=2^%i: %1.16f\n", i, pi*pi/6.0-sum);
		}
	}
	
	//exact value of pi
	
	
	printf("sum: %f\n", sum);
	
	
}
	
	

  
	
	
