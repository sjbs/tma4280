#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "common.h"		


int main(int argc, char** argv)
{
	// if (argc < 2) {
	// 	printf("Need n, number of vector elements. Tard.\n");
	// 	return 1;
	// }
	
	//int N=atoi(argv[1]);
	
	int N = pow(2,14);
	
	//exact value of pi
	double pi=(4.0*atan(1.0));
	// if (N<=0) {
	// 	printf("Need a positive vector length, mormon.\n");
	// 	return 1;
	// }
	
	double sum=0;
	int k=4;
	
	Vector v = createVector(N);
	
	for (int i=0;i<N;++i) {
		v->data[i] = 1.0/((double)(i+1)*(double)(i+1));	
		sum += v->data[i];
		//if ( ((i+1 != 0) && !(i+1 & (i+1 - 1)))&&i+1>=16){
		if (i+1==pow(2,k)){
			printf("difference at i=2^%2.0f: %1.16f\n", log(i)/log(2), pi*pi/6.0-sum);
			k++;
		}
	}
	
	
	
//	printf("sum: %f\n", sum);
	
	
}
	
	

  
	
	
