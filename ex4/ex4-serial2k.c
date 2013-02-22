#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "common.h"		


int main(int argc, char** argv)
{
	//maximum size of vector - 2^14 for Exercise 4
	int maxN = pow(2,14);
	
	//value of sum as n->infinity
	double exactsum=pow((4.0*atan(1.0)),2)/6.0;
	
	//initialize some variables	
	double sum=0;
	int k=4;
	int nextN=pow(2,k);
	
	//make the vector
	Vector v = createVector(maxN);
	
	//fill the vector
	for (int i=0;i<maxN;++i) {
		v->data[i] = 1.0/((double)(i+1)*(double)(i+1));	
		
		//calculate the sum on the fly - saves a second for loop
		sum += v->data[i];
		
		//print difference at 2^k, k=4...14
		//for the price of a logical every i, only have to run through vector
		//once instead of length(k) times
		if (i+1==nextN) {
			printf("difference at i=2^%2i: %1.16f\n", k, exactsum-sum);
			k++;
			nextN=pow(2,k);
		}
	}
}
	
	

  
	
	
