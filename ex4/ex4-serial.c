#include <stdlib.h>
#include <stdio.h>
#include <math.h>

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
	
	double sum=1;
	
	for (int i=2;i<N+1;++i) {
		sum += 1.0/(i*i);
	}
	
	double pi=(4.0*atan(1.0));
	
	printf("difference: %1.16f\n", pi*pi/6.0-sum);
//	printf("sum: %f\n", sum);
	
}
	
	

  
	
	
