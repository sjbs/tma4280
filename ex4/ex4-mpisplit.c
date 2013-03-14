#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "common.h"

int main(int argc, char** argv)
{
  int rank, size;
  init_app(argc, argv, &rank, &size);

  //check number of processors used is a power of two, otherwise getting the
  //the code for finding when the reduce happens gets fugly, and this is not
  //a programming course. in theory. but the the difference between theory
  //and practice is, in practice, larger than it is in theory. experiments
  //taught me that. so now i want to learn how to make pretty color plots and
  //never have anything to do with the real world ever again.
  if(size % 2){
    if(rank==0){
      printf("please choose a number of processors that is a multiple of 2, exiting\n");
    }
    close_app();
    return 0;
  }

  //maximum size of vector - 2^14 for Exercise 4
  int maxN = pow(2,14);
    
  //initialize some variables
  double mysum = 0;
  double sum=0;
  int k=4;
  int nextN=pow(2,k)/size;

  //split whole vector in parts, create subvector on each MPI process
  int *ofs, *cols;
  splitVector(maxN, size, &cols, &ofs);
  Vector v = createVector(cols[rank]);

  //fill vector
  for (int i=0;i<cols[rank];++i){
    //vector filling is interleaved - allows printing the difference
    //at each 2^k rather than recalculating entire vector each time
    v->data[i] = 1.0/pow((double)(i*size+1+rank),2);
    mysum += v->data[i];
  
	//print out difference at every 2^K n
    if (i+1==nextN){
      MPI_Reduce (&mysum, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      if (rank == 0) {
        double pi = 4.0 * atan(1.0);
        printf("difference at i=2^%2i: %1.16f\n", k, pi*pi/6.0-sum);
      }
      k++;
      nextN=pow(2,k)/size;
    }
  }

  //house cleaning. probably not needed (other than MPI_Finalize, but if you are 
  //going to cut and paste somebody else's code, you might as well go for broke.
  freeVector(v);
  free(cols);
  free(ofs);
  close_app();
  return 0;
}
