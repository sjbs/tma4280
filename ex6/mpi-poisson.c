/*
  C-program to solve the two-dimensional Poisson equation on 
  a unit square using one-dimensional eigenvalue decompositions
  and fast sine transforms

  einar m. ronquist
  ntnu, october 2000
  revised, october 2001

  revised to parallel version by Silas Spence
  March 2013
*/

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>

typedef double Real;

/* function prototypes */
Real *createRealArray (int n);
int *createIntArray (int n);
Real **createReal2DArray (int m, int n);
void transpose (Real **bt, Real **b, int m, int mglob, Real *sendbuf, 
    Real *recbuf, int *sendcnt, int *sdispl,int rank);
void fst_(Real *v, int *n, Real *w, int *nn);
void fstinv_(Real *v, int *n, Real *w, int *nn);
void splitVector(int globLen, int size, int** len, int** displ);

int main(int argc, char **argv )
{
  int rank, size, mglob, *ofs, *cols, *sendcnt, *sdispl;
  
  MPI_Init (&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &size);
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);


  Real *diag, **b, **bt, *z, *sendbuf, *recbuf;
  Real pi, h, umax;
  int i, j, n, m, nn;

  /* the total number of grid points in each spatial direction is (n+1) */
  /* the total number of degrees-of-freedom in each spatial direction is (n-1) */
  /* this version requires n to be a power of 2 */


  if( argc < 2 ) {
    if (rank==0){
    printf("need a problem size\n");
    }
    MPI_Finalize();
	  return(0);
  }

  n  = atoi(argv[1]);
  mglob = n-1;
  nn = 4*n;

  splitVector(n-1, size, &cols, &ofs);

  m=cols[rank];
  int localdof = m*mglob;
  printf("localdof[%i]=%i\n",rank,localdof);
  //nn=4*nlocal;

  //printf("%i\n",ofs[rank]);

  diag = createRealArray (mglob);
  b    = createReal2DArray (mglob,m);
  bt   = createReal2DArray (mglob,m);
  z    = createRealArray (nn);
  sendbuf = createRealArray (localdof);
  recbuf  = createRealArray (localdof);
  sendcnt = createIntArray (size);
  sdispl  = createIntArray (size);

  for (i=0; i<size; i++){
    sendcnt[i]=cols[rank]*cols[i];
    if (rank==1) printf("sendcnt=%i  ",sendcnt[i]);
  }
  if (rank==1) printf("\n");

  sdispl[0]=0;
  for (i=1; i<size; i++){
    sdispl[i]=sdispl[i-1]+sendcnt[i-1];
    if (rank==1) printf("sdispl=%i  ",sdispl[i]);
  }
  if (rank==1) printf("\n");
  
  h    = 1./(Real)n;
  pi   = 4.*atan(1.);

  for (i=0; i < mglob; i++) {
    diag[i] = 2.*(1.-cos((i+1)*pi/(Real)n));
  }
  for (j=0; j < mglob; j++) {
    for (i=0; i < m; i++) {
      b[j][i] = (i+1)*(j+2)-rank;
      if (rank==1) printf("b[%i][%i]=%f\t",j,i,b[j][i]);
    }
    if (rank==1) printf("\n");
  }
  
  for (j=0; j < m; j++) {
    //fst_(b[j], &n, z, &nn);
  }

  transpose (bt,b,m,mglob,sendbuf,recbuf,sendcnt,sdispl,rank);

  // for (i=0; i < m; i++) {
  //  fstinv_(bt[i], &n, z, &nn);
  // }
  
  // for (j=0; j < mglob; j++) {
  //  for (i=0; i < m; i++) {
  //    bt[j][i] = bt[j][i]/(diag[i]+diag[j]);
  //   }
  // }
  
  // for (i=0; i < m; i++) {
  //   fst_(bt[i], &n, z, &nn);
  // }

  // transpose (bt,b,m,mglob,sendbuf,recbuf,sendcnt,sdispl);

  // for (j=0; j < m; j++) {
  //   fstinv_(b[j], &n, z, &nn);
  // }

  // umax = 0.0;
  // for (j=0; j < mglob; j++) {
  //   for (i=0; i < m; i++) {
  //     if (b[j][i] > umax) umax = b[j][i];
  //   }
  // }
  // printf (" umax = %e \n",umax);

  MPI_Finalize();
  return 0;
}

void transpose (Real **bt, Real **b, int m, int mglob, Real *sendbuf,
  Real *recbuf,int *sendcnt, int *sdispl,int rank)
{
  for (int j=0; j < mglob; j++) {
    for (int i=0; i < m; i++) {
      int ind=j*m+i;
      sendbuf[ind]=b[j][i];
      if (rank==1) printf("%2.1f  ",sendbuf[ind]);
    }
  }
  if (rank==1) printf("\n\n");

  MPI_Alltoallv(sendbuf, sendcnt, sdispl,MPI_DOUBLE,
     recbuf, sendcnt, sdispl,MPI_DOUBLE,MPI_COMM_WORLD);
  
  for (int i=0;i<m*mglob;i++) {
    if (rank==1) printf("%2.1f  ",recbuf[i]);
  }
  if (rank==1) printf("\n\n");

  for (int j=0; j < mglob; j++) {
    for (int i=0; i < m; i++) {
      int ind=j*m+i;
      bt[j][i]=recbuf[ind];
      if (rank==1) printf("%2.1f  ",bt[j][i]);
    }
    if (rank==1) printf("\n");
  }

}

Real *createRealArray (int n)
{
  Real *a;
  int i;
  a = (Real *)malloc(n*sizeof(Real));
  for (i=0; i < n; i++) {
    a[i] = 0.0;
  }
  return (a);
}

int *createIntArray (int n)
{
  int *a;
  int i;
  a = (int *)malloc(n*sizeof(int));
  for (i=0; i < n; i++) {
    a[i] = 0;
  }
  return (a);
}

Real **createReal2DArray (int n1, int n2)
{
  int i, n;
  Real **a;
  a    = (Real **)malloc(n1   *sizeof(Real *));
  a[0] = (Real  *)malloc(n1*n2*sizeof(Real));
  for (i=1; i < n1; i++) {
    a[i] = a[i-1] + n2;
  }
  n = n1*n2;
  memset(a[0],0,n*sizeof(Real));
  return (a);
}

void splitVector(int globLen, int size, int** len, int** displ)
{
  *len = calloc(size,sizeof(int));
  *displ = calloc(size,sizeof(int));
  for (int i=0;i<size;++i) {
    (*len)[i] = globLen/size;
    if (globLen % size && i >= (size - globLen % size))
      (*len)[i]++;
    if (i < size-1)
      (*displ)[i+1] = (*displ)[i]+(*len)[i];
  }
}

  // for (i=0; i < m; i++) {
  //   diag[i] = 2.*(1.-cos((i+1+ofs[rank])*pi/(Real)n));
  //   //diag[i]=i+ofs[rank];
  //   //printf("diag[%i]=%i\n",rank,i+ofs[rank]); 
  //   printf("%2.4f\n",diag[i]);
  // }