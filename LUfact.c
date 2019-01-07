/* Erik Safford
   Program 2 - LU Factorization Engine
   CS 330
   October 2018 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "LUfact.h"

double **createMatrix(int N) { //Given N, creates an NxN matrix
  double **M = (double **) malloc(N*sizeof(double*));
  for (int i = 0; i < N; i++)
    M[i] = (double*) malloc(N*sizeof(double));
  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
      M[i][j] = (i == j) ? 1.0 : 0.0;
  return M;
}

void destroyMatrix(int N, double **M) { //Given pointer to 2d array (matrix) and N, frees allocated memory
  for (int i = 0; i < N; i++)
    free(M[i]);
  free(M);
}

int LUfactorize(double **A, int N, short *P) { 

    int i, j, k, imax;
    double maxA, *ptr, absA;

    for (i = 0; i <= N; i++)
        P[i] = i; //Unit permutation matrix, P[N] initialized with N

    for (i = 0; i < N; i++) {
        maxA = 0.0;
        imax = i;

        for (k = i; k < N; k++)
            if ((absA = fabs(A[k][i])) > maxA) {
                maxA = absA;
                imax = k;
            }

        if (imax != i) {
            //pivoting P
            j = P[i];
            P[i] = P[imax];
            P[imax] = j;

            //pivoting rows of A
            ptr = A[i];
            A[i] = A[imax];
            A[imax] = ptr;

            //counting pivots starting from N (for determinant)
            P[N]++;
        }

        for (j = i + 1; j < N; j++) {
            A[j][i] /= A[i][i];

            for (k = i + 1; k < N; k++)
                A[j][k] -= A[j][i] * A[i][k];
        }
    }

    return 1;  //factorization done
}

double LUdeterminant(double **A, int N, short *P) { //If det(A) == 0, A is singular 

    double det = A[0][0];

    for (int i = 1; i < N; i++)
        det *= A[i][i];

    if ((P[N] - N) % 2 == 0)
        return det;
    else
        return -det;
}

//Preforms factorization, and allocates, fills, and returns a LUfact object
LUfact *LUfactor(int N, const double **A) {
  LUfact *LU = (LUfact*) malloc(sizeof(LUfact));
  LU->N = N;
  LU->LU = createMatrix(N);
  LU->mutate = (short *) malloc(N*sizeof(short));

  // Clone A into LU
  double **A_ = LU->LU;
  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
      A_[i][j] = A[i][j];

  for (int i = 0; i < N; i++)
    LU->mutate[i] = (short) i;
  
  // actual factorizing happens here
  LUfactorize(A_,LU->N,LU->mutate);
  
  int check = LUdeterminant(A_,LU->N,LU->mutate);
  if (check == 0) //If det(A) == 0
    return NULL; //matrix is singular, return NULL

  return LU;
}

void LUdestroy(LUfact *fact) { //Deallocates the data allocated in LUfactor()
  //free fact->LU matrix
  destroyMatrix(fact->N,fact->LU);
  //free fact->mutate array
  free(fact->mutate);
  //free LUfact *fact
  free(fact);
}

void LUsolve(LUfact *fact, const double *b, double *x) { //Solves the system Ax=b for x
  double **A = fact->LU;
  short *P = fact->mutate;
  int N = fact->N;

  for (int i = 0; i < N; i++) {
    x[i] = b[P[i]];

    for (int k = 0; k < i; k++) {
      x[i] -= A[i][k] * x[k];
    }
  }

  for (int i = N - 1; i >= 0; i--) {
    for (int k = i + 1; k < N; k++) {
      x[i] -= A[i][k] * x[k];
    }
    x[i] = x[i] / A[i][i];
  }
}
