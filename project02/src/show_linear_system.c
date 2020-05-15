/*******************************************************************
 * This function prints the matrix and vector of the linear system
 * ****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "laplace.h"
#include "function.h"

void show_linear_system(Laplace* solver) {

    int N = solver->nx - 1;
    if(solver->dimension == 2) {
    	N = N*N;
    }

    double** A = (double **)malloc(N*sizeof(double *));
    for(int i = 0; i < N; i++) {
	A[i] = (double *)calloc(N, sizeof(double));
	for(int j = 0; j < solver->length[i]; j++) {
            A[i][solver->matrix_index[i][j]] = solver->matrix_value[i][j];
        }
    }
    double* b = (double *)calloc(N, sizeof(double));
    for(int i = 0; i < N; i++) {
        b[i] = solver->rht_vec[i];
    }

    // print matrix
    printf("A =\n\n");
    printf("[\n");
    for(int i = 0; i < N; i++) {
    	for(int j = 0; j < N; j++) {
	    printf("%3.0f ", A[i][j]);
	}
	printf("\n");
    }
    printf("];\n\n\n");

    // print right sided vector
    printf("b =\n\n");
    printf("[\n");
    for(int i = 0; i < N; i++) {
    	printf("%.2f\n", b[i]);
    }
    printf("];\n\n\n");

}
