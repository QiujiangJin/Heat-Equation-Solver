/******************************************************************************************
 * This function initializes all the other variables in the solver like matrix and vector
 * ***************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <grvy.h>
#include "laplace.h"
#include "function.h"

void initialize(Laplace* solver) {

    grvy_timer_begin(__func__);

    // print the information of masa solution
    if(solver->dimension == 1) {
    	printf("MASA :: Solution has 2 variables.\n");
    	printf("*-------------------------------------*\n");
    	printf("A_x is set to: 10.0\n");
    	printf("k_0 is set to: 1.0\n");
    	printf("*-------------------------------------*\n\n\n");
    } else {
    	printf("MASA :: Solution has 3 variables.\n");
        printf("*-------------------------------------*\n");
        printf("A_x is set to: 10.0\n");
        printf("B_y is set to: 10.0\n");
        printf("k_0 is set to: 1.0\n");
        printf("*-------------------------------------*\n\n\n");
    }

    if(solver->output_mode == 1) {
    	printf("debug modes is launched\n");
    	printf("[debug]: initialization - function begin\n\n\n");
    }

    printf("** Initializing data structures...\n\n\n");

    solver->h = (solver->xmax - solver->xmin)/solver->nx;
    int N = solver->nx - 1;
    if(solver->dimension == 2) {
	N = N*N;
    }

    // dynamically allocate memory for the pointers
    solver->matrix_value = (double **)malloc(N*sizeof(double *));
    solver->matrix_index = (int **)malloc(N*sizeof(int *));
    solver->length = (int *)calloc(N, sizeof(int));
    solver->rht_vec = (double *)malloc(N*sizeof(double));
    solver->sol_vec = (double *)calloc(N, sizeof(double));

    if(solver->output_mode == 1) {
        printf("[debug]: initialization - function end\n\n\n");
    }

    grvy_timer_end(__func__);

}
