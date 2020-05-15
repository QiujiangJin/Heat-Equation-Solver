/*************************************************************************************
 * This function solves the lienar system using Jacobi method or Gauss Seidel method
 * **********************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <grvy.h>
#include "laplace.h"
#include "function.h"

void solve_linear_system(Laplace* solver) {

    grvy_timer_begin(__func__);

    if(solver->output_mode == 1) {
	printf("[debug]: solve system - function begin\n\n\n");
    }

    printf("** Solving linear system...\n");

    int N = solver->nx;
    if(solver->dimension == 2) {
	N = (N-1)*(N-1) + 1;
    }
    int iter = 0;
    int max_iter = solver->max_iter;
    double eps = solver->eps;

    // this vector is to store the values of the vector form last iteration    
    double *x_vec = (double* )malloc((N-1)*sizeof(double));
    double err = 1.0;
    
    if(strcmp(solver->iter_method, "Jacobi") == 0) {
	// case of Jacobi method
	while(iter < max_iter && err >= eps) {
	    if(solver->output_mode == 1) {
		printf("current convergence error = %.12f\n", err);
	    }
	    memcpy(x_vec, solver->sol_vec, (N-1)*sizeof(double));
	    for(int i = 0; i < N-1; i++) {
		double s = 0.0;
		int diag = 0;
		for(int j = 0; j < solver->length[i]; j++) {
		    int index = solver->matrix_index[i][j];
		    if(index != i) {
			s += solver->matrix_value[i][j]*x_vec[index];
		    } else {
			diag = j;
		    }
		}
		solver->sol_vec[i] = (solver->rht_vec[i] - s)/solver->matrix_value[i][diag];
	    }	
	    iter += 1;
            err = error_l2(solver, x_vec);
	}

    } else {
	// case of Gauss Seidel method
	while(iter < max_iter && err >= eps) {                  
	    if(solver->output_mode == 1) {
                printf("current convergence error = %.12f\n", err);
            }
	    memcpy(x_vec, solver->sol_vec, (N-1)*sizeof(double));
	    for(int i = 0; i < N-1; i++) {
                double s = 0.0;
                int diag = 0;
                for(int j = 0; j < solver->length[i]; j++) {
                    int index = solver->matrix_index[i][j];
                    if(index < i) {
                        s += solver->matrix_value[i][j]*solver->sol_vec[index];
                    } else if(index == i){
                        diag = j;
                    } else {
			s += solver->matrix_value[i][j]*x_vec[index];
		    }
                }
                solver->sol_vec[i] = (solver->rht_vec[i] - s)/solver->matrix_value[i][diag];
            }
            iter += 1;
            err = error_l2(solver, x_vec);
        }
    
    }

    printf("--> Terminated at iteration: %i\n", iter);
    printf("--> The error norm: %.12f\n", error_l2(solver, x_vec));
    printf("--> The residual norm: %.12f\n\n\n", res_norm(solver));

    free(x_vec);
    
    if(solver->output_mode == 1) {
        printf("[debug]: solve system - function end\n\n\n");
    }

    grvy_timer_end(__func__);

}
