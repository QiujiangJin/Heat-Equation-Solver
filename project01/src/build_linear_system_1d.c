/******************************************************************
 * This function builds the linear system in one dimensional case
 * ***************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <grvy.h>
#include <masa.h>
#include "laplace.h"
#include "function.h"

void build_linear_system_1d(Laplace* solver) {

    grvy_timer_begin(__func__);

    printf("** Building one dimensional linear system...\n");
    printf("--> Enforcing analytic Dirichlet BCs using MASA (1D)\n\n\n");

    if(solver->output_mode == 1) {
        printf("[debug]: build linear system - function begin\n\n\n");
    }

    int N = solver->nx;
    double h = solver->h;
    double k = solver->k;
    double c = h*h/k;
    double xmin = solver->xmin;
    double xmax = solver->xmax;

    // initialize and set parameters of the masa. 
    // This can give the value of the right sided function and the boundary condition
    masa_init("nick","heateq_1d_steady_const");
    masa_set_param("A_x", 10.0);
    masa_set_param("k_0", k);
    
    // set the matrix values, matrix index, non-zero numbers of each row, right sided vector values and boundary conditions
    if(solver->fd_method == 2) {
	// case of second order scheme
	//
	// set the parameters for first row
        solver->matrix_value[0] = (double *)malloc(2*sizeof(double));
        memcpy(solver->matrix_value[0], (double [2]){ 2, -1}, 2*sizeof(double));
        solver->matrix_index[0] = (int *)malloc(2*sizeof(int));
        memcpy(solver->matrix_index[0], (int [2]){ 0, 1}, 2*sizeof(int));
        solver->length[0] = 2;
        solver->rht_vec[0] = masa_eval_1d_source_t(xmin+h)*c + masa_eval_1d_exact_t(xmin);

	// set the parameters for inner rows
        for(int i = 1; i < N-2; i++) {
	    solver->matrix_value[i] = (double *)malloc(3*sizeof(double));       
            memcpy(solver->matrix_value[i], (double [3]){ -1, 2, -1}, 3*sizeof(double));
            solver->matrix_index[i] = (int *)malloc(3*sizeof(int));
            memcpy(solver->matrix_index[i], (int [3]){ i-1, i, i+1}, 3*sizeof(int));
            solver->length[i] = 3;
            solver->rht_vec[i] = masa_eval_1d_source_t(xmin+(i+1)*h)*c;
        }

	// set the parameters for last row
	solver->matrix_value[N-2] = (double *)malloc(2*sizeof(double));
        memcpy(solver->matrix_value[N-2], (double [2]){ -1, 2}, 2*sizeof(double));
        solver->matrix_index[N-2] = (int *)malloc(2*sizeof(int));
        memcpy(solver->matrix_index[N-2], (int [2]){N-3, N-2}, 2*sizeof(int));
        solver->length[N-2] = 2;
        solver->rht_vec[N-2] = masa_eval_1d_source_t(xmin+(N-1)*h)*c + masa_eval_1d_exact_t(xmax);
    
    } else {
	// case of fourth oeder scheme
	//
	// set the parameters for first row
	solver->matrix_value[0] = (double *)malloc(3*sizeof(double));
        memcpy(solver->matrix_value[0], (double [3]){ 30, -16, 1}, 3*sizeof(double));
        solver->matrix_index[0] = (int *)malloc(3*sizeof(int));
        memcpy(solver->matrix_index[0], (int [3]){ 0, 1, 2}, 3*sizeof(int));
        solver->length[0] = 3;
        solver->rht_vec[0] = masa_eval_1d_source_t(xmin+1*h)*12*c + 16*masa_eval_1d_exact_t(xmin) - masa_eval_1d_exact_t(xmin-h);	

	// set the parameters for second row
	solver->matrix_value[1] = (double *)malloc(4*sizeof(double));
        memcpy(solver->matrix_value[1], (double [4]){ -16, 30, -16, 1}, 4*sizeof(double));
        solver->matrix_index[1] = (int *)malloc(4*sizeof(int));
        memcpy(solver->matrix_index[1], (int [4]){ 0, 1, 2, 3}, 4*sizeof(int));
        solver->length[1] = 4;
        solver->rht_vec[1] = masa_eval_1d_source_t(xmin+2*h)*12*c - masa_eval_1d_exact_t(xmin); 

	// set the parameters for inner rows
        for(int i = 2; i < N-3; i++) {
            solver->matrix_value[i] = (double *)malloc(5*sizeof(double));
            memcpy(solver->matrix_value[i], (double [5]){ 1, -16, 30, -16, 1}, 5*sizeof(double));
            solver->matrix_index[i] = (int *)malloc(5*sizeof(int));
            memcpy(solver->matrix_index[i], (int [5]){ i-2, i-1, i, i+1, i+2}, 5*sizeof(int));
            solver->length[i] = 5;
            solver->rht_vec[i] = masa_eval_1d_source_t(xmin+(i+1)*h)*12*c;
        }

	// set the parameters for last second row
        solver->matrix_value[N-3] = (double *)malloc(4*sizeof(double));
        memcpy(solver->matrix_value[N-3], (double [4]){ 1, -16, 30, -16}, 4*sizeof(double));
        solver->matrix_index[N-3] = (int *)malloc(4*sizeof(int));
        memcpy(solver->matrix_index[N-3], (int [4]){ N-5, N-4, N-3, N-2}, 4*sizeof(int));
        solver->length[N-3] = 4;
        solver->rht_vec[N-3] = masa_eval_1d_source_t(xmin+(N-2)*h)*12*c - masa_eval_1d_exact_t(xmax); 

	// set the parameters for last row
	solver->matrix_value[N-2] = (double *)malloc(3*sizeof(double));
        memcpy(solver->matrix_value[N-2], (double [3]){ 1, -16, 30}, 3*sizeof(double));
        solver->matrix_index[N-2] = (int *)malloc(3*sizeof(int));
        memcpy(solver->matrix_index[N-2], (int [3]){ N-4, N-3, N-2}, 3*sizeof(int));
        solver->length[N-2] = 3;
        solver->rht_vec[N-2] = masa_eval_1d_source_t(xmin+(N-1)*h)*12*c + 16*masa_eval_1d_exact_t(xmax) - masa_eval_1d_exact_t(xmax+h);
    
    }

    if(solver->output_mode == 1) {
        show_linear_system(solver);
        printf("[debug]: build linear system - function end\n\n\n");
    }

    grvy_timer_end(__func__);

}
