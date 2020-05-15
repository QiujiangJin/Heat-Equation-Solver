/******************************************************************
 * This function builds the linear system in two dimensional case
 * ***************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <grvy.h>
#include <masa.h>
#include "laplace.h"
#include "function.h"

void build_linear_system_2d(Laplace* solver) {

    grvy_timer_begin(__func__);

    printf("** Building two dimensional linear system...\n");
    printf("--> Enforcing analytic Dirichlet BCs using MASA (2D)\n\n\n");

    if(solver->output_mode == 1) {
        printf("[debug]: build linear system - function begin\n\n\n");
    }

    int N = solver->nx;
    double h = solver->h;
    double k = solver->k;
    double c = h*h/k;
    double xmin = solver->xmin;
    double xmax = solver->xmax;
    double ymin = solver->ymin;
    double ymax = solver->ymax;

    // initialize and set parameters of the masa. 
    // This can give the value of the right sided function and the boundary condition
    masa_init("nick","heateq_2d_steady_const");
    masa_set_param("A_x", 10.0);
    masa_set_param("B_y", 10.0);
    masa_set_param("k_0", k);

    // set the matrix values, matrix index, non-zero numbers of each row, right sided vector values and boundary conditions
    if(solver->fd_method == 2) {
	// case of second order scheme
	// 
	// set the parameters of first block row
	// set the parameters for fisrt row of first block row
	solver->matrix_value[0] = (double *)malloc(3*sizeof(double));
        memcpy(solver->matrix_value[0], (double [3]){ 4, -1, -1}, 3*sizeof(double));
        solver->matrix_index[0] = (int *)malloc(3*sizeof(int));
        memcpy(solver->matrix_index[0], (int [3]){ 0, 1, N-1}, 3*sizeof(int));
        solver->length[0] = 3;
        solver->rht_vec[0] = masa_eval_2d_source_t(xmin+h, ymin+h)*c + masa_eval_2d_exact_t(xmin, ymin+h) + masa_eval_2d_exact_t(xmin+h, ymin);

	// set the parameters for inner rows of first blcok row
        for(int j = 1; j < N-2; j++) {
            solver->matrix_value[j] = (double *)malloc(4*sizeof(double));
            memcpy(solver->matrix_value[j], (double [4]){ -1, 4, -1, -1}, 4*sizeof(double));
            solver->matrix_index[j] = (int *)malloc(4*sizeof(int));
            memcpy(solver->matrix_index[j], (int [4]){ j-1, j, j+1, N-1+j}, 4*sizeof(int));
            solver->length[j] = 4;
            solver->rht_vec[j] = masa_eval_2d_source_t(xmin+(j+1)*h, ymin+h)*c + masa_eval_2d_exact_t(xmin+(j+1)*h, ymin);
        }

	// set the parameters for last row of first block row
        solver->matrix_value[N-2] = (double *)malloc(3*sizeof(double));
        memcpy(solver->matrix_value[N-2], (double [3]){ -1, 4, -1}, 3*sizeof(double));
        solver->matrix_index[N-2] = (int *)malloc(3*sizeof(int));
        memcpy(solver->matrix_index[N-2], (int [3]){ N-3, N-2, 2*N-3}, 3*sizeof(int));
        solver->length[N-2] = 3;
        solver->rht_vec[N-2] = masa_eval_2d_source_t(xmin+(N-1)*h, ymin+h)*c + masa_eval_2d_exact_t(xmin+(N-1)*h, ymin) + masa_eval_2d_exact_t(xmax, ymin+h);
    
	//
	// set the parameters of inner block rows
        for(int i = 1; i < N-2; i++) {
            // set the parameters for first row of inner block rows
            solver->matrix_value[i*(N-1)] = (double *)malloc(4*sizeof(double));
            memcpy(solver->matrix_value[i*(N-1)], (double [4]){ -1, 4, -1, -1}, 4*sizeof(double));
            solver->matrix_index[i*(N-1)] = (int *)malloc(4*sizeof(int));
            memcpy(solver->matrix_index[i*(N-1)], (int [4]){ (i-1)*(N-1), i*(N-1), i*(N-1)+1, (i+1)*(N-1)}, 4*sizeof(int));
            solver->length[i*(N-1)] = 4;
            solver->rht_vec[i*(N-1)] = masa_eval_2d_source_t(xmin+h, ymin+(i+1)*h)*c + masa_eval_2d_exact_t(xmin, ymin+(i+1)*h);

	    // set the parameters for inner rows of inner block rows
            for(int j = 1; j < N-2; j++) {
                solver->matrix_value[i*(N-1)+j] = (double *)malloc(5*sizeof(double));
                memcpy(solver->matrix_value[i*(N-1)+j], (double [5]){ -1, -1, 4, -1, -1}, 5*sizeof(double));
                solver->matrix_index[i*(N-1)+j] = (int *)malloc(5*sizeof(int));
                memcpy(solver->matrix_index[i*(N-1)+j], (int [5]){ (i-1)*(N-1)+j, i*(N-1)+j-1, i*(N-1)+j, i*(N-1)+j+1, (i+1)*(N-1)+j}, 5*sizeof(int));
                solver->length[i*(N-1)+j] = 5;
                solver->rht_vec[i*(N-1)+j] = masa_eval_2d_source_t(xmin+(j+1)*h, ymin+(i+1)*h)*c;
            }

	    // set the parameters for last row of inner block rows
            solver->matrix_value[i*(N-1)+N-2] = (double *)malloc(4*sizeof(double));
            memcpy(solver->matrix_value[i*(N-1)+N-2], (double [4]){ -1, -1, 4, -1}, 4*sizeof(double));
            solver->matrix_index[i*(N-1)+N-2] = (int *)malloc(4*sizeof(int));
            memcpy(solver->matrix_index[i*(N-1)+N-2], (int [4]){ (i-1)*(N-1)+N-2, i*(N-1)+N-3, i*(N-1)+N-2, (i+1)*(N-1)+N-2}, 4*sizeof(int));
            solver->length[i*(N-1)+N-2] = 4;
            solver->rht_vec[i*(N-1)+N-2] = masa_eval_2d_source_t(xmin+(N-1)*h, ymin+(i+1)*h)*c + masa_eval_2d_exact_t(xmax, ymin+(i+1)*h);

        }

	//
	// set the parameters of last block row
	// set the parameters for first row of last block row
	solver->matrix_value[(N-2)*(N-1)] = (double *)malloc(3*sizeof(double));
        memcpy(solver->matrix_value[(N-2)*(N-1)], (double [3]){ -1, 4, -1}, 3*sizeof(double));
        solver->matrix_index[(N-2)*(N-1)] = (int *)malloc(3*sizeof(int));
        memcpy(solver->matrix_index[(N-2)*(N-1)], (int [3]){ (N-1)*(N-1)+2-2*N, (N-1)*(N-1)+1-N, (N-1)*(N-1)+2-N}, 3*sizeof(int));
        solver->length[(N-2)*(N-1)] = 3;
        solver->rht_vec[(N-2)*(N-1)] = masa_eval_2d_source_t(xmin+h, ymin+(N-1)*h)*c + masa_eval_2d_exact_t(xmin, ymin+(N-1)*h) + masa_eval_2d_exact_t(xmin+h, ymax);

	// set the parameters for inner rows of last block row
        for(int j = 1; j < N-2; j++) {
            solver->matrix_value[(N-2)*(N-1)+j] = (double *)malloc(4*sizeof(double));
            memcpy(solver->matrix_value[(N-2)*(N-1)+j], (double [4]){ -1, -1, 4, -1}, 4*sizeof(double));
            solver->matrix_index[(N-2)*(N-1)+j] = (int *)malloc(4*sizeof(int));
            memcpy(solver->matrix_index[(N-2)*(N-1)+j], (int [4]){ (N-1)*(N-1)+2-2*N+j, (N-1)*(N-1)-N+j, (N-1)*(N-1)-N+j+1, (N-1)*(N-1)-N+j+2}, 4*sizeof(int));
            solver->length[(N-2)*(N-1)+j] = 4;
            solver->rht_vec[(N-2)*(N-1)+j] = masa_eval_2d_source_t(xmin+(j+1)*h, ymin+(N-1)*h)*c + masa_eval_2d_exact_t(xmin+(j+1)*h, ymax);
        }

	// set the parameters for last row of last block row
        solver->matrix_value[N*(N-2)] = (double *)malloc(3*sizeof(double));
        memcpy(solver->matrix_value[N*(N-2)], (double [3]){ -1, -1, 4}, 3*sizeof(double));
        solver->matrix_index[N*(N-2)] = (int *)malloc(3*sizeof(int));
        memcpy(solver->matrix_index[N*(N-2)], (int [3]){ (N-1)*(N-1)-N, (N-1)*(N-1)-2, (N-1)*(N-1)-1}, 3*sizeof(int));
        solver->length[N*(N-2)] = 3;
        solver->rht_vec[N*(N-2)] = masa_eval_2d_source_t(xmin+(N-1)*h, ymin+(N-1)*h)*c + masa_eval_2d_exact_t(xmax, ymin+(N-1)*h) + masa_eval_2d_exact_t(xmin+(N-1)*h, ymax);

    } else {
	// case of fourth order scheme
	// 
	// set the parameters of first block row
	// set the parameters for first row of first block row
	solver->matrix_value[0] = (double *)malloc(5*sizeof(double));
        memcpy(solver->matrix_value[0], (double [5]){ 60, -16, 1, -16, 1}, 5*sizeof(double));
        solver->matrix_index[0] = (int *)malloc(5*sizeof(int));
        memcpy(solver->matrix_index[0], (int [5]){ 0, 1, 2, N-1, 2*(N-1)}, 5*sizeof(int));
        solver->length[0] = 5;
        solver->rht_vec[0] = masa_eval_2d_source_t(xmin+h, ymin+h)*12*c + 16*masa_eval_2d_exact_t(xmin, ymin+h) + 16*masa_eval_2d_exact_t(xmin+h, ymin) - masa_eval_2d_exact_t(xmin-h, ymin+h) - masa_eval_2d_exact_t(xmin+h, ymin-h);

	// set the parameters for second row of first block row
	solver->matrix_value[1] = (double *)malloc(6*sizeof(double));
        memcpy(solver->matrix_value[1], (double [6]){ -16, 60, -16, 1, -16, 1}, 6*sizeof(double));
        solver->matrix_index[1] = (int *)malloc(6*sizeof(int));
        memcpy(solver->matrix_index[1], (int [6]){ 0, 1, 2, 3, N, 2*N-1}, 6*sizeof(int));
        solver->length[1] = 6;
        solver->rht_vec[1] = masa_eval_2d_source_t(xmin+2*h, ymin+h)*12*c + 16*masa_eval_2d_exact_t(xmin+2*h, ymin) - masa_eval_2d_exact_t(xmin, ymin+h) - masa_eval_2d_exact_t(xmin+2*h, ymin-h);

	// set the parameters for inner rows of first block row
        for(int j = 2; j < N-3; j++) {
            solver->matrix_value[j] = (double *)malloc(7*sizeof(double));
            memcpy(solver->matrix_value[j], (double [7]){ 1, -16, 60, -16, 1, -16, 1}, 7*sizeof(double));
            solver->matrix_index[j] = (int *)malloc(7*sizeof(int));
            memcpy(solver->matrix_index[j], (int [7]){ j-2, j-1, j, j+1, j+2, N-1+j, 2*N-2+j}, 7*sizeof(int));
            solver->length[j] = 7;
            solver->rht_vec[j] = masa_eval_2d_source_t(xmin+(j+1)*h, ymin+h)*12*c + 16*masa_eval_2d_exact_t(xmin+(j+1)*h, ymin) - masa_eval_2d_exact_t(xmin+(j+1)*h, ymin-h);
        }

	// set the parameters for second last row of first block row
	solver->matrix_value[N-3] = (double *)malloc(6*sizeof(double));
        memcpy(solver->matrix_value[N-3], (double [6]){ 1, -16, 60, -16, -16, 1}, 6*sizeof(double));
        solver->matrix_index[N-3] = (int *)malloc(6*sizeof(int));
        memcpy(solver->matrix_index[N-3], (int [6]){ N-5, N-4, N-3, N-2, 2*N-4, 3*N-5}, 6*sizeof(int));
        solver->length[N-3] = 6;
        solver->rht_vec[N-3] = masa_eval_2d_source_t(xmin+(N-2)*h, ymin+h)*12*c + 16*masa_eval_2d_exact_t(xmin+(N-2)*h, ymin) - masa_eval_2d_exact_t(xmax, ymin+h) - masa_eval_2d_exact_t(xmin+(N-2)*h, ymin-h);

	// set the parameters for last row of first block row
        solver->matrix_value[N-2] = (double *)malloc(5*sizeof(double));
        memcpy(solver->matrix_value[N-2], (double [5]){ 1, -16, 60, -16, 1}, 5*sizeof(double));
        solver->matrix_index[N-2] = (int *)malloc(5*sizeof(int));
        memcpy(solver->matrix_index[N-2], (int [5]){ N-4, N-3, N-2, 2*N-3, 3*N-4}, 5*sizeof(int));
        solver->length[N-2] = 5;
        solver->rht_vec[N-2] = masa_eval_2d_source_t(xmin+(N-1)*h, ymin+h)*12*c + 16*masa_eval_2d_exact_t(xmin+(N-1)*h, ymin) + 16*masa_eval_2d_exact_t(xmax, ymin+h) - masa_eval_2d_exact_t(xmin+(N-1)*h, ymin-h) - masa_eval_2d_exact_t(xmax+h, ymin+h);	

	// 
	// set the parameters of second block row
	// set the parameters for first row of second block row
	solver->matrix_value[N-1] = (double *)malloc(6*sizeof(double));
        memcpy(solver->matrix_value[N-1], (double [6]){ -16, 60, -16, 1, -16, 1}, 6*sizeof(double));
        solver->matrix_index[N-1] = (int *)malloc(6*sizeof(int));
        memcpy(solver->matrix_index[N-1], (int [6]){ 0, N-1, N, N+1, 2*(N-1), 3*(N-1)}, 6*sizeof(int));
        solver->length[N-1] = 6;
        solver->rht_vec[N-1] = masa_eval_2d_source_t(xmin+h, ymin+2*h)*12*c + 16*masa_eval_2d_exact_t(xmin, ymin+2*h) - masa_eval_2d_exact_t(xmin-h, ymin+2*h) - masa_eval_2d_exact_t(xmin+h, ymin);

	// set the parameters for second row of second block row
        solver->matrix_value[N] = (double *)malloc(7*sizeof(double));
        memcpy(solver->matrix_value[N], (double [7]){ -16, -16, 60, -16, 1, -16, 1}, 7*sizeof(double));
        solver->matrix_index[N] = (int *)malloc(7*sizeof(int));
        memcpy(solver->matrix_index[N], (int [7]){ 1, N-1, N, N+1, N+2, 2*N-1, 3*N-2}, 7*sizeof(int));
        solver->length[N] = 7;
        solver->rht_vec[N] = masa_eval_2d_source_t(xmin+2*h, ymin+2*h)*12*c - masa_eval_2d_exact_t(xmin+2*h, ymin) - masa_eval_2d_exact_t(xmin, ymin+2*h);

	// set the parameters for inner rows of second block row
        for(int j = 2; j < N-3; j++) {
            solver->matrix_value[N-1+j] = (double *)malloc(8*sizeof(double));
            memcpy(solver->matrix_value[N-1+j], (double [8]){ -16, 1, -16, 60, -16, 1, -16, 1}, 8*sizeof(double));
            solver->matrix_index[N-1+j] = (int *)malloc(8*sizeof(int));
            memcpy(solver->matrix_index[N-1+j], (int [8]){ j, N+j-3, N+j-2, N+j-1, N+j, N+j+1, 2*N+j-2, 3*N+j-3}, 8*sizeof(int));
            solver->length[N-1+j] = 8;
            solver->rht_vec[N-1+j] = masa_eval_2d_source_t(xmin+(j+1)*h, ymin+2*h)*12*c - masa_eval_2d_exact_t(xmin+(j+1)*h, ymin);
        }

	// set the parameters for second last row of second block row
	solver->matrix_value[2*N-4] = (double *)malloc(7*sizeof(double));
        memcpy(solver->matrix_value[2*N-4], (double [7]){ -16, 1, -16, 60, -16, -16, 1}, 7*sizeof(double));
        solver->matrix_index[2*N-4] = (int *)malloc(7*sizeof(int));
        memcpy(solver->matrix_index[2*N-4], (int [7]){ N-3, 2*N-6, 2*N-5, 2*N-4, 2*N-3, 3*N-5, 4*N-6}, 7*sizeof(int));
        solver->length[2*N-4] = 7;
        solver->rht_vec[2*N-4] = masa_eval_2d_source_t(xmin+(N-2)*h, ymin+2*h)*12*c - masa_eval_2d_exact_t(xmax, ymin+2*h) - masa_eval_2d_exact_t(xmin+(N-2)*h, ymin);

	// set the parameters for last row of second block row
        solver->matrix_value[2*N-3] = (double *)malloc(6*sizeof(double));
        memcpy(solver->matrix_value[2*N-3], (double [6]){ -16, 1, -16, 60, -16, 1}, 6*sizeof(double));
        solver->matrix_index[2*N-3] = (int *)malloc(6*sizeof(int));
        memcpy(solver->matrix_index[2*N-3], (int [6]){ N-2, 2*N-5, 2*N-4, 2*N-3, 3*N-4, 4*N-5}, 6*sizeof(int));
        solver->length[2*N-3] = 6;
        solver->rht_vec[2*N-3] = masa_eval_2d_source_t(xmin+(N-1)*h, ymin+2*h)*12*c + 16*masa_eval_2d_exact_t(xmax, ymin+2*h) - masa_eval_2d_exact_t(xmax+h, ymin+2*h) - masa_eval_2d_exact_t(xmin+(N-1)*h, ymin);

	// 
	// set the parameters of inner block rows
	for(int i = 2; i < N-3; i++) {
	    // set the parameters for first row of inner block rows
	    solver->matrix_value[i*(N-1)] = (double *)malloc(7*sizeof(double));
            memcpy(solver->matrix_value[i*(N-1)], (double [7]){ 1, -16, 60, -16, 1, -16, 1}, 7*sizeof(double));
            solver->matrix_index[i*(N-1)] = (int *)malloc(7*sizeof(int));
            memcpy(solver->matrix_index[i*(N-1)], (int [7]){ (i-2)*(N-1), (i-1)*(N-1),i*(N-1), i*(N-1)+1, i*(N-1)+2, (i+1)*(N-1), (i+2)*(N-1)}, 7*sizeof(int));
            solver->length[i*(N-1)] = 7;
            solver->rht_vec[i*(N-1)] = masa_eval_2d_source_t(xmin+h, ymin+(i+1)*h)*12*c + 16*masa_eval_2d_exact_t(xmin, ymin+(i+1)*h) - masa_eval_2d_exact_t(xmin-h, ymin+(i+1)*h);
	    // set the parameters for second row of inner block rows
            solver->matrix_value[i*(N-1)+1] = (double *)malloc(8*sizeof(double));
            memcpy(solver->matrix_value[i*(N-1)+1], (double [8]){ 1, -16, -16, 60, -16, 1, -16, 1}, 8*sizeof(double));
            solver->matrix_index[i*(N-1)+1] = (int *)malloc(8*sizeof(int));
            memcpy(solver->matrix_index[i*(N-1)+1], (int [8]){ (i-2)*(N-1)+1, (i-1)*(N-1)+1, i*(N-1), i*(N-1)+1, i*(N-1)+2, i*(N-1)+3, (i+1)*(N-1)+1, (i+2)*(N-1)+1}, 8*sizeof(int));
            solver->length[i*(N-1)+1] = 8;
            solver->rht_vec[i*(N-1)+1] = masa_eval_2d_source_t(xmin+2*h, ymin+(i+1)*h)*12*c - masa_eval_2d_exact_t(xmin, ymin+(i+1)*h);

	    // set the parameters for inner rows of inner block rows
            for(int j = 2; j < N-3; j++) {
                solver->matrix_value[i*(N-1)+j] = (double *)malloc(9*sizeof(double));
                memcpy(solver->matrix_value[i*(N-1)+j], (double [9]){ 1, -16, 1, -16, 60, -16, 1, -16, 1}, 9*sizeof(double));
                solver->matrix_index[i*(N-1)+j] = (int *)malloc(9*sizeof(int));
                memcpy(solver->matrix_index[i*(N-1)+j], (int [9]){ (i-2)*(N-1)+j, (i-1)*(N-1)+j, i*(N-1)+j-2, i*(N-1)+j-1, i*(N-1)+j, i*(N-1)+j+1, i*(N-1)+j+2, (i+1)*(N-1)+j, (i+2)*(N-1)+j}, 9*sizeof(int));
                solver->length[i*(N-1)+j] = 9;
                solver->rht_vec[i*(N-1)+j] = masa_eval_2d_source_t(xmin+(j+1)*h, ymin+(i+1)*h)*12*c;
            }

	    // set the parameters for second last row of inner block rows
	    solver->matrix_value[i*(N-1)+N-3] = (double *)malloc(8*sizeof(double));
            memcpy(solver->matrix_value[i*(N-1)+N-3], (double [8]){ 1, -16, 1, -16, 60, -16, -16, 1}, 8*sizeof(double));
            solver->matrix_index[i*(N-1)+N-3] = (int *)malloc(8*sizeof(int));
            memcpy(solver->matrix_index[i*(N-1)+N-3], (int [8]){ (i-2)*(N-1)+N-3, (i-1)*(N-1)+N-3, i*(N-1)+N-5, i*(N-1)+N-4, i*(N-1)+N-3, i*(N-1)+N-2, (i+1)*(N-1)+N-3, (i+2)*(N-1)+N-3}, 8*sizeof(int));
            solver->length[i*(N-1)+N-3] = 8;
            solver->rht_vec[i*(N-1)+N-3] = masa_eval_2d_source_t(xmin+(N-2)*h, ymin+(i+1)*h)*12*c - masa_eval_2d_exact_t(xmax, ymin+(i+1)*h);

	    // set the parameters for last row of inner block rows
            solver->matrix_value[i*(N-1)+N-2] = (double *)malloc(7*sizeof(double));
            memcpy(solver->matrix_value[i*(N-1)+N-2], (double [7]){ 1, -16, 1, -16, 60, -16, 1}, 7*sizeof(double));
            solver->matrix_index[i*(N-1)+N-2] = (int *)malloc(7*sizeof(int));
            memcpy(solver->matrix_index[i*(N-1)+N-2], (int [7]){ (i-2)*(N-1)+N-2, (i-1)*(N-1)+N-2, i*(N-1)+N-4, i*(N-1)+N-3, i*(N-1)+N-2, (i+1)*(N-1)+N-2, (i+2)*(N-1)+N-2}, 7*sizeof(int));
            solver->length[i*(N-1)+N-2] = 7;
            solver->rht_vec[i*(N-1)+N-2] = masa_eval_2d_source_t(xmin+(N-1)*h, ymin+(i+1)*h)*12*c + 16*masa_eval_2d_exact_t(xmax, ymin+(i+1)*h) - masa_eval_2d_exact_t(xmax+h, ymin+(i+1)*h);

	}

	// 
	// set the parameters of second last block row
	// set the parameters for first row of second last block row
	solver->matrix_value[(N-3)*(N-1)] = (double *)malloc(6*sizeof(double));
        memcpy(solver->matrix_value[(N-3)*(N-1)], (double [6]){ 1, -16, 60, -16, 1, -16}, 6*sizeof(double));
        solver->matrix_index[(N-3)*(N-1)] = (int *)malloc(6*sizeof(int));
        memcpy(solver->matrix_index[(N-3)*(N-1)], (int [6]){ (N-5)*(N-1), (N-4)*(N-1), (N-3)*(N-1), (N-3)*(N-1)+1, (N-3)*(N-1)+2, (N-2)*(N-1)}, 6*sizeof(int));
        solver->length[(N-3)*(N-1)] = 6;
        solver->rht_vec[(N-3)*(N-1)] = masa_eval_2d_source_t(xmin+h, ymin+(N-2)*h)*12*c + 16*masa_eval_2d_exact_t(xmin, ymin+(N-2)*h) - masa_eval_2d_exact_t(xmin-h, ymin+(N-2)*h) - masa_eval_2d_exact_t(xmin+h, ymax);

	// set the parameters for second row of second last block row
        solver->matrix_value[(N-3)*(N-1)+1] = (double *)malloc(7*sizeof(double));
        memcpy(solver->matrix_value[(N-3)*(N-1)+1], (double [7]){ 1, -16, -16, 60, -16, 1, -16}, 7*sizeof(double));
        solver->matrix_index[(N-3)*(N-1)+1] = (int *)malloc(7*sizeof(int));
        memcpy(solver->matrix_index[(N-3)*(N-1)+1], (int [7]){ (N-5)*(N-1)+1, (N-4)*(N-1)+1, (N-3)*(N-1), (N-3)*(N-1)+1, (N-3)*(N-1)+2, (N-3)*(N-1)+3, (N-2)*(N-1)+1}, 7*sizeof(int));
        solver->length[(N-3)*(N-1)+1] = 7;
        solver->rht_vec[(N-3)*(N-1)+1] = masa_eval_2d_source_t(xmin+2*h, ymin+(N-2)*h)*12*c - masa_eval_2d_exact_t(xmin, ymin+(N-2)*h) - masa_eval_2d_exact_t(xmin+2*h, ymax);

	// set the parameters for inner rows of second last block row
        for(int j = 2; j < N-3; j++) {
            solver->matrix_value[(N-3)*(N-1)+j] = (double *)malloc(8*sizeof(double));
            memcpy(solver->matrix_value[(N-3)*(N-1)+j], (double [8]){ 1, -16, 1, -16, 60, -16, 1, -16}, 8*sizeof(double));
            solver->matrix_index[(N-3)*(N-1)+j] = (int *)malloc(8*sizeof(int));
            memcpy(solver->matrix_index[(N-3)*(N-1)+j], (int [8]){ (N-5)*(N-1)+j, (N-4)*(N-1)+j, (N-3)*(N-1)+j-2, (N-3)*(N-1)+j-1, (N-3)*(N-1)+j, (N-3)*(N-1)+j+1, (N-3)*(N-1)+j+2, (N-2)*(N-1)+j}, 8*sizeof(int));
            solver->length[(N-3)*(N-1)+j] = 8;
            solver->rht_vec[(N-3)*(N-1)+j] = masa_eval_2d_source_t(xmin+(j+1)*h, ymin+(N-2)*h)*12*c - masa_eval_2d_exact_t(xmin+(j+1)*h, ymax);
        }

	// set the parameters for second last row of second last block row
	solver->matrix_value[(N-3)*N] = (double *)malloc(7*sizeof(double));
        memcpy(solver->matrix_value[(N-3)*N], (double [7]){ 1, -16, 1, -16, 60, -16, -16}, 7*sizeof(double));
        solver->matrix_index[(N-3)*N] = (int *)malloc(7*sizeof(int));
        memcpy(solver->matrix_index[(N-3)*N], (int [7]){ (N-5)*(N-1)+N-3, (N-4)*(N-1)+N-3, (N-3)*(N-1)+N-5, (N-3)*(N-1)+N-4, (N-3)*(N-1)+N-3, (N-3)*(N-1)+N-2, (N-2)*(N-1)+N-3}, 7*sizeof(int));
        solver->length[(N-3)*N] = 7;
        solver->rht_vec[(N-3)*N] = masa_eval_2d_source_t(xmin+(N-2)*h, ymin+(N-2)*h)*12*c - masa_eval_2d_exact_t(xmax, ymin+(N-2)*h) - masa_eval_2d_exact_t(xmin+(N-2)*h, ymax);

	// set the parameters for last row of second last block row
        solver->matrix_value[(N-3)*N+1] = (double *)malloc(6*sizeof(double));
        memcpy(solver->matrix_value[(N-3)*N+1], (double [6]){ 1, -16, 1, -16, 60, -16}, 6*sizeof(double));
        solver->matrix_index[(N-3)*N+1] = (int *)malloc(6*sizeof(int));
        memcpy(solver->matrix_index[(N-3)*N+1], (int [6]){ (N-5)*(N-1)+N-2, (N-4)*(N-1)+N-2, (N-3)*(N-1)+N-4, (N-3)*(N-1)+N-3, (N-3)*(N-1)+N-2, (N-2)*(N-1)+N-2}, 6*sizeof(int));
        solver->length[(N-3)*N+1] = 6;
        solver->rht_vec[(N-3)*N+1] = masa_eval_2d_source_t(xmin+(N-1)*h, ymin+(N-2)*h)*12*c + 16*masa_eval_2d_exact_t(xmax, ymin+(N-2)*h) - masa_eval_2d_exact_t(xmax+h, ymin+(N-2)*h) - masa_eval_2d_exact_t(xmin+(N-1)*h, ymax);

	// 
	// set the parameters of last block row
	// set the parameters for first row of last block row
	solver->matrix_value[(N-2)*(N-1)] = (double *)malloc(5*sizeof(double));
        memcpy(solver->matrix_value[(N-2)*(N-1)], (double [5]){ 1, -16, 60, -16, 1}, 5*sizeof(double));
        solver->matrix_index[(N-2)*(N-1)] = (int *)malloc(5*sizeof(int));
        memcpy(solver->matrix_index[(N-2)*(N-1)], (int [5]){ (N-4)*(N-1), (N-3)*(N-1), (N-2)*(N-1), (N-2)*(N-1)+1, (N-2)*(N-1)+2}, 5*sizeof(int));
        solver->length[(N-2)*(N-1)] = 5;
        solver->rht_vec[(N-2)*(N-1)] = masa_eval_2d_source_t(xmin+h, ymin+(N-1)*h)*12*c + 16*masa_eval_2d_exact_t(xmin, ymin+(N-1)*h) + 16*masa_eval_2d_exact_t(xmin+h, ymax) - masa_eval_2d_exact_t(xmin-h, ymin+(N-1)*h) - masa_eval_2d_exact_t(xmin+h, ymax+h);

	// set the parameters for second row of last block row
        solver->matrix_value[(N-2)*(N-1)+1] = (double *)malloc(6*sizeof(double));
        memcpy(solver->matrix_value[(N-2)*(N-1)+1], (double [6]){ 1, -16, -16, 60, -16, 1}, 6*sizeof(double));
        solver->matrix_index[(N-2)*(N-1)+1] = (int *)malloc(6*sizeof(int));
        memcpy(solver->matrix_index[(N-2)*(N-1)+1], (int [6]){ (N-4)*(N-1)+1, (N-3)*(N-1)+1, (N-2)*(N-1), (N-2)*(N-1)+1, (N-2)*(N-1)+2, (N-2)*(N-1)+3}, 6*sizeof(int));
        solver->length[(N-2)*(N-1)+1] = 6;
        solver->rht_vec[(N-2)*(N-1)+1] = masa_eval_2d_source_t(xmin+2*h, ymin+(N-1)*h)*12*c + 16*masa_eval_2d_exact_t(xmin+2*h, ymax) - masa_eval_2d_exact_t(xmin, ymin+(N-1)*h) - masa_eval_2d_exact_t(xmin+2*h, ymax+h);

	// set the parameters for inner rows of last block row
        for(int j = 2; j < N-3; j++) {
            solver->matrix_value[(N-2)*(N-1)+j] = (double *)malloc(7*sizeof(double));
            memcpy(solver->matrix_value[(N-2)*(N-1)+j], (double [7]){ 1, -16, 1, -16, 60, -16, 1}, 7*sizeof(double));
            solver->matrix_index[(N-2)*(N-1)+j] = (int *)malloc(7*sizeof(int));
            memcpy(solver->matrix_index[(N-2)*(N-1)+j], (int [7]){ (N-4)*(N-1)+j, (N-3)*(N-1)+j, (N-2)*(N-1)+j-2, (N-2)*(N-1)+j-1, (N-2)*(N-1)+j, (N-2)*(N-1)+j+1, (N-2)*(N-1)+j+2}, 7*sizeof(int));
            solver->length[(N-2)*(N-1)+j] = 7;
            solver->rht_vec[(N-2)*(N-1)+j] = masa_eval_2d_source_t(xmin+(j+1)*h, ymin+(N-1)*h)*12*c + 16*masa_eval_2d_exact_t(xmin+(j+1)*h, ymax) - masa_eval_2d_exact_t(xmin+(j+1)*h, ymax+h);
        }

	// set the parameters for second last row of last block row
	solver->matrix_value[(N-2)*(N-1)+N-3] = (double *)malloc(6*sizeof(double));
        memcpy(solver->matrix_value[(N-2)*(N-1)+N-3], (double [6]){ 1, -16, 1, -16, 60, -16}, 6*sizeof(double));
        solver->matrix_index[(N-2)*(N-1)+N-3] = (int *)malloc(6*sizeof(int));
        memcpy(solver->matrix_index[(N-2)*(N-1)+N-3], (int [6]){ (N-4)*(N-1)+N-3, (N-3)*(N-1)+N-3, (N-2)*(N-1)+N-5, (N-2)*(N-1)+N-4, (N-2)*(N-1)+N-3, (N-2)*(N-1)+N-2}, 6*sizeof(int));
        solver->length[(N-2)*(N-1)+N-3] = 6;
        solver->rht_vec[(N-2)*(N-1)+N-3] = masa_eval_2d_source_t(xmin+(N-2)*h, ymin+(N-1)*h)*12*c + 16*masa_eval_2d_exact_t(xmin+(N-2)*h, ymax) - masa_eval_2d_exact_t(xmax, ymin+(N-1)*h) - masa_eval_2d_exact_t(xmin+(N-2)*h, ymax+h);

	// set the parameters for last row of last block row
        solver->matrix_value[(N-2)*N] = (double *)malloc(5*sizeof(double));
        memcpy(solver->matrix_value[(N-2)*N], (double [5]){ 1, -16, 1, -16, 60}, 5*sizeof(double));
        solver->matrix_index[(N-2)*N] = (int *)malloc(5*sizeof(int));
        memcpy(solver->matrix_index[(N-2)*N], (int [5]){ (N-4)*(N-1)+N-2, (N-3)*(N-1)+N-2, (N-2)*(N-1)+N-4, (N-2)*(N-1)+N-3, (N-2)*(N-1)+N-2}, 5*sizeof(int));
        solver->length[(N-2)*N] = 5;
        solver->rht_vec[(N-2)*N] = masa_eval_2d_source_t(xmin+(N-1)*h, ymin+(N-1)*h)*12*c + 16*masa_eval_2d_exact_t(xmin+(N-1)*h, ymax) + 16*masa_eval_2d_exact_t(xmax, ymin+(N-1)*h) - masa_eval_2d_exact_t(xmin+(N-1)*h, ymax+h) - masa_eval_2d_exact_t(xmax+h, ymin+(N-1)*h);

    }

    if(solver->output_mode == 1) {
	show_linear_system(solver);
        printf("[debug]: build linear system - function end\n\n\n");
    }

    grvy_timer_end(__func__);

}
