/********************************************************
 * This function outputs the result to the output_file
 * *****************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <grvy.h>
#include <masa.h>
#include "laplace.h"
#include "function.h"

void output(Laplace* solver) {

    grvy_timer_begin(__func__);

    if(solver->output_mode == 1) {
        printf("[debug]: output - function begin\n\n\n");
    }

    printf("** Writing solution to %s\n\n\n", solver->output_file);

    int N = solver->nx;
    double xmin = solver->xmin;
    double xmax = solver->xmax;
    double ymin = solver->ymin;
    double ymax = solver->ymax;
    double k = solver->k;
    double h = solver->h;

    // output the numerical solution to the file
    FILE *f = fopen(solver->output_file, "w");
    if(f == NULL)  {
        printf("Error opening file!\n");
        exit(1);
    }

    if(solver->dimension == 1) {
	masa_init("nick","heateq_1d_steady_const");
        masa_set_param("A_x", 10.0);
        masa_set_param("k_0", k);
	fprintf(f, "%f\n", masa_eval_1d_exact_t(xmin));
	for(int i = 0; i < N-1; i++) {
	    fprintf(f, "%f\n", solver->sol_vec[i]);
	}
	fprintf(f, "%f\n", masa_eval_1d_exact_t(xmax));
    } else {
	masa_init("nick","heateq_2d_steady_const");
        masa_set_param("A_x", 10.0);
        masa_set_param("B_y", 10.0);
        masa_set_param("k_0", k);
	for(int j = -1; j < N; j++) {
	    fprintf(f, "%f\n", masa_eval_2d_exact_t(xmin+(j+1)*h, ymin));
	}
	for(int i = 0; i < N-1; i++) {
	    fprintf(f, "%f\n", masa_eval_2d_exact_t(xmin, ymin+(i+1)*h));
            for(int j = 0; j < N-1; j++) {
                fprintf(f, "%f\n", solver->sol_vec[i*(N-1)+j]);
            }
	    fprintf(f, "%f\n", masa_eval_2d_exact_t(xmax, ymin+(i+1)*h));
        }
	for(int j = -1; j < N; j++) {
            fprintf(f, "%f\n", masa_eval_2d_exact_t(xmin+(j+1)*h, ymax));
        }
    }

    fclose(f);

    // output the analytical solution to the file
    f = fopen(solver->masa_file, "w");
    if(f == NULL)  {
        printf("Error opening file!\n");
        exit(1);
    }

    if(solver->dimension == 1) {
        masa_init("nick","heateq_1d_steady_const");
        masa_set_param("A_x", 10.0);
        masa_set_param("k_0", k);
        for(int i = 0; i < N+1; i++) {
            fprintf(f, "%f\n", masa_eval_1d_exact_t(xmin + i*h));
        }
    } else {
        masa_init("nick","heateq_2d_steady_const");
        masa_set_param("A_x", 10.0);
        masa_set_param("B_y", 10.0);
        masa_set_param("k_0", k);
        for(int i = 0; i < N+1; i++) {
            for(int j = 0; j < N+1; j++) {
                fprintf(f, "%f\n", masa_eval_2d_exact_t(xmin+j*h, ymin+i*h));
            }
        }
    }

    fclose(f);

    // verification
    if(solver->verify_mode == 0) {
	double err = 0.0;
	double temp = 0.0;
        if(solver->dimension == 1) {
            for(int i = 0; i < N-1; i++) {
                temp = solver->sol_vec[i] - masa_eval_1d_exact_t(xmin+(i+1)*h);
		err += temp*temp;
	    }
	    err = sqrt(err/(N+1));
        } else {
	    for(int i = 0; i < N-1; i++) {
	        for(int j = 0; j < N-1; j++) {
	            temp = solver->sol_vec[i*(N-1)+j] - masa_eval_2d_exact_t(xmin+(j+1)*h, ymin+(i+1)*h);
	            err += temp*temp;
		}
	    }
	    err = sqrt(err/((N+1)*(N+1)));
        }
	
	printf("** The verification mode is launched\n");
	printf("** Computing l2 error norm...\n");
	printf("--> l2 norm of the error between numerical solution and the analytic solution= %.12f\n\n\n", err);
    }

    // free the dynamically allocated memory
    int N_s = N - 1;
    if(solver->dimension == 2) {
	N_s = (N-1)*(N-1);
    }
    for(int i = 0; i < N_s; i++) {
        free(solver->matrix_value[i]);
        free(solver->matrix_index[i]);
    }
    free(solver->matrix_value);
    free(solver->matrix_index);
    free(solver->length);
    free(solver->rht_vec);
    free(solver->sol_vec);

    if(solver->output_mode == 1) {
        printf("[debug]: output - function end\n\n\n");
    }

    grvy_timer_end(__func__);

}
