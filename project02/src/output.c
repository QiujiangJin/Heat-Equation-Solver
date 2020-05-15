/********************************************************
 * This function outputs the result to the output_file
 * *****************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <grvy.h>
#include <masa.h>
#include <hdf5.h>
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
 
    hid_t file, dataset, dataspace;
    hsize_t dimsf[1];
    herr_t status;
    file = H5Fcreate(solver->output_file, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    if(solver->dimension == 1) {
	// case of 1 dimension
	dimsf[0] = N-1;
	dataspace = H5Screate_simple(1, dimsf, NULL);
	// store x coordinate
	double *cor_x = (double *)malloc((N-1)*sizeof(double));
	for(int i = 0; i < N-1; i++) {
	    cor_x[i] = xmin + (i+1)*h;
	}
	dataset = H5Dcreate2(file, "coordinate_x", H5T_NATIVE_DOUBLE, dataspace,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, cor_x);
	free(cor_x);
	// store numerical solution
	dataset = H5Dcreate2(file, "numerical_solution", H5T_NATIVE_DOUBLE, dataspace,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, solver->sol_vec);
	// store analytical solution
	masa_init("nick","heateq_1d_steady_const");
        masa_set_param("A_x", 10.0);
        masa_set_param("k_0", k);
	double *masa_dat = (double *)malloc((N-1)*sizeof(double));
        for(int i = 0; i < N-1; i++) {
            masa_dat[i] = masa_eval_1d_exact_t(xmin + (i+1)*h);
        }
        dataset = H5Dcreate2(file, "analytical_solution", H5T_NATIVE_DOUBLE, dataspace,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, masa_dat);
        free(masa_dat);

    } else {
	// case of 2 dimension
	dimsf[0] = (N-1)*(N-1);
	dataspace = H5Screate_simple(1, dimsf, NULL);
	// store x coordinate
        double *cor_x = (double *)malloc((N-1)*(N-1)*sizeof(double));
        for(int i = 0; i < N-1; i++) {
	    for(int j = 0; j < N-1; j++) {
		cor_x[i*(N-1)+j] = xmin + (j+1)*h;
	    }
        }
        dataset = H5Dcreate2(file, "coordinate_x", H5T_NATIVE_DOUBLE, dataspace,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, cor_x);
        free(cor_x);
	// store y coordinate
	double *cor_y = (double *)malloc((N-1)*(N-1)*sizeof(double));
        for(int i = 0; i < N-1; i++) {
            for(int j = 0; j < N-1; j++) {
                cor_y[i*(N-1)+j] = ymin + (i+1)*h;
            }
        }
        dataset = H5Dcreate2(file, "coordinate_y", H5T_NATIVE_DOUBLE, dataspace,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, cor_y);
        free(cor_y);
	// store numerical solution
        dataset = H5Dcreate2(file, "numerical_solution", H5T_NATIVE_DOUBLE, dataspace,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, solver->sol_vec);
	// store analytical solution
        masa_init("nick","heateq_2d_steady_const");
        masa_set_param("A_x", 10.0);
        masa_set_param("B_y", 10.0);
        masa_set_param("k_0", k);
	double *masa_dat = (double *)malloc((N-1)*(N-1)*sizeof(double));
        for(int i = 0; i < N-1; i++) {
	    for(int j = 0; j < N-1; j++) {
		masa_dat[i*(N-1)+j] = masa_eval_2d_exact_t(xmin + (j+1)*h, ymin + (i+1)*h);	
	    }
        }
        dataset = H5Dcreate2(file, "analytical_solution", H5T_NATIVE_DOUBLE, dataspace,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, masa_dat);
        free(masa_dat);

    }

    H5Sclose(dataspace);
    H5Dclose(dataset);
    H5Fclose(file);

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
