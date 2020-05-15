/***************************************************************************************
 * This funciton inputs all the data from input_file into the parameters in the solver
 * It also gives errors or warnings for inputs that are wrong
 * ************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <grvy.h>
#include "laplace.h"
#include "function.h"

void input(Laplace* solver, const char* input_file) {

    grvy_timer_begin(__func__);

    printf("\n\n\n**  Finite-difference based Heat Equation Solver (steady-state)\n\n\n");

    if(!grvy_input_fopen(input_file)) {
        printf("Input Error! Cannot open the input.dat. Please try again.\n");
	exit(1);
    }

    printf("**  Parsing runtime options from the file input.dat\n");
    printf("**  Runtime mesh settings:\n");

    if(grvy_input_fread_int("mesh/dimension", &(solver->dimension))) {
        printf("--> %-25s = %i\n","dimension", solver->dimension);
    }

    if(solver->dimension != 1 && solver->dimension != 2) {
       printf("Input Error! The dimension must be 1 or 2. Please try again.\n");
       exit(1);
    }

    if(grvy_input_fread_double("mesh/xmin", &(solver->xmin))) {
        printf("--> %-25s = %.12f\n","xmin", solver->xmin);
    }

    if(grvy_input_fread_double("mesh/xmax", &(solver->xmax))) {
        printf("--> %-25s = %.12f\n","xmax", solver->xmax);
    }

    if(solver->xmax < solver->xmin) {
        printf("Input Error! xmax cannot be smaller than xmin. Please try again.\n");
        exit(1);
    }
    
    if(solver->dimension == 2) {
        
	if(grvy_input_fread_double("mesh/ymin", &(solver->ymin))) {
            printf("--> %-25s = %.12f\n","ymin", solver->ymin);
        }

        if(grvy_input_fread_double("mesh/ymax", &(solver->ymax))) {
	    printf("--> %-25s = %.12f\n","ymax", solver->ymax);
        }

        if(solver->ymax < solver->ymin) {
            printf("Input Error! ymax cannot be smaller than ymin. Please try again.\n");
            exit(1);
        }

        if((solver->xmax - solver->xmin) != (solver->ymax - solver->ymin)) {
            printf("Input Error! (xmax - xmin) must be equal to (ymax -ymin). Please try again.\n");
            exit(1);
        }

    }

    if(grvy_input_fread_int("mesh/nx", &(solver->nx))) {
        printf("--> %-25s = %i\n","nx", solver->nx);
    }

    if(solver->dimension == 2) {

        if(grvy_input_fread_int("mesh/ny", &(solver->ny))) {
            printf("--> %-25s = %i\n","ny", solver->ny);
        }

        if(solver->nx != solver->ny) {
            printf("Input Error! nx must be equal to ny. Please try again.\n");
            exit(1);
        }

    }

    printf("\n\n** Runtime solver settings:\n");

    if(grvy_input_fread_int("solver/fd_method", &(solver->fd_method))) {
	if(solver->fd_method == 2) {
            printf("--> %-25s = %ind\n","finite difference order", solver->fd_method);
	} else {
            printf("--> %-25s = %ith\n","finite difference order", solver->fd_method);
	}
    }

    if(solver->fd_method != 2 && solver->fd_method != 4) {
       printf("Input Error! The fd_method must be 2 or 4. Please try again.\n");
       exit(1);
    }

    if(grvy_input_fread_char("solver/iter_method", &(solver->iter_method))) {
        printf("--> %-25s = %s\n","iteration method", solver->iter_method);
    }

    if(strcmp(solver->iter_method, "Jacobi") != 0 && strcmp(solver->iter_method, "Gauss-Seidel") != 0) {
       printf("Input Error! The iter_method must be \"Jacobi\" or \"Gauss-Seidel\". Please try again.\n");
       exit(1);
    }

    if(grvy_input_fread_int("solver/verify_mode", &(solver->verify_mode))) {
        if(solver->verify_mode == 0) {
            printf("--> %-25s = %s\n","verification mode", "verify");
        } else {
            printf("--> %-25s = %s\n","verification mode", "unverify");
        }
    }

    if(solver->verify_mode != 0 && solver->verify_mode != 1) {
       printf("Input Error! The verify_mode must be 0 or 1. Please try again.\n");
       exit(1);
    }

    if(grvy_input_fread_int("solver/output_mode", &(solver->output_mode))) {
        if(solver->output_mode == 0) {
            printf("--> %-25s = %s\n","output mode", "standard");
        } else {
            printf("--> %-25s = %s\n","output mode", "debug");
        }
    }

    if(solver->output_mode != 0 && solver->output_mode != 1) {
       printf("Input Error! The output_mode must be 0 or 1. Please try again.\n");
       exit(1);
    }

    if(grvy_input_fread_double("solver/k", &(solver->k))) {
        printf("--> %-25s = %.12f\n","thermal conductivity", solver->k);
    }

    if(grvy_input_fread_double("solver/eps", &(solver->eps))) {
        printf("--> %-25s = %.12f\n","convergence tolerance", solver->eps);
    }

    if(grvy_input_fread_int("solver/max_iter", &(solver->max_iter))) {
        printf("--> %-25s = %i\n","max iterations", solver->max_iter);
    }

    if(grvy_input_fread_char("solver/output_file", &(solver->output_file))) {
        printf("--> %-25s = %s\n","numerical solution file", solver->output_file);
    }

    if(grvy_input_fread_char("solver/masa_file", &(solver->masa_file))) {
        printf("--> %-25s = %s\n","analytical solution file", solver->masa_file);
    }

    if(solver->fd_method == 4 & strcmp(solver->iter_method, "Jacobi") == 0) {
	printf("Input Warning! The Jacobi method doesn't converge for the fourth order scheme. Please change your input.\n");
        exit(1);
    }

    grvy_input_fclose();

    printf("\n\n");

    grvy_timer_end(__func__);

}
