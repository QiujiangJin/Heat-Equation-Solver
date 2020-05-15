/*********************************************************************
 * This is the main function of this project. It creates a solver of 
 * type Laplace and use different functions to solve the problem
 * ******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <grvy.h>
#include "laplace.h"
#include "function.h" 

int main(int argc, char *argv[]) {

    if(argc != 2) {
        printf("Command Error!. The corrct command should be: \"./solver input.dat\". Please try again.\n");
        exit(1);
    }

    Laplace solver;

    grvy_timer_init("Heat Equation Solver");
 
    input(&solver, argv[1]);
    initialize(&solver);
    if(solver.dimension == 1) {
        build_linear_system_1d(&solver);
    } else {
        build_linear_system_2d(&solver);
    }
    if(strcmp(solver.iter_method, "PETSC") == 0) {
	solve_linear_system_PETSC(&solver);
    } else {
	solve_linear_system(&solver);
    }
    output(&solver);
    
    grvy_timer_finalize();
    grvy_timer_summarize();

    return 0;

}

