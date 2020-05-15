/****************************************************************************************
 * This is the header file of the declaration of all the functions used in this project
 * *************************************************************************************/

#include <laplace.h>
#ifndef FUNCTION_H
#define FUNCTION_H

// This funciton inputs all the data from input_file into the parameters in the solver
void input(Laplace* solver, const char* input_file);

// This function initializes all the other variables in the solver like matrix and vector
void initialize(Laplace* solver);

// This function builds the linear system in one dimensional case
void build_linear_system_1d(Laplace* solver);

// This function builds the linear system in two dimensional case
void build_linear_system_2d(Laplace* solver);

// This function solves the lienar system using Jacobi method or Gauss Seidel method
void solve_linear_system(Laplace* solver);

// This function solves the lienar system using PETSC
int solve_linear_system_PETSC(Laplace* solver);

// This function outputs the result to the output_file
void output(Laplace* solver);

// This fucntion calculates the 2 norm of b-Ax
double res_norm(Laplace* solver);

// This fucntion calculates the discreted l2 norm of the error between the current and previous iteration
double error_l2(Laplace* solver, double* x_vec);

// This function prints the matrix and vector of the linear system
void show_linear_system(Laplace* solver);

#endif

