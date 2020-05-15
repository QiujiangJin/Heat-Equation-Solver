/**********************************************************************************************************
 * This fucntion calculates the discreted l2 norm of the error between the current and previous iteration
 * * *****************************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "laplace.h"
#include "function.h"

double error_l2(Laplace* solver, double* x_vec) {

    double res =0.0;
    int N = solver->nx - 1;
    if(solver->dimension == 2) {
	N = (N-1)*(N-1);
    }

    for(int i = 0; i < N; i++) {
	res += (solver->sol_vec[i] - x_vec[i])*(solver->sol_vec[i] - x_vec[i]);
    }
    
    return sqrt(res/N);

}
