/***********************************************
 * This fucntion calculates the 2 norm of b-Ax 
 * * ******************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "laplace.h"
#include "function.h"

double res_norm(Laplace* solver) {

    double res =0.0;
    int N = solver->nx;
    if(solver->dimension == 2) {
	N = (N-1)*(N-1) + 1;
    }

    for(int i = 0; i < N-1; i++) {
	double s = 0.0;
	for(int j = 0;j < solver->length[i]; j++) {
	    s += solver->matrix_value[i][j]*solver->sol_vec[solver->matrix_index[i][j]];
	}
	res += (solver->rht_vec[i] - s)*(solver->rht_vec[i] - s);
    }
    
    return sqrt(res);

}
