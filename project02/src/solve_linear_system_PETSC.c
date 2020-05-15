/*************************************************************************************
 * This function solves the lienar system using PETSC
 * **********************************************************************************/
#ifdef INCLUDE_PETSC

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <grvy.h>
#include "petsc.h"
#include "laplace.h"
#include "function.h"

int solve_linear_system_PETSC(Laplace* solver) {
    
    grvy_timer_begin(__func__);

    if(solver->output_mode == 1) {
	printf("[debug]: solve system - function begin\n\n\n");
    }

    printf("** Solving linear system...\n");

    // initialize
    PetscErrorCode ierr = PetscInitializeNoArguments(); CHKERRQ(ierr);

    Mat A;
    KSP Solver;
    PC Prec;
    Vec Rhs, Sol, Res;
    PetscInt N = solver->nx - 1;
    PetscReal norm;
    if(solver->dimension == 2) {
	N = N*N;
    }    

    // create and set values of matrix
    ierr = MatCreateSeqAIJ(PETSC_COMM_SELF, N, N, 9, solver->length, &A); CHKERRQ(ierr);
    for(PetscInt i = 0; i < N; i++) {
	for(PetscInt j = 0; j < solver->length[i]; j++) {
	    ierr = MatSetValue(A, i, solver->matrix_index[i][j], solver->matrix_value[i][j], INSERT_VALUES); CHKERRQ(ierr);
	}
    }	
    ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

    // create and set solver
    ierr = KSPCreate(PETSC_COMM_SELF, &Solver); CHKERRQ(ierr);
    ierr = KSPSetOperators(Solver, A, A); CHKERRQ(ierr);
    ierr = KSPSetType(Solver, KSPCGS); CHKERRQ(ierr);

    ierr = KSPGetPC(Solver, &Prec); CHKERRQ(ierr);
    ierr = PCSetType(Prec, PCJACOBI); CHKERRQ(ierr);

    // create and set values of right sided vector
    ierr = VecCreateSeqWithArray(PETSC_COMM_SELF, 1, N, solver->rht_vec, &Rhs); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(Rhs); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(Rhs); CHKERRQ(ierr);

    // create and set values of solution vector
    ierr = VecCreateSeqWithArray(PETSC_COMM_SELF, 1, N, solver->sol_vec, &Sol); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(Sol); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(Sol); CHKERRQ(ierr);
   
    // solve the system
    ierr = KSPSolve(Solver, Rhs, Sol); CHKERRQ(ierr);
    
    // calculate the residual norm
    ierr = VecDuplicate(Rhs, &Res); CHKERRQ(ierr);
    ierr = MatMult(A, Sol, Res); CHKERRQ(ierr);
    ierr = VecAXPY(Res, -1, Rhs); CHKERRQ(ierr);
    ierr = VecNorm(Res, NORM_2, &norm); CHKERRQ(ierr);
    ierr = PetscPrintf(MPI_COMM_WORLD, "residual norm: %e\n"); CHKERRQ(ierr); 

    ierr = VecGetArray(Sol, &solver->sol_vec); CHKERRQ(ierr);

    // free memory
    ierr = VecDestroy(&Res); CHKERRQ(ierr);
    ierr = KSPDestroy(&Solver); CHKERRQ(ierr);
    ierr = VecDestroy(&Rhs); CHKERRQ(ierr);
    ierr = VecDestroy(&Sol); CHKERRQ(ierr);

    // finalize
    ierr = PetscFinalize(); CHKERRQ(ierr);
    
    if(solver->output_mode == 1) {
        printf("[debug]: solve system - function end\n\n\n");
    }

    grvy_timer_end(__func__);

}

#endif
