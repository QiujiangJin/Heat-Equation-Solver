/**********************************************************
 * This is the header file of the type Laplace. It stores
 * all the parameters of the problem
 * *******************************************************/

#ifndef LAPLACE_H
#define LAPLACE_H

typedef struct Laplace {

    int dimension;
    double xmin;
    double xmax;
    double ymin;
    double ymax;
    int nx;
    int ny;
    
    // h is the distance between two neighbour points
    double h;
    
    double k;
    int fd_method;
    char* iter_method;
    int verify_mode;
    int output_mode;
    double eps;
    int max_iter;
    char* output_file;
    char* masa_file;

    // matrix_value stores the value of the non-zero element of the sparse matrix in each row
    double** matrix_value;
    
    // matrix_index stores the index of the non-zero element of the sparse matrix in each row
    int** matrix_index;
 
    // length stores the number of the non-zero element of the sparse matrix in each row
    int* length;
   
    // rht_vex stores the value of the right side vector 
    double* rht_vec;
    
    // sol_vex stores the value of the solution vector
    double* sol_vec;

} Laplace;

#endif
