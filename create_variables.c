#include "structs.h"
#include "variables.h"
#include "functions.h"
#include <stdio.h>
#include <stdlib.h>

double** create_matrix(int n, int m){
    // Create null matrix (N.x+2)x(N.y+2).
    double** mat;
    int i, j, k;

    mat = (double**) malloc((n+2)*sizeof(double*));
    for(i=0; i<= n+1; i++){
            mat[i] = (double*) malloc((m+2)*sizeof(double));
    }

    for(i=0; i<= n; i++){
        for(j=0; j<= m+1; j++){
            mat[i][j] = 0;
        }
    }

    return mat;
}

int** create_matrix_int(int n, int m){
    // Create null matrix (N.x+2)x(N.y+2).
    int** mat;
    int i, j, k;

    mat = (int**) malloc((n+2)*sizeof(int*));
    for(i=0; i<= n+1; i++){
            mat[i] = (int*) malloc((m+2)*sizeof(int));
    }

    for(i=0; i<= n+1; i++){
        for(j=0; j<= m+1; j++){
            mat[i][j] = 0;
        }
    }

    return mat;
}

double* create_vector(int n){
    // Create null vector (n+2).
    double* vec;
    int i;

    vec = (double*) malloc(sizeof(double)*(n+2));

    for(i=0; i<= n+1; i++){
        vec[i] = 0;
    }
    
    return vec;
}

int* create_int(int n){
    int* vec;
    int i;

    vec = (int*) malloc(sizeof(double)*(n+2));

    for(i=0; i<= n+1; i++){
        vec[i] = 0;
    }
    
    return vec;
}




