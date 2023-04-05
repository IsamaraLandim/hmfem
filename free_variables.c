#include "structs.h"
#include "variables.h"
#include "functions.h"
#include <stdio.h>
#include <stdlib.h>

void libera(double **mat, int N, int M)
{
    // Free matrix memory
    int i;

    for(i=0; i<= N+1; i++){
        free(mat[i]);
    }
    free(mat);

    return;
}



