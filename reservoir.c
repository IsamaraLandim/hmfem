#include "structs.h"
#include "variables.h"
#include "functions.h"
#include <stdio.h>
#include <stdlib.h>


void reservoir(Dim N, Perm *K,double *WI, double *KI, double **pc_old, IND *indcol, Bound *B, ELE *perm, grid h, FACE boundary, double viscosidade, double *pw, DimK *M, double *p)
{
int i,aux,dir,j,dirx,diry;

double *b, *pc_aux;

 //x = create_vector(N.t); /* x: inicial approximation for AMG and solution of pressure linear system - pressure */
 b = create_vector(N.t);
 pc_aux=create_vector(N.t);
 


for(j=0;j<number_karst;j++){
    if(M[j].pos==1){
        dir=M[j].sxi;
        for(i=0;i<M[j].s;i++){
        aux=pindex(dir,M[j].syi,M[j].szi,N);
        pc_aux[aux]=pc_old[i][j];
        dir++;  
    }
    }

    else if(M[j].pos==2){
        dir=M[j].syi;
        for(i=0;i<M[j].s;i++){
            aux=pindex(M[j].sxi,dir,M[j].szi,N);
            pc_aux[aux]=pc_old[i][j];
            dir++; 
        }
    }
    else if(M[j].pos==3){
        dir=M[j].szi;
        for(i=0;i<M[j].s;i++){
            aux=pindex(M[j].sxi,M[j].syi,dir,N);
            pc_aux[aux]=pc_old[i][j];
            dir++; 
        }
    }  
}



permeabilidade_arestas(K, perm, indcol, h, N, boundary,WI,KI, viscosidade); 
source_term(B, perm, h, N, boundary,pc_aux,pw, KI,WI,viscosidade,b);
solve_linsys(perm,indcol,b,p,N);

free(b);
free(pc_aux);
    
}
