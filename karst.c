#include "structs.h"
#include "variables.h"
#include "functions.h"
#include <stdio.h>
#include <stdlib.h>


void conduit(double **p_aux,double **cont_karst, double **KI_c, double viscosidade, DimK *M, double *hs,IND *indcol1, FACEK *BoundKarst, ELEK *keff, Perm *K, ELE *perm, double *p_old, double **pc)

{
int i,j;

double *x1,*f;

for(j=0;j<number_karst;j++){

x1= create_vector(M[j].s); 
f=create_vector(M[j].s);
permeabilidade_karst(K,keff,indcol1,hs,M,BoundKarst,viscosidade,KI_c,j);


source_term1d(keff, hs, M, BoundKarst,cont_karst,KI_c,p_aux,K,perm,p_old,j,f);


solve_linsys_1d(keff,indcol1,f,x1,M,j);
for(i=0;i<M[j].s;i++){
     pc[i][j]=x1[i];
}
free(f);
free(x1);

} 


     
//return(pressure);
//free(pressure);

  
 

}
