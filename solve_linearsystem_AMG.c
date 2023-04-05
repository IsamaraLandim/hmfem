#include "structs.h"
#include "variables.h"
#include "functions.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int ind(int i, int j, int k, int N, int M)
{
/* Counts, left -> right (i), down -> up (j), bottom -> top (k) */

    return j*N + k*N*M + i;
}

int pindex(int i, int j, int k, Dim N)
{
/* Counts, left -> right (i), down -> up (j), bottom -> top (k) */

    return j*N.x + k*N.xy + i;
}

int pxyindex(int i, int j, Dim N)
{
    return j*N.x + i;
}

int pxzindex(int i, int k, Dim N)
{
    return k*N.x + i;
}

int pyzindex(int j, int k, Dim N)
{
    return k*N.y + j;
}

void solve_linsys(ELE *perm, IND *indcol, double *b, double *x, Dim N)
{
/* Solve the pressure linear system by AMG toolbox. 
    store matrix by Compressed Row Storage.

    cnt[size=number of elements = N.t]: vector that stores the number of nonzero (NNZ) elements in each column; 

    hdr[size=2]: hdr[0] = number of nonzeros (NNZ) of matrix, hdr[1] = N.t;

    col[size=NNZ]: stores the column pindex of each nonzero element;
    ele[size=NNZ]: stores the actual value of each nonzero element;

*/

    int i,j,k;
    int l,s;
    int total_non_zero_elements = 0;
int *cnt, *hdr, *col;
double *ele;

    cnt = (int *) malloc ((N.t+1)*sizeof(int));

    for (s=0; s < N.t; s++)
    {
        cnt[s] = 0;
        if (indcol[s].d !=  -1){
          cnt[s] = cnt[s] + 1;
          total_non_zero_elements = total_non_zero_elements + 1;
        }
        if (indcol[s].l !=  -1){
          cnt[s] = cnt[s] + 1;
          total_non_zero_elements = total_non_zero_elements + 1;
        }
        if (indcol[s].c !=  -1){
          cnt[s] = cnt[s] + 1;
          total_non_zero_elements = total_non_zero_elements + 1;
        }
        if (indcol[s].r !=  -1){
          cnt[s] = cnt[s] + 1;
          total_non_zero_elements = total_non_zero_elements + 1;
        }
        if (indcol[s].u !=  -1){
          cnt[s] = cnt[s] + 1;
          total_non_zero_elements = total_non_zero_elements + 1;
        }
        if (indcol[s].t !=  -1){
          cnt[s] = cnt[s] + 1;
          total_non_zero_elements = total_non_zero_elements + 1;
        }
        if (indcol[s].b !=  -1){
          cnt[s] = cnt[s] + 1;
          total_non_zero_elements = total_non_zero_elements + 1;
        }
    }

    if(total_non_zero_elements == 0)
    {
      printf("ERROR: No nonzero elements. \n");
      exit(0);
    }

    hdr = (int *)malloc(2*sizeof(int));
    col = (int *)malloc(total_non_zero_elements*sizeof(int));
    ele = (double *)malloc(total_non_zero_elements*sizeof(double));

    hdr[0] = N.t;
    hdr[1] = total_non_zero_elements;
    
    l = 0;
    for (s=0; s < N.t; s++)
    {
        if (indcol[s].d !=  -1){
          col[l] = indcol[s].d;
          ele[l] = -perm[s].d;
          l ++;
        }
        if (indcol[s].l !=  -1){
          col[l] = indcol[s].l;
          ele[l] = -perm[s].l;
          l ++;
        }
        if (indcol[s].c !=  -1){
          col[l] = indcol[s].c;
          ele[l] = perm[s].c;
          l ++;
        }
        if (indcol[s].r !=  -1){
          col[l] = indcol[s].r;
          ele[l] = -perm[s].r;
          l ++;
        }
        if (indcol[s].u !=  -1){
          col[l] = indcol[s].u;
          ele[l] = -perm[s].u;
          l ++;
        }
        if (indcol[s].t !=  -1){
          col[l] = indcol[s].t;
          ele[l] = -perm[s].t;
          l ++;
        }
        if (indcol[s].b !=  -1){
          col[l] = indcol[s].b;
          ele[l] = -perm[s].b;
          l ++;
        }
    }

    /* compute(.) is a function that calls the AMG solver.
       AMG solver is written in C++.
       DO NOT EDIT THIS FUNCTION.
   */
    
    compute(hdr,cnt,col,ele,b,x);

    free(hdr);
    free(cnt);
    free(col);
    free(ele);

}


void prod_mat_vec(ELE *perm, IND *indcol, double *x, double *res, Dim N)
{
/* Calculates the product of the pressure linear system matrix by a vector x.
    Stores the anwser in vector res.
*/

    int i, j, k;
    int s;

    for (s=0; s < N.t; s++)
    {
        if (indcol[s].d !=  -1){
          res[s] = res[s] - perm[s].d*x[indcol[s].d];
        }
        if (indcol[s].l !=  -1){
          res[s] = res[s] - perm[s].l*x[indcol[s].l];
        }
        if (indcol[s].c !=  -1){
          res[s] = res[s] + perm[s].c*x[indcol[s].c];
        }
        if (indcol[s].r !=  -1){
          res[s] = res[s] - perm[s].r*x[indcol[s].r];
        }
        if (indcol[s].u !=  -1){
          res[s] = res[s] - perm[s].u*x[indcol[s].u];
        }
        if (indcol[s].t !=  -1){
          res[s] = res[s] - perm[s].t*x[indcol[s].t];
        }
        if (indcol[s].b !=  -1){
          res[s] = res[s] - perm[s].b*x[indcol[s].b];
        }
    }

}



/* CASE 1D */

void solve_linsys_1d(ELEK *keff, IND *indcol1, double *f, double *x1, DimK *M,int j)
{
/* Solve the pressure linear system by AMG toolbox. 
    store matrix by Compressed Row Storage.

    cnt[size=number of elements = N.t]: vector that stores the number of nonzero (NNZ) elements in each column; 

    hdr[size=2]: hdr[0] = number of nonzeros (NNZ) of matrix, hdr[1] = N.t;

    col[size=NNZ]: stores the column pindex of each nonzero element;
    ele[size=NNZ]: stores the actual value of each nonzero element;

*/
int *cnt1, *hdr1, *col1;
double *ele1;

    int i;
    int l,s;
    int total_non_zero_elements = 0;

    cnt1 = (int *) malloc ((Mmax+1)*sizeof(int));


    for (s=0; s < M[j].s; s++)
    {
        cnt1[s] = 0;
       
        if (indcol1[s].l !=  -1){
          cnt1[s] = cnt1[s] + 1;
          total_non_zero_elements = total_non_zero_elements + 1;
        }
       
        if (indcol1[s].r !=  -1){
          cnt1[s] = cnt1[s] + 1;
          total_non_zero_elements = total_non_zero_elements + 1;
        }
        
        if (indcol1[s].c !=  -1){
          cnt1[s] = cnt1[s] + 1;
          total_non_zero_elements = total_non_zero_elements + 1;
        }
    }

    if(total_non_zero_elements == 0)
    {
      printf("ERROR: No nonzero elements. \n");
      exit(0);
    }

    hdr1 = (int *)malloc(2*sizeof(int));
    col1 = (int *)malloc(total_non_zero_elements*sizeof(int));
    ele1 = (double *)malloc(total_non_zero_elements*sizeof(double));

    hdr1[0] = M[j].s;
    hdr1[1] = total_non_zero_elements;
    
    l = 0;
    for (s=0; s < M[j].s; s++)
    {
          if (indcol1[s].l !=  -1){
          col1[l] = indcol1[s].l;
          ele1[l] = -keff->l[s][j];
          l ++;
        }
        if (indcol1[s].c !=  -1){
          col1[l] = indcol1[s].c;
          ele1[l] = keff->c[s][j];
          l ++;
        }
        if (indcol1[s].r !=  -1){
          col1[l] = indcol1[s].r;
          ele1[l] = -keff->r[s][j];
          l ++;
        }
        
        
    }

    /* compute(.) is a function that calls the AMG solver.
       AMG solver is written in C++.
       DO NOT EDIT THIS FUNCTION.
   */
    
    compute(hdr1,cnt1,col1,ele1,f,x1);

     free(hdr1);
     free(cnt1);
     free(col1);
     free(ele1);
     //free(indcol1);

}


void prod_mat_vec1d(ELE *keff, IND *indcol1, double *x, double *res, DimK *M,int j)
{
/* Calculates the product of the pressure linear system matrix by a vector x.
    Stores the anwser in vector res.
*/

    int i;
    int s;

    for (s=0; s < M[j].s; s++)
    {
          if (indcol1[s].l !=  -1){
          res[s] = res[s] - keff[s].l*x[indcol1[s].l];
        }
        if (indcol1[s].c !=  -1){
          res[s] = res[s] + keff[s].c*x[indcol1[s].c];
        }
        if (indcol1[s].r !=  -1){
          res[s] = res[s] - keff[s].r*x[indcol1[s].r];
        }
        
    }

}




