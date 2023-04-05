#include "structs.h"
#include "variables.h"
#include "functions.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void create_boundary(double val_bound[7],Bound *B, Dim N,FACE boundary)
{
/* Set boundary values at each face.
    val_bound[1]: value of boundary (Dirichlet or Neumann) at Left boundary;
    val_bound[2]: value of boundary (Dirichlet or Neumann) at Right boundary;
    val_bound[3]: value of boundary (Dirichlet or Neumann) at Upper boundary;
    val_bound[4]: value of boundary (Dirichlet or Neumann) at Down boundary;
    val_bound[5]: value of boundary (Dirichlet or Neumann) at Top boundary;
    val_bound[6]: value of boundary (Dirichlet or Neumann) at Bottom boundary;
*/

    /* Create boundary.
        set_boundary(A,B,C,D);
        A: actual value of boundary (Dirichlet or Neumamm);
        B: B struct that will carry the fine vectors containing each boundary condition on each face;
        C: N struct - number of elements;
        D: numeration of faces
            1: Left
            2: Right
            3: Up 
            4: Down 
            5: Top 
            6: Bottom 
    */

    set_boundary(val_bound[1],B,N,1);
    set_boundary(val_bound[2],B,N,2);
    set_boundary(val_bound[3],B,N,3);
    set_boundary(val_bound[4],B,N,4);
    set_boundary(val_bound[5],B,N,5);
    set_boundary(val_bound[6],B,N,6);

} 

void set_boundary(double value, Bound *B, Dim N, int face)
{
/* Stores boundary values in struct B */

    int i, j, k;
    int s;

    switch(face)
    {
        case 1:
            // CASE Left (L):
            for(k=0;k < N.z;k++)
            {
                for(j=0;j < N.y;j++)
                {
                    s = pyzindex(j,k,N);
                    B->l[s] = value;
                }
            }

        break;
        
         case 2:
            // CASE Right (R):
            for(k=0;k < N.z;k++)
            {
                for(j=0;j < N.y;j++)
                {
                    s = pyzindex(j,k,N);
                    B->r[s] = value;
                }
            }

        break;
          case 3:
            // CASE Upper (U):
            for(k=0;k < N.z;k++)
            {
                for(i=0; i < N.x;i++)
                {
                    s = pxzindex(i,k,N);
                    B->u[s] = value;
                }
            }

        break;
           case 4:
            // case down (d):
            for(k=0;k < N.z;k++)
            {
                for(i=0; i < N.x;i++)
                {
                    s = pxzindex(i,k,N);
                    B->d[s] = value;
                }
            }

        break;
           case 5:
            // case top (t):
            for(j=0;j < N.y;j++)
            {
                for(i=0; i < N.x;i++)
                {
                    s = pxyindex(i,j,N);
                    B->t[s] = value;
                }
            }

        break;
           case 6:
            // case bottom (b):
            for(j=0;j < N.y;j++)
            {
                for(i=0; i < N.x;i++)
                {
                    s = pxyindex(i,j,N);
                    B->b[s] = value;
                }
            }

        break;
    }
}
