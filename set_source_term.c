#include "structs.h"
#include "variables.h"
#include "functions.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


void source_term(Bound *B, ELE *perm, grid h, Dim N, FACE boundary, double *pc_aux, double *pw, double *KI, double *WI, double viscosidade,double *b)
{                                                        
/* Set source term vector: RHS of pressure linear system */

    int i, j, k;
    int s, l, interface;
   

    /* IF: in case of a complete Neumann boundary:
        Specifically here, a classical 1/4 of 5-spot problem.
        CHANGE values and location of source and sinks
    */


    for(interface=1;interface<7;interface++)
    {
        switch(interface)
        {
            case 1: // LEFT (L)

                if(boundary.l == 0)
                { // Dirichlet
                    for(j=0; j< N.y; j++)
                    {
                        for(k=0; k< N.z; k++){
                            s = pindex(0,j,k,N);
                            l = pyzindex(j,k,N);
                            b[s] = b[s]+ (perm[s].l)*B->l[l];
                        }
                    }
                }

                if(boundary.l == 1)
                { // Neumann 
                    for(j=0; j< N.y; j++)
                    {
                        for(k=0; k< N.z; k++){
                            s = pindex(0,j,k,N);
                            l = pyzindex(j,k,N);
                            b[s]= b[s]+ B->l[l]/h.x;
                        }
                    }
                }

                /* Case all Neumann boundary condition */
                if(boundary.l == 1 && boundary.r == 1 && boundary.u == 1 && boundary.d == 1 && boundary.t == 1 && boundary.b == 1)
                {

                    /* set pressure = 1.0 if Neumann problem at point (N.x-1, 0, N.z/2) */
                    s = pindex(0,N.y-1,N.z/2,N);
                    b[s]= b[s]+ (perm[s].l)*1.0;
                }

            break;

            case 2: // RIGHT (R) 

                if(boundary.r == 0)
                { // Dirichlet
                    for(j=0; j< N.y; j++)
                    {
                        for(k=0; k< N.z; k++){
                            s = pindex(N.x-1,j,k,N);
                            l = pyzindex(j,k,N);
                            b[s] =b[s]+ perm[s].r*B->r[l];
                        }
                    }
                }

                if(boundary.r == 1)
                { // Neumann 
                    for(j=0; j< N.y; j++)
                    {
                        for(k=0; k< N.z; k++){
                            s = pindex(N.x-1,j,k,N);
                            l = pyzindex(j,k,N);
                            b[s] = b[s]- B->r[l]/h.x;
                        }
                    }
                }

            break;

            case 3: // DOWN (D) 

                if(boundary.d == 0)
                { // Dirichlet
                    for(i=0; i< N.x; i++)
                    {
                        for(k=0; k< N.z; k++){
                            s = pindex(i,0,k,N);
                            l = pxzindex(i,k,N);
                            b[s] = b[s] + perm[s].d*B->d[l];
                        }
                    }
                }

                if(boundary.d == 1)
                { // Neumann 
                    for(i=0; i< N.x; i++)
                    {
                        for(k=0; k< N.z; k++){
                            s = pindex(i,0,k,N);
                            l = pxzindex(i,k,N);
                            b[s] = b[s] + B->d[l]/h.y;
                        }
                    }
                }

            break;

            case 4: // UPPER (U) 

                if(boundary.u == 0)
                { // Dirichlet
                    for(i=0; i< N.x; i++)
                    {
                        for(k=0; k< N.z; k++){
                            s = pindex(i,N.y-1,k,N);
                            l = pxzindex(i,k,N);
                            b[s] = b[s] + perm[s].u*B->u[l];
                        }
                    }
                }

                if(boundary.u == 1)
                { // Neumann 
                    for(i=0; i< N.x; i++)
                    {
                        for(k=0; k< N.z; k++){
                            s = pindex(i,N.y-1,k,N);
                            l = pxzindex(i,k,N);
                            b[s]= b[s] - B->u[l]/h.y;
                        }
                    }
                }

            break;

            case 5: // TOP (T) 

                if(boundary.t == 0)
                { // Dirichlet
                    for(i=0; i< N.x; i++)
                    {
                        for(j=0; j< N.y; j++){
                            s = pindex(i,j,N.z-1,N);
                            l = pxyindex(i,j,N);
                            b[s] = b[s] + perm[s].t*B->t[l];
                        }
                    }
                }

                if(boundary.t == 1)
                { // Neumann 
                    for(i=0; i< N.x; i++)
                    {
                        for(j=0; j< N.y; j++){
                            s = pindex(i,j,N.z-1,N);
                            l = pxyindex(i,j,N);
                            b[s] = b[s] - B->t[l]/h.z;
                        }
                    }
                }

            break;

            case 6: // BOTTOM (B) 

                if(boundary.b == 0)
                { // Dirichlet
                    for(i=0; i< N.x; i++)
                    {
                        for(j=0; j< N.y; j++){
                            s = pindex(i,j,0,N);
                            l = pxyindex(i,j,N);
                            b[s]= b[s] + perm[s].b*B->b[l];
                        }
                    }
                }

                if(boundary.b == 1)
                { // Neumann 
                    for(i=0; i< N.x; i++)
                    {
                        for(j=0; j< N.y; j++){
                            s = pindex(i,j,0,N);
                            l = pxyindex(i,j,N);
                            b[s] = b[s] + B->b[l]/h.z;
                        }
                    }
                }

            break;
        }
    }

double Ve;
Ve=h.x*h.y*h.z;
double a,c,d,e;

    for (s=0;s<N.t;s++){


               
b[s] =b[s]+(KI[s]*pc_aux[s])/(viscosidade*Ve)+(WI[s]*pw[s])/(viscosidade*Ve);
 

//printf("b (%d)=%e, a=%e, c=%e, d=%e, e=%e\n",s, b[s], a,c,d,e);
            
        
    }


}


void set_source_WI(Dim N)
{
/*
Read the absolute permeability from a file and set the permeability;
*/
    int err_chk, chk_scan; 
    double bb;
    FILE *fwi;
    int maxchar = 2048;
    char strg[maxchar]; 
    char *delim;
    int i, j, k, s, r, indx, indy, indz, s_wi; 
    int nx, ny, nz;
    double hx_wi, hy_wi, hz_wi;
    double posx, posy, posz;
    int ii,jj,kk,ss; 
    int iii,jjj,kkk,sss;
    
    for (k=0;k<N.z;k++){

        s_wi = pindex(ix0,jy0,k,N);
        wi[s_wi] = karst;
        IDwi[s_wi] = 1; 
        Nwi = Nwi + 1;
    }
}


/* CASE 1D */

void source_term1d(ELEK *keff, double *hs, DimK *M, FACEK *BoundKarst,double **cont_karst, double **KI_c, double **p_aux,Perm *K,ELE *perm, double *p_old, int j,double *f)                                                   
/* Set source term vector: RHS of pressure linear system */
{
    int i;
    int s, l, interface,saux;
    double c;
    
    for(i=0;i<M[j].s;i++)
    {
    f[i]=(KI_c[i][j]*p_aux[i][j])/viscosidade;
    }
    
// Elementos com brechas
double cR,uR,uL;
double R;

for(i=0;i<n_brechas[j];i++){
s=brechas[i][j];
cR=hm(K->S[s][j],K->S[s+1][j]);

if(ucR_old[i][j]<0){
    R=-C[i][j]*beta[i][j];
}
else{
    R=C[i][j]*beta[i][j];
}

uR=(cR/(viscosidade*hs[j]))*(pc_old[s][j]-pc_old[s+1][j]-rho*ucR_old[i][j]*ucR_old[i][j]*R);
uL=(cR/(viscosidade*hs[j]))*(pc_old[s+1][j]-pc_old[s][j]+rho*ucL_old[i][j]*ucL_old[i][j]*R);

ucR_old[i][j]=uR;
ucL_old[i][j]=uL;

f[s]=f[s]+(cR/(viscosidade*hs[j]))*(rho*uR*uR*R);
f[s+1]=f[s+1]-(cR/(viscosidade*hs[j]))*(rho*uL*uL*R);
   
}

    for(interface=1;interface<3;interface++)
    {
        switch(interface)
        {
            case 1: // LEFT (L)
                if(BoundKarst->l[j] == 0)
                { // Dirichlet
                f[0] = f[0] + (keff->l[0][j])*cont_karst[0][j];

                }
                if(BoundKarst->l[j] == 1)
                { // Neumann
                f[0] = f[0] +cont_karst[0][j]/(hs[j]);
                     
                }

                if(BoundKarst->l[j]==3)
                {//Robin nao linear
                double font;                
                
                if(M[j].pos==1){
                s = pindex(M[j].sxi-1,M[j].syi,M[j].szi,N);
                saux = pindex(M[j].sxi,M[j].syi,M[j].szi,N);
                double ur = h.x*perm[s].r*(p_old[s] - p_old[saux]);
                double lR=(K->X[s]*p_old[s]+K->X[saux]*p_old[saux])/(K->X[s]+K->X[saux]);
                    if(ur<0){
                       c=-C_contorno[0][j]*(cont_karst[0][j]-1);}
                    else{
                        c=C_contorno[0][j]*(cont_karst[0][j]-1);}
                    
                   
                font=-lR+rho*ur*ur*c;

                }

                else if(M[j].pos==2){
                s = pindex(M[j].sxi,M[j].syi-1,M[j].szi,N);
                saux = pindex(M[j].sxi,M[j].syi,M[j].szi,N);
                double uU = h.y*perm[s].u*(p_old[s] - p_old[saux]);
                double lU=(K->Y[s]*p_old[s]+K->Y[saux]*p_old[saux])/(K->Y[s]+K->Y[saux]);
                    if(uU<0){
                       c=-C_contorno[0][j]*(cont_karst[0][j]-1);}
                    else{
                        c=C_contorno[0][j]*(cont_karst[0][j]-1);}

                font=-lU+rho*uU*uU*c;
                }

                else if(M[j].pos==3){
                s = pindex(M[j].sxi,M[j].syi,M[j].szi-1,N);
                saux = pindex(M[j].sxi,M[j].syi,M[j].szi,N);
                double ut = h.z*perm[s].t*(p_old[s] - p_old[saux]);
                double lT=(K->Z[s]*p_old[s]+K->Z[saux]*p_old[saux])/(K->Z[s]+K->Z[saux]);
                    if(ut<0){
                      c=-C_contorno[0][j]*(cont_karst[0][j]-1);}
                    else{
                        c=C_contorno[0][j]*(cont_karst[0][j]-1);}

                font=-lT+rho*ut*ut*c;
                }
            

                f[0]=f[0]-2*K->S[0][j]/(viscosidade*hs[j]*hs[j])*font;

                }
                

                break;

            case 2: // RIGHT (R) 

                if(BoundKarst->r[j] == 0)                       
                { // Dirichlet
                    f[M[j].s-1] = f[M[j].s-1] + keff->r[M[j].s-1][j]*cont_karst[1][j];
                }

                if(BoundKarst->r[j] == 1)
                { // Neumann
                    f[M[j].s-1] = f[M[j].s-1] -cont_karst[1][j]/(hs[j]);
                     
                }

                                if(BoundKarst->r[j]==3)
                {
                double font;

                
                if(M[j].pos==1){
                s = pindex(M[j].sxi+M[j].s,M[j].syi,M[j].szi,N);
                saux = pindex(M[j].sxi+M[j].s-1,M[j].syi,M[j].szi,N);
                double ul = h.x*perm[s].l*(p_old[s] - p_old[saux]);
                double lL= (K->X[s]*p_old[s]+K->X[saux]*p_old[saux])/(K->X[s]+K->X[saux]);
                if(ul<0){
                    c=-C_contorno[1][j]*(cont_karst[1][j]-1);}
                else{
                    c=C_contorno[1][j]*(cont_karst[1][j]-1);}
                font=-lL+rho*ul*ul*c;
                }

                else if(M[j].pos==2){
                s = pindex(M[j].sxi,M[j].syi+M[j].s,M[j].szi,N);
                saux = pindex(M[j].sxi,M[j].syi+M[j].s-1,M[j].szi,N);
                double ud = h.y*perm[s].d*(p_old[s] - p_old[saux]);
                double ld= (K->Y[s]*p_old[s]+K->Y[saux]*p_old[saux])/(K->Y[s]+K->Y[saux]);
                if(ud<0){
                       c=-C_contorno[1][j]*(cont_karst[1][j]-1);}
                else{
                    c=C_contorno[1][j]*(cont_karst[1][j]-1);}
                font=-ld+rho*ud*ud*c;
                }

                            
                f[M[j].s-1]=f[M[j].s-1]-font*2*K->S[M[j].s-1][j]/(viscosidade*hs[j]*hs[j]);
                }

              
            break;          

}
}
}
