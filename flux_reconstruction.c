#include "structs.h"
#include "variables.h"
#include "functions.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


int fluxos(double *p, ELE *perm, ELE *v, Bound *B, grid h, Dim N)
{
/* Reconstruction of fluxes using:

   (uR + uL)/h.x + (uU + uD)/h.y + (uT + uB)/h.z = b

    where, in the interior elements we have,
    ui = (2Keff_(i,ii)/h.x)*(p - pii), i = R,L and ii = L,R;
    ui = (2Keff_(i,ii)/h.y)*(p - pii), i = U,D and ii = D,U;
    ui = (2Keff_(i,ii)/h.z)*(p - pii), i = T,B and ii = B,T;

    and in the boundary we have,
    ui = (2K/h.x)*(p - li), i = R,L;
    ui = (2K/h.y)*(p - li), i = U,D;
    ui = (2K/h.z)*(p - li), i = T,B;

*/

     /* CALCULATE THE FLUXES */ 
    int i, j, k, s, l, saux, interface;

    /* Internal elements */ 
    // Case UP (U):
	for(i=0; i< N.x; i++){
		for(j=0; j< (N.y-1); j++){
            for(k=0; k < N.z; k++){
                s = pindex(i,j,k,N);
                saux = pindex(i,j+1,k,N);
                v[s].u = h.y*perm[s].u*(p[s] - p[saux]);
            }
		}
	}

	// Case RIGHT (R):
	for(i=0; i< (N.x-1); i++){
		for(j=0; j< N.y; j++){
            for(k=0; k < N.z; k++){
                s = pindex(i,j,k,N);
                saux = pindex(i+1,j,k,N);
                v[s].r = h.x*perm[s].r*(p[s] - p[saux]); 
            }
		}
	}

	// Case DOWN (D):
	for(i=0; i< N.x; i++){
		for(j=1; j< N.y; j++){
            for(k=0; k < N.z; k++){           
                s = pindex(i,j,k,N);
                saux = pindex(i,j-1,k,N);
                v[s].d = h.y*perm[s].d*(p[s] - p[saux]); 
            }
		}
	}

	// Case LEFT (L):
	for(i=1; i< N.x; i++){
		for(j=0; j< N.y; j++){
            for(k=0; k < N.z; k++){
                s = pindex(i,j,k,N);
                saux = pindex(i-1,j,k,N);
                v[s].l = h.x*perm[s].l*(p[s] - p[saux]); 
            }
		}
	}


	// Case TOP (T):
	for(i=0; i< N.x; i++){
		for(j=0; j< N.y; j++){
            for(k=0; k < (N.z-1); k++){
                s = pindex(i,j,k,N);
                saux = pindex(i,j,k+1,N);
                v[s].t = h.z*perm[s].t*(p[s] - p[saux]); 
            }
		}
	}


	// Case BOTTOM (B):
	for(i=0; i< N.x; i++){
		for(j=0; j< N.y; j++){
            for(k=1; k < N.z; k++){
                s = pindex(i,j,k,N);
                saux = pindex(i,j,k-1,N);
                v[s].b = h.z*perm[s].b*(p[s] - p[saux]); 
            }
		}
	}

    /* Boundary elements */ 

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
                            v[s].l = h.x*perm[s].l*(p[s] - B->l[l]);
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
                        v[s].l = -B->l[l];
                        }
                    }
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
                            v[s].r = h.x*perm[s].r*(p[s] - B->r[l]);
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
                            v[s].r = B->r[l];
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
                            v[s].d = h.y*perm[s].d*(p[s] - B->d[l]);
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
                            v[s].d = -B->d[l];
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
                            v[s].u = h.y*perm[s].u*(p[s] - B->u[l]);
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
                            v[s].u = B->u[l];
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
                            v[s].t = h.z*perm[s].t*(p[s] - B->t[l]);
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
                            v[s].t = B->t[l];
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
                            v[s].b = h.z*perm[s].b*(p[s] - B->b[l]);
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
                            v[s].b = -B->b[l];
                        }
                    }
                }

            break;
        }
    }

    divvel_pontual(v,h,N);
    
    return 0;
}

void divvel_pontual(ELE *v, grid h, Dim N)
{
// Calculate the divergent ( div v = 0) max and min
    int i, j, k, s;
    double divvel;
    double maxx, minn;

    maxx = 0.0;
    minn = 100.0;

    for(i=0; i< N.x; i++)
    {
        for(j=0; j< N.y; j++)
        {
            for(k=0; k< N.z; k++)
            {
                s = pindex(i,j,k,N);

                divvel = v[s].r/h.x + v[s].l/h.x + v[s].d/h.y + v[s].u/h.y + v[s].b/h.z + v[s].t/h.z;

                if(fabs(divvel) > maxx)
                {
                    maxx = fabs(divvel);
                }

                if(fabs(divvel) < minn)
                {
                    minn = fabs(divvel);
                }

            }
        }
    }

    maxx = fabs(maxx);
    minn = fabs(minn);

    printf("Divv_min = %e, Divv_max = %e\n", minn,maxx);
}

//case 1d
int fluxo1d(double **x, ELEK *keff, ELEK *uc, double *hs, DimK *M, double **cont_karst, FACEK *BoundKarst, Perm *K, ELE *perm, double *p_old)
{
    /* CALCULATE THE FLUXES */ 
    int i, s, l,saux, interface,j;
    double c,cR;

    /* Internal elements */ 
   	// Case RIGHT (R):
    for(j=0;j<number_karst;j++){
        for(i=0; i<(M[j].s-1); i++){
        s = i;
        saux = i+1;
        uc->r[s][j] = hs[j]*keff->r[s][j]*(x[s][j] - x[saux][j]);
        }

        // Case LEFT (L):
        for(i=1; i< M[j].s; i++){
        s = i;
        saux = i-1;
        uc->l[s][j] = hs[j]*keff->l[s][j]*(x[s][j] - x[saux][j]); 
        }


         //ELEMENTOS COM BRECHAS
        
        for(i=0;i<n_brechas[j];i++){
            if(ucR_old[i][j]<0){
        c=-C[i][j]*beta[i][j];
        }
        else{
        c=C[i][j]*beta[i][j];
        }


        s=brechas[i][j];
        cR=hm(K->S[s][j], K->S[s+1][j]);
        uc->r[s][j]=(cR/(viscosidade*hs[j]))*(x[s][j]-x[s+1][j]-0.5*rho*ucR_old[i][j]*ucR_old[i][j]*c);
        uc->l[s+1][j]=(cR/(viscosidade*hs[j]))*(x[s+1][j]-x[s][j]+0.5*rho*ucL_old[i][j]*ucL_old[i][j]*c);
        }

    /* Boundary elements */ 

    for(interface=1;interface<3;interface++)
    {
        switch(interface)
        {
            case 1: // LEFT (L)

                if(BoundKarst->l[j] == 0)
                { // Dirichlet
                    uc->l[0][j] = hs[j]*keff->l[0][j]*(x[0][j] - cont_karst[0][j]);
                    
                }

                if(BoundKarst->l[j] == 1)
                { // neumann
                    uc->l[0][j] = -cont_karst[1][j];
                    
                }
               
                if(BoundKarst->l[j]==3)
                {

                double c;               
                if(M[j].pos==1){
                s = pindex(M[j].sxi-1,M[j].syi,M[j].szi,N);
                saux = pindex(M[j].sxi,M[j].syi,M[j].szi,N);
                double ur = h.x*perm[s].r*(p_old[s] - p_old[saux]);
                double lR=(K->X[s]*p_old[s]+K->X[saux]*p_old[saux])/(K->X[s]+K->X[saux]);
                    if(ur<0){
                        c=-C_contorno[0][j]*(cont_karst[0][j]-1);}
                    else{
                        c=C_contorno[0][j]*(cont_karst[0][j]-1);}
                uc->l[0][j]=(2*K->S[0][j])/(viscosidade*hs[j])*(x[0][j]-lR+rho*ur*ur*c);
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
                uc->l[0][j]=(2*K->S[0][j])/(viscosidade*hs[j])*(x[0][j]-lU+rho*uU*uU*c);
                }

                else if(M[j].pos==3){
                s = pindex(M[j].sxi,M[j].syi,M[j].szi-1,N);
                saux = pindex(M[j].sxi,M[j].syi,M[j].szi,N);
                double ut = h.z*perm[s].t*(p_old[s] - p_old[saux]);
                double lt=(K->Z[s]*p_old[s]+K->Z[saux]*p_old[saux])/(K->Z[s]+K->Z[saux]);
                    if(ut<0){
                       c=-C_contorno[0][j]*(cont_karst[0][j]-1);}
                    else{
                        c=C_contorno[0][j]*(cont_karst[0][j]-1);}
                uc->l[0][j]=(2*K->S[0][j])/(viscosidade*hs[j])*(x[0][j]-lt+rho*ut*ut*c);
                }              
            
            }

            break;

            case 2: // RIGHT (R) 

                if(BoundKarst->r[j] == 0)
                { // Dirichlet
                    s = i;
                    l = i+1;
                    uc->r[M[j].s-1][j] = hs[j]*keff->r[M[j].s-1][j]*(x[M[j].s-1][j] - cont_karst[1][j]);
                  
                }
                if(BoundKarst->r[j] == 1)
                { // neumann
                    uc->r[M[j].s-1][j] =  cont_karst[1][j];
                    
                }

                
                if(BoundKarst->r[j]==3)
                {

                double c;

            
                if(M[j].pos==1){
                s = pindex(M[j].sxi+M[j].s,M[j].syi,M[j].szi,N);
                saux = pindex(M[j].sxi+M[j].s-1,M[j].syi,M[j].szi,N);
                double ul=h.x*perm[s].l*(p_old[s]-p_old[saux]);
                double lL=(K->X[s]*p_old[s]+K->X[saux]*p_old[saux])/(K->X[s]+K->X[saux]);
                if(ul<0){
                    c=-C_contorno[1][j]*(cont_karst[1][j]-1);}
                else{
                    c=C_contorno[1][j]*(cont_karst[1][j]-1);}
                uc->r[M[j].s-1][j]=(2*K->S[M[j].s-1][j])/(viscosidade*hs[j])*(x[M[j].s-1][j]-lL+rho*ul*ul*c);
                }

                else if(M[j].pos==2){
                s = pindex(M[j].sxi,M[j].syi+M[j].s,M[j].szi,N);
                saux = pindex(M[j].sxi,M[j].syi+M[j].s-1,M[j].szi,N);
                double ud=h.y*perm[s].d*(p_old[s]-p_old[saux]);
                double ld=(K->Y[s]*p_old[s]+K->Y[saux]*p_old[saux])/(K->Y[s]+K->Y[saux]);
                if(ud<0){
                    c=-C_contorno[1][j]*(cont_karst[1][j]-1);}
                else{
                    c=C_contorno[1][j]*(cont_karst[1][j]-1);}
                uc->r[M[j].s-1][j]=(2*K->S[M[j].s-1][j])/(viscosidade*hs[j])*(x[M[j].s-1][j]-ld+rho*ud*ud*c);
                }

                else if(M[j].pos==3){
                s = pindex(M[j].sxi,M[j].syi,M[j].szi+M[j].s,N);
                saux = pindex(M[j].sxi,M[j].syi,M[j].szi+M[j].s-1,N);
                double ub=h.z*perm[s].b*(p_old[s]-p_old[saux]);
                double lb=(K->Z[s]*p_old[s]+K->Z[saux]*p_old[saux])/(K->Z[s]+K->Z[saux]);
                if(ub<0){
                    c=-C_contorno[1][j]*(cont_karst[1][j]-1);}
                else{
                    c=C_contorno[1][j]*(cont_karst[1][j]-1);}
                uc->r[M[j].s-1][j]=(2*K->S[M[j].s-1][j])/(viscosidade*hs[j])*(x[M[j].s-1][j]-lb+rho*ub*ub*c);
                } 

                
            }             


            
            break;           

        }
    }

}    
    return 0;
}



void divvel_pontual_1d(ELE *uc, grid h, Dim N)
{
// Calculate the divergent ( div v = 0) max and min
    int i, s;
    double divvel;
    double maxx, minn;

    maxx = 0.0;
    minn = 100.0;

    for(i=0; i< N.t; i++)
    {
                s = i;

                divvel = uc[s].r/h.x + uc[s].l/h.x;

                if(fabs(divvel) > maxx)
                {
                    maxx = fabs(divvel);
                }

                if(fabs(divvel) < minn)
                {
                    minn = fabs(divvel);
                }

    }

    maxx = fabs(maxx);
    minn = fabs(minn);

    printf("Divv_min = %e, Divv_max = %e\n", minn,maxx);
}




