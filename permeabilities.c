#include "structs.h"
#include "variables.h"
#include "functions.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


void set_perm_faces(Perm *K, grid I, grid h, Dim N,int Mmax, double *Kc)
{
/* Set type of permeability
    P->permeability = 1: homogeneous permeability
    P->permeability = 2: read permeability of a file: SPE 3D
    P->permeability = 3: Gaussian permeability tensor 
    P->permeability = 4: read permeability of a file: SPE 2D layer 
*/

    if((K->permeability > 7) || (K->permeability <= 0))
    {
        printf("ERROR! Type of permeability can be 1, 2 or 3, and P->permeability = %d \n", K->permeability);
        printf("Exit run... \n");
        exit(0);
    }

    /* Create permeability */

    switch(K->permeability){

        case 1:
            perm_homogeneous(K,N); 
        break;

        case 2:
            read_set_perm_field(K,N,h,I); // SPE
        break;

        case 3:
            gaussian_permfield(K, N, h, I);
        break;
        case 4:
            read_layer_perm_field(K,N,h,I); // SPE layer
        break;
    }

    /* Calculate the permeability in the conduit */ 
int j,i;
    for(j=0;j<number_karst;j++){
	    for(i=0;i<Mmax;i++){
		    K->S[i][j]=Kc[j];
	    }
    }

   

}

int read_layer_perm_field(Perm *K, Dim N, grid h, grid I)
{
/*
Reads and Construct the heterogeneous permeability
*/

    FILE *fperm;		   /* file pointer for perm field in x direction*/
    double bb;	   /* aux. varibles */

    int i, j, k, s, r, indx, indy, indz, s_geo; 
    int nx_geo, ny_geo, nz_geo, nt_geo;
    double hx_geo, hy_geo, hz_geo, fact_x, fact_y, fact_z;
    double Lx,Ly,Lz;
    double posx, posy, posz;
    int ii,jj,kk,ss; 
    int iii,jjj,kkk,sss;
    int chk_scan;

    Lx = 1200; 
    Ly = 2200;
    Lz = 170;
    nx_geo = 60;
    ny_geo = 220;
    nz_geo = 85;
    
    nt_geo = nx_geo*ny_geo*1;

    Perm *P_geo = (Perm*)malloc(sizeof(Perm));
    P_geo->X = (double*)malloc((nt_geo)*sizeof(double));
    P_geo->Y = (double*)malloc((nt_geo)*sizeof(double));
    P_geo->Z = (double*)malloc((nt_geo)*sizeof(double));

    /* Set geological parameters */
    hx_geo = I.x/(double)(nx_geo);
    hy_geo = I.y/(double)(ny_geo);

    /* Read permeability file */
    
    fperm = fopen(fn,"r");

    if (fperm == NULL)
        exit(EXIT_FAILURE);
    

    for (k=nz_geo-1;k>=0;k--){
    	for (j=0;j<ny_geo;j++){
        	for (i=0;i<nx_geo;i++){

                chk_scan = fscanf(fperm,"%lf",&bb);

                if (chk_scan != 1)
                    exit(EXIT_FAILURE);
                s = ind(nx_geo-i-1,ny_geo-j-1,0,nx_geo,ny_geo);
                
                if (k == layer) 
                    P_geo->X[s] = bb; 
            }
        }
    }

    for (k=nz_geo-1;k>=0;k--){
    	for (j=0;j<ny_geo;j++){
        	for (i=0;i<nx_geo;i++){
            
                chk_scan = fscanf(fperm,"%lf",&bb);

                if (chk_scan != 1)
                    exit(EXIT_FAILURE);

                s = ind(nx_geo-i-1,ny_geo-j-1,0,nx_geo,ny_geo);
                
                if (k == layer) 
                    P_geo->Y[s] = bb; 
            }
        }
    }

    for (k=nz_geo-1;k>=0;k--){
    	for (j=0;j<ny_geo;j++){
        	for (i=0;i<nx_geo;i++){
            
                chk_scan = fscanf(fperm,"%lf",&bb);

                if (chk_scan != 1)
                    exit(EXIT_FAILURE);

                s = ind(nx_geo-i-1,ny_geo-j-1,0,nx_geo,ny_geo);
                
                if (k == layer) 
                    P_geo->Z[s] = 1e-8; //bb; 
            }
        }
    }

    fclose(fperm);

    for (ii=0;ii<N.x;ii++){
        for (jj=0;jj<N.y;jj++){
            for (kk=0;kk<N.z;kk++){

                posx = (ii + 0.5)*h.x;
                posy = (jj + 0.5)*h.y;
                posz = (kk + 0.5)*h.z;

                indx = floor(posx/hx_geo);
                indy = floor(posy/hy_geo);
                indz = 0;
                    
                s_geo = ind(indx,indy,indz,nx_geo,ny_geo);

                ss = pindex(ii,jj,0,N);

                K->X[ss] = P_geo->X[s_geo]; 
                K->Y[ss] = P_geo->Y[s_geo]; 
                K->Z[ss] = P_geo->Z[s_geo]; 

            }
        }
    }

    if(P_geo !=NULL){
        free(P_geo->X);
        free(P_geo->Y);
        free(P_geo->Z);
        free(P_geo);
    }


    return 0;
}

double hm(double var1, double var2)
{
/* Calculates the harmonic mean */

    return 2*var1*var2/(var1 + var2);
}


int permeabilidade_arestas(Perm *K, ELE *perm, IND *indcol, grid h, Dim N, FACE boundary, double *WI, double *KI, double viscosidade)
{
/* Calculate the harmonic mean of the permeabilities on the edges of the elements.
    perm[s]: is an array of structs for Keff that travels each element using the pindex() function
    (go to function file for more information).
    The struct contains each edge + a variable on the "center" that corresponds to the sum
    of the keff in each edge.

    indcol[s]: is an array of structs that keeps an pindex for to indicate the location of
    the correspondent value of perm[s] in the pressure linear system global matrix.
    Used in the solution of the linear system by the AMG solver.
    indcol[s] = -1 when keff is on the boundary, and do not count for the matrix construction
    (it will go to the RHS).

*/ 
    int i, j, k, s, saux;
    int interface;

    /* Internal elements */ 
    /* Case UP (U) */ 
	for(i=0; i< N.x; i++){
		for(j=0; j< (N.y-1); j++){
            for(k=0; k < N.z; k++){
                s = pindex(i,j,k,N);
                saux = pindex(i,j+1,k,N);
                perm[s].u = hm(K->Y[s],K->Y[saux]);
                perm[s].u = perm[s].u/(h.y*h.y*viscosidade); // ok
                indcol[s].u = s + N.x;
            }
		}
	}

	/* Case RIGHT (R) */ 
	for(i=0; i< (N.x-1); i++){
		for(j=0; j< N.y; j++){
            for(k=0; k < N.z; k++){
                s = pindex(i,j,k,N);
                saux = pindex(i+1,j,k,N);
                perm[s].r = hm(K->X[s],K->X[saux]); 
                perm[s].r = perm[s].r/(h.x*h.x*viscosidade);// ok
                indcol[s].r = s + 1;
            }
		}
	}

	/* Case DOWN (D) */ 
	for(i=0; i< N.x; i++){
		for(j=1; j< N.y; j++){
            for(k=0; k < N.z; k++){           
                s = pindex(i,j,k,N);
                saux = pindex(i,j-1,k,N);
                perm[s].d = hm(K->Y[s],K->Y[saux]); 
                perm[s].d = perm[s].d/(h.y*h.y*viscosidade); // ok
                indcol[s].d = s - N.x;
            }
		}
	}

	/* Case LEFT (L) */ 
	for(i=1; i< N.x; i++){
		for(j=0; j< N.y; j++){
            for(k=0; k < N.z; k++){
                s = pindex(i,j,k,N);
                saux = pindex(i-1,j,k,N);
                perm[s].l = hm(K->X[s],K->X[saux]); 
                perm[s].l = perm[s].l/(h.x*h.x*viscosidade);// ok
                indcol[s].l = s - 1;
            }
		}
	}


	/* Case TOP (T) */ 
	for(i=0; i< N.x; i++){
		for(j=0; j< N.y; j++){
            for(k=0; k < (N.z-1); k++){
                s = pindex(i,j,k,N);
                saux = pindex(i,j,k+1,N);
                perm[s].t = hm(K->Z[s],K->Z[saux]); 
                perm[s].t = perm[s].t/(h.z*h.z*viscosidade);// ok
                indcol[s].t = s + N.xy;
            }
		}
	}


	/* Case BOTTOM (B) */ 
	for(i=0; i< N.x; i++){
		for(j=0; j< N.y; j++){
            for(k=1; k < N.z; k++){
                s = pindex(i,j,k,N);
                saux = pindex(i,j,k-1,N);
                perm[s].b = hm(K->Z[s],K->Z[saux]); 
                perm[s].b = perm[s].b/(h.z*h.z*viscosidade);// ok
                indcol[s].b = s - N.xy;
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
                            perm[s].l = 2.0*K->X[s]/(h.x*h.x*viscosidade);
                            indcol[s].l = -1;
                        }
                    }
                }

                if(boundary.l == 1)
                { // Neumann 
                    for(j=0; j< N.y; j++)
                    {
                        for(k=0; k< N.z; k++){
                        s = pindex(0,j,k,N);
                        perm[s].l = 0;
                        indcol[s].l = -1;
                        }
                    }
                }

                /* Case all Neumann boundary condition */
                if(boundary.l == 1 && boundary.r == 1 && boundary.u == 1 && boundary.d == 1 && boundary.t == 1 && boundary.b == 1)
                {
                    s = pindex(0,N.y-1,N.z/2,N);
                    perm[s].l = 2.0*K->X[s]/(h.x*h.x*viscosidade);
                    indcol[s].l = -1;
                }

            break;

            case 2: // RIGHT (R) 

                if(boundary.r == 0)
                { // Dirichlet
                    for(j=0; j< N.y; j++)
                    {
                        for(k=0; k< N.z; k++){
                            s = pindex(N.x-1,j,k,N);
                            perm[s].r = 2.0*K->X[s]/(h.x*h.x*viscosidade);
                            indcol[s].r = -1;
                        }
                    }
                }

                if(boundary.r == 1)
                { // Neumann 
                    for(j=0; j< N.y; j++)
                    {
                        for(k=0; k< N.z; k++){
                            s = pindex(N.x-1,j,k,N);
                            perm[s].r = 0;
                            indcol[s].r = -1;
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
                            perm[s].d = 2.0*K->Y[s]/(h.y*h.y*viscosidade);
                            indcol[s].d = -1;
                        }
                    }
                }

                if(boundary.d == 1)
                { // Neumann 
                    for(i=0; i< N.x; i++)
                    {
                        for(k=0; k< N.z; k++){
                            s = pindex(i,0,k,N);
                            perm[s].d = 0;
                            indcol[s].d = -1;
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
                            perm[s].u = 2.0*K->Y[s]/(h.y*h.y*viscosidade);
                            indcol[s].u = -1;
                        }
                    }
                }

                if(boundary.u == 1)
                { // Neumann 
                    for(i=0; i< N.x; i++)
                    {
                        for(k=0; k< N.z; k++){
                            s = pindex(i,N.y-1,k,N);
                            perm[s].u = 0;
                            indcol[s].u = -1;
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
                            perm[s].t = 2.0*K->Z[s]/(h.z*h.z*viscosidade);
                            indcol[s].t = -1;
                        }
                    }
                }

                if(boundary.t == 1)
                { // Neumann 
                    for(i=0; i< N.x; i++)
                    {
                        for(j=0; j< N.y; j++){
                            s = pindex(i,j,N.z-1,N);
                            perm[s].t = 0;
                            indcol[s].t = -1;
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
                            perm[s].b = 2.0*K->Z[s]/(h.z*h.z*viscosidade);
                            indcol[s].b = -1;
                        }
                    }
                }

                if(boundary.b == 1)
                { // Neumann 
                    for(i=0; i< N.x; i++)
                    {
                        for(j=0; j< N.y; j++){
                            s = pindex(i,j,0,N);
                            perm[s].b = 0;
                            indcol[s].b = -1;
                        }
                    }
                }

            break;
        }
    }


	/* Case CENTER (C): sum of the keff in each edge */ 
	for(i=0; i< N.x; i++){
		for(j=0; j< N.y; j++){
            for(k=0; k< N.z; k++){
                s = pindex(i,j,k,N);
                perm[s].c =perm[s].l + perm[s].r + perm[s].u + perm[s].d + perm[s].t + perm[s].b+KI[s]/(viscosidade*h.x*h.y*h.z)+WI[s]/(viscosidade*h.x*h.y*h.z); 
                indcol[s].c = s;
            }
		}
	}
    
    return 0;
}

int read_set_perm_field(Perm *K, Dim N, grid h, grid I)
{
/*
Reads and Construct the heterogeneous permeability
*/

    FILE *fperm;		   /* file pointer for perm field in x direction*/
    double bb;	   /* aux. varibles */

    int i, j, k, s, r, indx, indy, indz, s_geo; 
    int nx_geo, ny_geo, nz_geo, nt_geo;
    double hx_geo, hy_geo, hz_geo, fact_x, fact_y, fact_z;
    double Lx,Ly,Lz;
    double posx, posy, posz;
    int ii,jj,kk,ss; 
    int iii,jjj,kkk,sss;
    int chk_scan;

    Lx = 1200; 
    Ly = 2200;
    Lz = 170;
    nx_geo = 60;
    ny_geo = 220;
    nz_geo = 85;
    
    //if (lay == 1)
    //nt_geo = nx_geo*ny_geo*1;
    //else 
    nt_geo = nx_geo*ny_geo*nz_geo;

    Perm *P_geo = (Perm*)malloc(sizeof(Perm));
    P_geo->X = (double*)malloc((nt_geo)*sizeof(double));
    P_geo->Y = (double*)malloc((nt_geo)*sizeof(double));
    P_geo->Z = (double*)malloc((nt_geo)*sizeof(double));


    /* Set geological parameters */
    hx_geo = I.x/(double)(nx_geo);
    hy_geo = I.y/(double)(ny_geo);
    hz_geo = I.z/(double)(nz_geo);

    /* Read permeability file */
    
    fperm = fopen(fn,"r");

    if (fperm == NULL)
        exit(EXIT_FAILURE);
    

    for (k=nz_geo-1;k>=0;k--){
    	for (j=0;j<ny_geo;j++){
        	for (i=0;i<nx_geo;i++){

                chk_scan = fscanf(fperm,"%lf",&bb);

                if (chk_scan != 1)
                    exit(EXIT_FAILURE);
                /*
                if (lay == 1 || k==layer){
                s = ind(nx_geo-i-1,ny_geo-j-1,0,nx_geo,ny_geo);
                
                P_geo->X[s] = bb; 
               	}else{*/
		s = ind(nx_geo-i-1,ny_geo-j-1,k,nx_geo,ny_geo);
                
                P_geo->X[s] = bb;

		//} 
            }
        }
    }

    for (k=nz_geo-1;k>=0;k--){
    	for (j=0;j<ny_geo;j++){
        	for (i=0;i<nx_geo;i++){
            
                chk_scan = fscanf(fperm,"%lf",&bb);

                if (chk_scan != 1)
                    exit(EXIT_FAILURE);

                /*
                if (lay == 1 || k==layer){
                s = ind(nx_geo-i-1,ny_geo-j-1,0,nx_geo,ny_geo);
                
                P_geo->Y[s] = bb; 
               	}else{*/
		s = ind(nx_geo-i-1,ny_geo-j-1,k,nx_geo,ny_geo);
                
                P_geo->Y[s] = bb;
		//} 
            }
        }
    }

    for (k=nz_geo-1;k>=0;k--){
    	for (j=0;j<ny_geo;j++){
        	for (i=0;i<nx_geo;i++){
            
                chk_scan = fscanf(fperm,"%lf",&bb);

                if (chk_scan != 1)
                    exit(EXIT_FAILURE);
                /*
                if (lay == 1 || k==layer){
                //s = ind(nx_geo-i-1,ny_geo-j-1,0,nx_geo,ny_geo);
                
                //P_geo->Z[s] = bb; 
               	}else{*/
                s = ind(nx_geo-i-1,ny_geo-j-1,k,nx_geo,ny_geo);
                
                P_geo->Z[s] = bb;
 		//} 
            }
        }
    }

    fclose(fperm);

    for (ii=0;ii<N.x;ii++){
        for (jj=0;jj<N.y;jj++){
            for (kk=0;kk<N.z;kk++){
                /*
                if (lay == 1 || k==layer){
                posx = (ii + 0.5)*h.x;
                posy = (jj + 0.5)*h.y;
                //posz = (kk + 0.5)*h.z;

                indx = floor(posx/hx_geo);
                indy = floor(posy/hy_geo);
                indz = 0;
                    
                s_geo = ind(indx,indy,indz,nx_geo,ny_geo);

                ss = pindex(ii,jj,0,N);

                K->X[ss] = P_geo->X[s_geo]; 
                K->Y[ss] = P_geo->Y[s_geo]; 
                K->Z[ss] = P_geo->Z[s_geo]; 

               	}else{*/

                /* Global numeration of cells. */
                posx = (ii + 0.5)*h.x;
                posy = (jj + 0.5)*h.y;
                posz = (kk + 0.5)*h.z;

                indx = floor(posx/hx_geo);
                indy = floor(posy/hy_geo);
                indz = floor(posz/hz_geo);
                    
                s_geo = ind(indx,indy,indz,nx_geo,ny_geo);

                /* Local numeration of cells in subdomain s. */
                ss = pindex(ii,jj,kk,N);

                K->X[ss] = P_geo->X[s_geo]; 
                K->Y[ss] = P_geo->Y[s_geo]; 
                K->Z[ss] = P_geo->Z[s_geo];
		//}

            }
        }
    }

    if(P_geo !=NULL){
        free(P_geo->X);
        free(P_geo->Y);
        free(P_geo->Z);
        free(P_geo);
    }


    return 0;
}




void gaussian_permfield(Perm *K, Dim N, grid h, grid I)
{
/*
Read the absolute permeability from a file and set the permeability;
*/
    int err_chk, chk_scan; 
    double bb;
    FILE *fperm;
 
    int i, j, k, s, col; 
    int indx, indy, indz, s_geo; 
    int nx, ny, nz;
    double hx_geo, hy_geo, hz_geo;
    double posx, posy, posz;
    int ii,jj,kk,ss; 
    int iii,jjj,kkk,sss;
    int Ngeox, Ngeoy, Ngeoz, nt_geo;

    Ngeox = nxg; 
    Ngeoy = nyg;
    Ngeoz = nzg;
    
    nt_geo = Ngeox*Ngeoy*Ngeoz;


    Perm *P_geo = (Perm*)malloc(sizeof(Perm));
    P_geo->X = (double*)malloc((nt_geo)*sizeof(double));
    P_geo->Y = (double*)malloc((nt_geo)*sizeof(double));
    P_geo->Z = (double*)malloc((nt_geo)*sizeof(double));

    /* Set geological parameters */
    hx_geo = I.x/(double)(Ngeox);
    hy_geo = I.y/(double)(Ngeoy);
    hz_geo = I.z/(double)(Ngeox);

    /* Read permeability file */
    
    fperm = fopen(fn,"r");

    if (fperm == NULL)
        exit(EXIT_FAILURE);
    
    for(i = 0;  i < Ngeox; i++){
        for(j = 0; j < Ngeoy; j++){
            for(k = 0; k < Ngeoz; k++){	

                s = ind(i,j,k,Ngeox,Ngeoy);

                for(col = 0; col < 3; col++){	
                
                    chk_scan = fscanf(fperm,"%lf",&bb);

                    if (chk_scan != 1)
                        exit(EXIT_FAILURE);

                    if (col == 0){ 
                            P_geo->X[s] = bb; 
                            P_geo->Y[s] = bb; 
                            P_geo->Z[s] = bb; 
                    }
                } 
            }
        }
    }

    fclose(fperm);

    for (ii=0;ii<N.x;ii++){
        for (jj=0;jj<N.y;jj++){
            for (kk=0;kk<N.z;kk++){

                /* Global numeration of cells. */
                posx = (ii + 0.5)*h.x;
                posy = (jj + 0.5)*h.y;
                posz = (kk + 0.5)*h.z;

                indx = floor(posx/hx_geo);
                indy = floor(posy/hy_geo);
                indz = floor(posz/hz_geo);
                    
                s_geo = ind(indx,indy,indz,Ngeox,Ngeoy);

                /* Local numeration of cells in subdomain s. */
                ss = pindex(ii,jj,kk,N);

                K->X[ss] = P_geo->X[s_geo]; 
                K->Y[ss] = P_geo->Y[s_geo]; 
                K->Z[ss] = P_geo->Z[s_geo]; 

            }
        }
    }
   
    if(P_geo !=NULL){
        free(P_geo->X);
        free(P_geo->Y);
        free(P_geo->Z);
        free(P_geo);
    }

    return ;
}

void perm_homogeneous(Perm *K, Dim N)
{
/*
Define the constant permeability
*/
    int s, i, j, k; 

    for(s=0; s<N.t; s++){
        K->X[s] = 1.0e-13;
        K->Y[s] = 1.0e-13;
        K->Z[s] = 1.0e-13;
    }

}


// CASE 1D


int permeabilidade_karst(Perm *K, ELEK *keff,IND *indcol1, double *hs, DimK *M, FACEK *BoundKarst, double viscosidade, double **KI_c,int j)
{
int i, s, saux;
//internal elements
for (i=0; i<(M[j].s-1);i++){
	s=i;
	saux=i+1;
	keff->r[s][j]=hm(K->S[s][j],K->S[saux][j]);
	keff->r[s][j]=keff->r[s][j]/(hs[j]*hs[j]*viscosidade);
    indcol1[s].r=s+1;
  
}

for(i=1; i<M[j].s;i++){
	s=i;
	saux=i-1;
	keff->l[i][j]=hm(K->S[i][j],K->S[i-1][j]);
	keff->l[s][j]=keff->l[s][j]/(hs[j]*hs[j]*viscosidade);
    indcol1[s].l=s-1;
        
}


//Boundary elements

if(BoundKarst->l[j]==0)
{
	keff->l[0][j]=2*(K->S[0][j])/(hs[j]*hs[j]*viscosidade); 
    indcol1[0].l=-1;
}

if(BoundKarst->l[j]==1)
{
keff->l[0][j]=0;
indcol1[0].l=-1;
}


if(BoundKarst->l[j]==3) //Robin não linear
{
keff->l[0][j]=2*K->S[0][j]/(viscosidade*hs[j]*hs[j]);
indcol1[0].l=-1;
}


if(BoundKarst->r[j]==0)
{
	keff->r[M[j].s-1][j]=2*(K->S[M[j].s-1][j])/(hs[j]*hs[j]*viscosidade);
    indcol1[M[j].s-1].r=-1;
       
}
if(BoundKarst->r[j]==1)
{
keff->r[M[j].s-1][j]=0;
indcol1[M[j].s-1].r=-1;
}

if(BoundKarst->r[j]==3)//Robin não linear
{
keff->r[M[j].s-1][j]=2*K->S[M[j].s-1][j]/(viscosidade*hs[j]*hs[j]);
indcol1[M[j].s-1].r=-1;
}

// Case CENTER(C): sum of the Keff in each face//
for(i=0;i<M[j].s;i++){
s=i;
keff->c[s][j]=keff->l[s][j]+keff->r[s][j]+KI_c[s][j]/viscosidade;
indcol1[s].c=s;
}

return 0;

}




