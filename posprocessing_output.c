#include "structs.h"
#include "variables.h"
#include "functions.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void save_flux_Right(ELE *v, Dim N)
{
/* save right flux in a file called: R_flux.dat.
*/   

    int i,j,k,s;
    FILE *fr;

   /* Right Flux */ 
    sprintf(fnam,"R_flux.dat");
    fr = fopen(fnam, "w");
    if(fr == NULL){
        printf("Unable to open file %s. Terminating run... \n", fnam);
        exit(0);
    }    

    for(i=0; i < N.x; i++){
        for(j=0; j < N.y; j++){
            for(k=0; k< N.z; k++){
                
                s = pindex(i,j,k,N);
                fprintf(fr, "%g ", v[s].r);

            }
        fprintf(fr, " \n");
        }
	fprintf(fr, " \n");
    }

    fclose(fr);
}

void save_flux_Left(ELE *v, Dim N)
{
/* save left flux in a file called: L_flux.dat.
*/   

    int i,j,k,s;
    FILE *fl;

   /* Left Flux */ 
    sprintf(fnam,"L_flux.dat");
    fl = fopen(fnam, "w");
    if(fl == NULL){
        printf("Unable to open file %s. Terminating run... \n", fnam);
        exit(0);
    }    

    for(i=0; i < N.x; i++){
        for(j=0; j < N.y; j++){
            for(k=0; k< N.z; k++){
                
                s = pindex(i,j,k,N);
                fprintf(fl, "%g ", v[s].l);

            }
        fprintf(fl, " \n");
        }
	fprintf(fl, " \n");
    }

    fclose(fl);
}

void save_flux_Down(ELE *v, Dim N)
{
/* save down flux in a file called: D_flux.dat.
*/   

    int i,j,k,s;
    FILE *fd;

   /* Down Flux */ 
    sprintf(fnam,"D_flux.dat");
    fd = fopen(fnam, "w");
    if(fd == NULL){
        printf("Unable to open file %s. Terminating run... \n", fnam);
        exit(0);
    }    

    for(i=0; i < N.x; i++){
        for(j=0; j < N.y; j++){
            for(k=0; k< N.z; k++){
                
                s = pindex(i,j,k,N);
                fprintf(fd, "%g ", v[s].d);

            }
        fprintf(fd, " \n");
        }
	fprintf(fd, " \n");
    }

    fclose(fd);
}

void save_flux_Up(ELE *v, Dim N)
{
/* save up flux in a file called: U_flux.dat.
*/   

    int i,j,k,s;
    FILE *fu;

   /* Up Flux */ 
    sprintf(fnam,"U_flux.dat");
    fu = fopen(fnam, "w");
    if(fu == NULL){
        printf("Unable to open file %s. Terminating run... \n", fnam);
        exit(0);
    }    

    for(i=0; i < N.x; i++){
        for(j=0; j < N.y; j++){
            for(k=0; k< N.z; k++){
                
                s = pindex(i,j,k,N);
                fprintf(fu, "%g ", v[s].u);

            }
        fprintf(fu, " \n");
        }
	fprintf(fu, " \n");
    }

    fclose(fu);
}

void save_flux_Bottom(ELE *v, Dim N)
{
/* save bottom flux in a file called: B_flux.dat.
*/   

    int i,j,k,s;
    FILE *fb;

   /* Bottom Flux */ 
    sprintf(fnam,"B_flux.dat");
    fb = fopen(fnam, "w");
    if(fb == NULL){
        printf("Unable to open file %s. Terminating run... \n", fnam);
        exit(0);
    }    

    for(i=0; i < N.x; i++){
        for(j=0; j < N.y; j++){
            for(k=0; k< N.z; k++){
                
                s = pindex(i,j,k,N);
                fprintf(fb, "%g ", v[s].b);

            }
        fprintf(fb, " \n");
        }
	fprintf(fb, " \n");
    }

    fclose(fb);
}

void save_flux_Top(ELE *v, Dim N)
{
/* save top flux in a file called: T_flux.dat.
*/   

    int i,j,k,s;
    FILE *ft;

   /* Top Flux */ 
    sprintf(fnam,"T_flux.dat");
    ft = fopen(fnam, "w");
    if(ft == NULL){
        printf("Unable to open file %s. Terminating run... \n", fnam);
        exit(0);
    }    

    for(i=0; i < N.x; i++){
        for(j=0; j < N.y; j++){
            for(k=0; k< N.z; k++){
                
                s = pindex(i,j,k,N);
                fprintf(ft, "%g ", v[s].t);

            }
        fprintf(ft, " \n");
        }
	fprintf(ft, " \n");
    }

    fclose(ft);
}

void save_pressure(double *x, Dim N, grid h)
{
/* save pressure in a file called: pressure.dat.
    write in 4 columns:

    x(i) | y(j) | z(k) | pressure(i,j,k)

*/   

    int i,j,k,s;
    FILE *saida_pressao;

   /* Pressure */ 
    sprintf(fnam,"%s/pressure_reservoir.dat",dirname);
    saida_pressao = fopen(fnam, "w");
    if(saida_pressao == NULL){
        printf("Unable to open file %s. Terminating run... \n", fnam);
        exit(0);
    }    

    for(i=0; i < N.x; i++){
        for(j=0; j < N.y; j++){
            for(k=0; k< N.z; k++){
                
                s = pindex(i,j,k,N);
                fprintf(saida_pressao, "%.6lf %.6lf %.6lf %1.9e\n", h.x*(i + 0.5),h.y*(j + 0.5),h.z*(k + 0.5), x[s]);

            }
        fprintf(saida_pressao, " \n");
        }
	fprintf(saida_pressao, " \n");
    }

    fclose(saida_pressao);
}

void save_flux(ELE *v, Dim N, grid h)
{
/* save flux in a file called: flux.dat.
    write in 6 columns and for the flux at the center of elements:

    x(i) | y(j) | z(k) | vx(i,j,k) | vy(i,j,k) | vz(i,j,k)

    and

    vx = vR - vL;
    vy = vU - vD;
    vz = vT - vB;

    we also normalize the vector by norm 2.
*/

    int i,j,k,s;
    double vx, vy, vz, v_norm;
    FILE *saida_fluxo;

    /* Flux */ 
    sprintf(fnam,"%s/flux_reservoir.dat",dirname);
    saida_fluxo = fopen(fnam, "w");
    if(saida_fluxo == NULL){
        printf("Unable to open file %s. Terminating run... \n", fnam);
        exit(0);
    }
    for(i=0; i< N.x; i++){
        for(j=0; j< N.y; j++){
            for(k=0; k< N.z; k++){
            
                s = pindex(i,j,k,N);

                vx = 0.5*(v[s].r - v[s].l);
                vy = 0.5*(v[s].u - v[s].d);
                vz = 0.5*(v[s].t - v[s].b);

                //v_norm = sqrt(vx*vx + vy*vy + vz*vz);
                v_norm = 1.0;

                if (v_norm > 0){
                    /*
                    vx = 1e-1*vx/v_norm;
                    vy = 1e-1*vy/v_norm;
                    vz = 1e-1*vz/v_norm;
                    */

                }
                else{
                    vx = 0.0;
                    vy = 0.0;
                    vz = 0.0;
                }

                //fprintf(saida_fluxo, "%.6lf %.6lf %.6lf %1.9e %1.9e %1.9e %1.9e %1.9e %1.9e %1.9e %1.9e %1.9e\n", h.x*(i + 0.5),h.y*(j + 0.5),h.z*(k + 0.5),vx,vy,vz, v[s].r,v[s].l,v[s].u,v[s].d,v[s].t,v[s].b);
                fprintf(saida_fluxo, "%.6lf %.6lf %.6lf %1.9e %1.9e %1.9e \n", h.x*(i + 0.5),h.y*(j + 0.5),h.z*(k + 0.5),vx,vy,vz);
            }
        }
    }
    fclose(saida_fluxo);
}    

void erros(ELE *perm, IND *indcol, double *x, double *b, double *r, double *res, grid h, Dim N)
{
/* Calculation and file saving of Residues and errors */ 

    int i,j,k,s;
    /* Norm/resuduals variables.
        rs1: vector norm 1;
        rs2: vector norm 2;
        rsmax: vector norm max;
        add: auxiliar variable;
        rs2_norm:  vector norm 2 divided by number of elements;
    */
    double rs1, rs2, rsmax, add, rs2_norm;
    FILE *tabela;

    /* initialize residual variables */
    rs1 = 0.0; rs2 = 0.0; rsmax = 0.0; add = 0.0; rs2_norm = 0.0;

    /* Residues on finest grid */ 
    prod_mat_vec(perm,indcol,x,res,N);
    
    /*
    for(i=0; i < N.x; i++){
        for(j=0; j< N.y; j++){
            for(k=0; k< N.z; k++){
                s = pindex(i,j,k,N);
                printf("x[%d] = %e, res[%d] = %e\n",s,x[s],s,res[s]);
            }
        }
    }
    printf("============================\n");


    for(i=0; i < N.x; i++){
        for(j=0; j< N.y; j++){
            for(k=0; k< N.z; k++){
                s = pindex(i,j,k,N);
                r[s] = b[s] - res[s];
                printf("b[%d] = %e, res[%d] = %e\n",s,b[s],s,res[s]);
            }
        }
    }
*/
    /* Calculation */ 
    for(i=0; i < N.x; i++){
        for(j=0; j< N.y; j++){
            for(k=0; k< N.z; k++){
                s = pindex(i,j,k,N);
                /* Norm-1/-2 */ 
                rs2 = rs2 + h.x*h.y*h.z*(r[s]*r[s]);
                add = fabs(r[s]);
                rs1 = rs1 + add;

                /* Norm-max */ 
                if(add > rsmax){
                    rsmax = add;
                }
                add = 0.0;
            }
        }
    }
    // Norm-2:
    rs2 = sqrt(rs2);

    // ||r||/sqrt(N.t):
    add = sqrt(N.t);
    rs2_norm = rs2/add;

    sprintf(fnam,"errors_table.dat");
    tabela = fopen(fnam, "w");
    if(tabela == NULL){
        printf("Unable to open file %s. Terminating run... \n", fnam);
        exit(0);
    }
    fprintf(tabela,"res 1 & res 2 & res max & res norma \n");
    fprintf(tabela, "%g & %g & %g & %g  \n \n \n", rs1, rs2, rsmax, rs2_norm);
    fprintf(tabela,"Dim x & Dim y & Dim z \n");
    fprintf(tabela, "%d & %d & %d \n", N.x, N.y, N.z);


    fclose(tabela);
}


//1d case

void save_flux_Right_1d(ELE *uc, Dim N)
{
/* save right flux in a file called: R_flux.dat.
*/   

    int i,s;
    FILE *fr;

   /* Right Flux */ 
    sprintf(fnam,"R_flux.dat");
    fr = fopen(fnam, "w");
    if(fr == NULL){
        printf("Unable to open file %s. Terminating run... \n", fnam);
        exit(0);
    }    

    for(i=0; i < N.t; i++){
                       
                s = i;
                fprintf(fr, "%g ", uc[s].r);

            }
 
    fclose(fr);
}

void save_flux_Left_1d(ELE *uc, Dim N)
{
/* save left flux in a file called: L_flux.dat.
*/   

    int i,s;
    FILE *fl;

   /* Left Flux */ 
    sprintf(fnam,"L_flux.dat");
    fl = fopen(fnam, "w");
    if(fl == NULL){
        printf("Unable to open file %s. Terminating run... \n", fnam);
        exit(0);
    }    

    for(i=0; i < N.t; i++){
                      
                s = i;
                fprintf(fl, "%g ", uc[s].l);

            }

    fclose(fl);
}


void save_pressure_1d(double **x, DimK *M, double *hs,grid h)
{
    int i,s,j, dirx,diry,dirz,saux;
    FILE *saida_pressao;

   /* Pressure */ 
   

    for(j=0;j<number_karst;j++){
         sprintf(fnam,"%s/pressure_karst_%i.dat",dirname,j);
    saida_pressao = fopen(fnam, "w");
    if(saida_pressao == NULL){
        printf("Unable to open file %s. Terminating run... \n", fnam);
        exit(0);
    } 

    if(M[j].pos==1){
        dirx=M[j].sxi;
         for(i=0;i<M[j].s;i++){
                    s = dirx;
                    fprintf(saida_pressao, "%.6lf %.6lf %.6lf %1.9e\n", h.x*(s + 0.5),h.y*(M[j].syi + 0.5),h.z*(M[j].szi + 0.5), x[i][j]);
                    dirx++;
                    }  
    } 

    else if(M[j].pos==2){
        diry=M[j].syi;
         for(i=0;i<M[j].s;i++){
                    s = diry;
                    fprintf(saida_pressao, "%.6lf %.6lf %.6lf %1.9e\n", h.x*(M[j].sxi+ 0.5),h.y*(s + 0.5),h.z*(M[j].szi + 0.5), x[i][j]);
                    diry++;
                    }  
    } 

    else if(M[j].pos==3){
        dirz=M[j].szi;
         for(i=0;i<M[j].s;i++){
                    s = dirz;
                    fprintf(saida_pressao, "%.6lf %.6lf %.6lf %1.9e\n", h.x*(M[j].sxi+ 0.5),h.y*(M[j].syi + 0.5),h.z*(s + 0.5), x[i][j]);
                    dirz++;
                    }  
    }


    else if(M[j].pos==4){
        dirx=M[j].sxi;
        diry=M[j].syi;
         for(i=0;i<M[j].s;i++){
                    s = dirx;
                    saux = diry;
                    fprintf(saida_pressao, "%.6lf %.6lf %.6lf %1.9e\n", h.x*(s+ 0.5),h.y*(saux+ 0.5),h.z*(M[j].szi + 0.5), x[i][j]);
                    dirx++;
                    diry++;
                    }  
    }

     else if(M[j].pos==5){
        dirx=M[j].sxi;
        diry=M[j].syi;
         for(i=0;i<M[j].s;i++){
                    s = dirx;
                    saux = diry;
                    fprintf(saida_pressao, "%.6lf %.6lf %.6lf %1.9e\n", h.x*(s+ 0.5),h.y*(saux+ 0.5),h.z*(M[j].szi + 0.5), x[i][j]);
                    dirx++;
                    diry--;
                    }  
    }

    else if(M[j].pos==6){
        dirx=M[j].sxi;
        diry=M[j].syi;
         for(i=0;i<M[j].s;i++){
                    s = dirx;
                    saux = diry;
                    fprintf(saida_pressao, "%.6lf %.6lf %.6lf %1.9e\n", h.x*(s+ 0.5),h.y*(saux+ 0.5),h.z*(M[j].szi + 0.5), x[i][j]);
                    dirx--;
                    diry++;
                    }  
    }

    else if(M[j].pos==7){
        dirx=M[j].sxi;
        diry=M[j].syi;
         for(i=0;i<M[j].s;i++){
                    s = dirx;
                    saux = diry;
                    fprintf(saida_pressao, "%.6lf %.6lf %.6lf %1.9e\n", h.x*(s+ 0.5),h.y*(saux+ 0.5),h.z*(M[j].szi + 0.5), x[i][j]);
                    dirx--;
                    diry--;
                    }  
    }
       
 
}


fclose(saida_pressao);
}

void save_flux_1d(ELEK *uc, DimK *M, double *hs,grid h)
{

    int i,s,dirx,diry,dirz,j,saux;
    double vx, v_norm;
    FILE *saida_fluxo;

    /* Flux */ 
    
    for(j=0;j<number_karst;j++){
        sprintf(fnam,"%s/flux_karst_%i.dat",dirname,j);
    saida_fluxo = fopen(fnam, "w");
    if(saida_fluxo == NULL){
        printf("Unable to open file %s. Terminating run... \n", fnam);
        exit(0);
    }


    // for(i=0;i<M[j].s;i++){
    //        vx = 0.5*(uc->r[i][j] - uc->l[i][j]);
    //         v_norm = 1.0;

    //             if (v_norm > 0){
                    
    //                 vx = 1e-1*vx/v_norm;
                   
                    

    //             }
    //             else{
    //                 vx = 0.0;
                 
    //             }
    // }

        if(M[j].pos==1){
        dirx=M[j].sxi;
         for(i=0;i<M[j].s;i++){
                    s = dirx;
                    vx = 0.5*(uc->r[i][j] - uc->l[i][j]);
                
                    fprintf(saida_fluxo, "%.6lf %.6lf %.6lf %1.9e\n", h.x*(s + 0.5),h.y*(M[j].syi + 0.5),h.z*(M[j].szi + 0.5), vx);
                    dirx++;
                    }  
    } 

    else if(M[j].pos==2){
        diry=M[j].syi;
         for(i=0;i<M[j].s;i++){
                    s = diry;
                    vx = 0.5*(uc->r[i][j] - uc->l[i][j]);
                    fprintf(saida_fluxo, "%.6lf %.6lf %.6lf %1.9e\n", h.x*(M[j].sxi+ 0.5),h.y*(s + 0.5),h.z*(M[j].szi + 0.5), vx);
                    diry++;
                    }  
    } 

    else if(M[j].pos==3){
        dirz=M[j].szi;
         for(i=0;i<M[j].s;i++){
                    s = dirz;
                    vx = 0.5*(uc->r[i][j] - uc->l[i][j]);
                    fprintf(saida_fluxo, "%.6lf %.6lf %.6lf %1.9e\n", h.x*(M[j].sxi+ 0.5),h.y*(M[j].syi + 0.5),h.z*(s + 0.5), vx);
                    dirz++;
                    }  
    }


    else if(M[j].pos==4){
        dirx=M[j].sxi;
        diry=M[j].syi;
         for(i=0;i<M[j].s;i++){
                    s = dirx;
                    saux = diry;
                    vx = 0.5*(uc->r[i][j] - uc->l[i][j]);
                    fprintf(saida_fluxo, "%.6lf %.6lf %.6lf %1.9e\n", h.x*(s+ 0.5),h.y*(saux+ 0.5),h.z*(M[j].szi + 0.5), vx);
                    dirx++;
                    diry++;
                    }  
    }

     else if(M[j].pos==5){
        dirx=M[j].sxi;
        diry=M[j].syi;
         for(i=0;i<M[j].s;i++){
                    s = dirx;
                    saux = diry;
                    vx = 0.5*(uc->r[i][j] - uc->l[i][j]);
                    fprintf(saida_fluxo, "%.6lf %.6lf %.6lf %1.9e\n", h.x*(s+ 0.5),h.y*(saux+ 0.5),h.z*(M[j].szi + 0.5), vx);
                    dirx++;
                    diry--;
                    }  
    }

    else if(M[j].pos==6){
        dirx=M[j].sxi;
        diry=M[j].syi;
         for(i=0;i<M[j].s;i++){
                    s = dirx;
                    saux = diry;
                    vx = 0.5*(uc->r[i][j] - uc->l[i][j]);
                    fprintf(saida_fluxo, "%.6lf %.6lf %.6lf %1.9e\n", h.x*(s+ 0.5),h.y*(saux+ 0.5),h.z*(M[j].szi + 0.5), vx);
                    dirx--;
                    diry++;
                    }  
    }

    else if(M[j].pos==7){
        dirx=M[j].sxi;
        diry=M[j].syi;
         for(i=0;i<M[j].s;i++){
                    s = dirx;
                    saux = diry;
                    vx = 0.5*(uc->r[i][j] - uc->l[i][j]);
                    fprintf(saida_fluxo, "%.6lf %.6lf %.6lf %1.9e\n", h.x*(s+ 0.5),h.y*(saux+ 0.5),h.z*(M[j].szi + 0.5), vx);
                    dirx--;
                    diry--;
                    }  
    }

 // fprintf(saida_fluxo, "%.6lf  %1.9e  \n", hs[j]*(s + 0.5),vx);
 //                dir++;

    }
   
    fclose(saida_fluxo);
}    

void erros_1d(ELE *keff_uni, IND *indcol1, double *x, double *b, double *r, double *res, grid h, Dim N)
{
/* Calculation and file saving of Residues and errors */ 

    int i,s;
    /* Norm/resuduals variables.
        rs1: vector norm 1;
        rs2: vector norm 2;
        rsmax: vector norm max;
        add: auxiliar variable;
        rs2_norm:  vector norm 2 divided by number of elements;
    */
    double rs1, rs2, rsmax, add, rs2_norm;
    FILE *tabela;

    /* initialize residual variables */
    rs1 = 0.0; rs2 = 0.0; rsmax = 0.0; add = 0.0; rs2_norm = 0.0;

    /* Residues on finest grid */ 
    //prod_mat_vec1d(keff_uni,indcol1,x,res,N);
    
    /*
    for(i=0; i < N.t; i++){
       
                s = i;
                printf("x[%d] = %e, res[%d] = %e\n",s,x[s],s,res[s]);
     
    }
    printf("============================\n");


    for(i=0; i < N.t; i++){
       
                s = i;
                r[s] = b[s] - res[s];
                printf("b[%d] = %e, res[%d] = %e\n",s,b[s],s,res[s]);
      
    }
*/
    /* Calculation */ 
    for(i=0; i < N.t; i++){
               s = i;
                /* Norm-1/-2 */ 
                rs2 = rs2 + h.x*(r[s]*r[s]);
                add = fabs(r[s]);
                rs1 = rs1 + add;

                /* Norm-max */ 
                if(add > rsmax){
                    rsmax = add;
                }
                add = 0.0;
      
    }
    // Norm-2:
    rs2 = sqrt(rs2);

    // ||r||/sqrt(N.t):
    add = sqrt(N.t);
    rs2_norm = rs2/add;

    sprintf(fnam,"errors_table.dat");
    tabela = fopen(fnam, "w");
    if(tabela == NULL){
        printf("Unable to open file %s. Terminating run... \n", fnam);
        exit(0);
    }
    fprintf(tabela,"res 1 & res 2 & res max & res norma \n");
    fprintf(tabela, "%g & %g & %g & %g  \n \n \n", rs1, rs2, rsmax, rs2_norm);
    fprintf(tabela,"Dim x & Dim y & Dim z \n");
    fprintf(tabela, "%d & %d & %d \n", N.x, N.y, N.z);


    fclose(tabela);
}
