#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "structs.h"

// Create arrays:

double** create_matrix(int n, int m);
double* create_vector(int n);
int* create_int(int n);
int** create_matrix_int(int n, int m);

// Free arrays:

void libera(double **mat, int N, int M);

// Enumeration of elements:

int pindex(int i, int j, int k, Dim N);
int pxyindex(int i, int j, Dim N);
int pxzindex(int i, int k, Dim N);
int pyzindex(int j, int k, Dim N);
int ind(int i, int j, int k, int N, int M);

// Define permeabilities and effective permeabilities:
void set_perm_faces(Perm *K, grid I, grid h, Dim N, int Mmax, double *Kc);
int permeabilidade_arestas(Perm *K, ELE *perm, IND *indcol, grid h, Dim N, FACE boundary, double *WI, double *KI, double viscosidade);
int read_set_perm_field(Perm *K, Dim N, grid h, grid I);
void perm_homogeneous(Perm *K, Dim N);
double hm(double var1, double var2);
void gaussian_permfield(Perm *K, Dim N, grid h, grid I);
int read_layer_perm_field(Perm *K, Dim N, grid h, grid I);
int permeabilidade_karst(Perm *K, ELEK *keff,IND *indcol1, double *hs, DimK *M, FACEK *BoundKarst, double viscosidade, double **KI_c,int j);

// Define Boundary:

void create_boundary(double val_bound[7],Bound *B, Dim N,FACE boundary);
void set_boundary(double value, Bound *B, Dim N, int face);
void set_source_WI(Dim N);


// Define Source Term:
void source_term(Bound *B, ELE *perm, grid h, Dim N, FACE boundary, double *pc_old, double *pw, double *KI, double *WI, double viscosidade,double *b);
void source_term1d(ELEK *keff, double *hs, DimK *M, FACEK *BoundKarst,double **cont_karst, double **KI_c, double **p,Perm *K,ELE *perm, double *p_old,int j, double *f);    
// Reconstruct velocity:

int fluxos(double *p, ELE *perm, ELE *v, Bound *B, grid h, Dim N);
void divvel_pontual(ELE *v, grid h, Dim N);
void divvel_pontual_1d(ELE *uc, grid h, Dim N);
int fluxo1d(double **x, ELEK *keff, ELEK *uc, double *hs, DimK *M, double **cont_karst, FACEK *BoundKarst, Perm *K, ELE *perm, double *p_old);
// Solve linear system by AMG:

int compute(int *hdr, int *cnt, int *col, double *ele, double *b, double *x);
void solve_linsys(ELE *perm, IND *indcol, double *b, double *x, Dim N);
void prod_mat_vec(ELE *perm, IND *indcol, double *x, double *res, Dim N);
void prod_mat_vec1d(ELE *keff, IND *indcol1, double *x, double *res, DimK *M,int j);
void solve_linsys_1d(ELEK *keff, IND *indcol1, double *f, double *x1, DimK *M,int j);

// Output functions:

void save_flux_Right(ELE *v, Dim N);
void save_flux_Left(ELE *v, Dim N);
void save_flux_Up(ELE *v, Dim N);
void save_flux_Down(ELE *v, Dim N);
void save_flux_Top(ELE *v, Dim N);
void save_flux_Bottom(ELE *v, Dim N);
void save_pressure(double *x, Dim N, grid h);
void save_flux(ELE *v, Dim N, grid h);
void erros(ELE *perm, IND *indcol, double *x, double *b, double *r, double *res, grid h, Dim N);

void erros_1d(ELE *keff_uni, IND *indcol1, double *x, double *b, double *r, double *res, grid h, Dim N);
void save_flux_1d(ELEK *uc, DimK *M, double *hs, grid h);
void save_pressure_1d(double **x, DimK *M, double *hs, grid h);
void save_flux_Left_1d(ELE *uc, Dim N);
void save_flux_Right_1d(ELE *uc, Dim N);

void reservoir(Dim N, Perm *K,double *WI, double *KI, double **pc_old, IND *indcol, Bound *B, ELE *perm, grid h, FACE boundary, double viscosidade, double *pw,DimK *M, double *p);
void conduit(double **p_aux,double **cont_karst, double **KI_c, double viscosidade, DimK *M, double *hs,IND *indcol1, FACEK *BoundKarst, ELEK *keff, Perm *K,ELE *perm, double *p_old,double **pc);

// VTK (.vts) print
int vts_StructuredGrid_pressure_file(double *x, Dim N, grid h); 
int vts_StructuredGrid_flux_file(ELE *v, Dim N, grid h);
int vts_StructuredGrid_file(ELE *v, double *x, Dim N, grid h);
#endif
