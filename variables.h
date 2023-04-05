#ifndef VARIABLE_H 
#define VARIABLE_H

#include "structs.h"

/* Structs */ 
extern Dim N;
extern grid I;
extern grid h;
extern Perm *P;
extern ELE *perm;
extern IND *indcol;
extern ELE *v;
extern Bound *B;
extern FACE boundary;
extern FACEK *BoundKarst;
extern IND *indcol1;
extern ELEK *uc;
extern Perm *KI;
extern DimK *M;


extern double *wi;
extern int *IDwi;
extern int Nwi;
extern int ix0, ixf, jy0, jyf;
extern double karst;
extern double viscosidade;
extern double **beta;
extern double rho;
extern int *n_brechas;
extern int **brechas;
extern double **pc_old;
extern double **ucL_old;
extern double **ucR_old;
extern double *hs;
extern double *Kc;
extern double **cont_karst;
extern int number_karst;
extern int Mmax;
extern double **pc;
extern double **C;
extern double **C_contorno;

extern double nxg, nyg, nzg;
/* vectors */


/* vectors for the linear system matrix construction and solution by AMG */


/* Input files */
extern char fnam[1024];

/* fn: name of file to read the absolute permeability 
*/

extern const char *fn;  
extern const char *fnwi;  

extern const char *dirname;  

extern double time_total;
extern double time_amg;
extern int iter_amg;

extern int layer;
extern int lay;
 
#endif
