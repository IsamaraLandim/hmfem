#include "structs.h"
#include "variables.h"
#include "functions.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>

/* Global variables declaration */
/* Structs */ 
Dim N;
grid I;
grid h;
Perm *K;
ELE *perm; //keff
IND *indcol;
IND *indcol1;
ELE *v;
ELEK *uc;
Bound *B;
FACE boundary;
FACEK *BoundKarst;
DimK *M;

double nxg, nyg, nzg;
double viscosidade;
double val_bound[7];
double *hs;
double *Kc;
int number_karst;
double **cont_karst;
double **beta;
double rho;
double **pc_old, **ucL_old, **ucR_old;
int **brechas;
int *n_brechas; 
double **pc;
double **C;
double **C_contorno;
int Mmax;
double *wi;
int *IDwi;
int Nwi;
int ix0, ixf, jy0, jyf;
/* Input files */
char fnam[1024];
char fname[1024];
FILE *fn_time;

const char *fn;
const char *fnwi;
const char *dirname;
int layer;
double karst;
double time_total;
double time_amg;
int iter_amg;
double timeW_amg;
double timeW_amg2;
double timeWtotal; 
 
double get_wall_time(){
    struct timeval time;

    if (gettimeofday(&time,NULL)) //  Handle error
        return 0;
                    
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}


double *atuaVector(double *v1, double *v2, int n)
{
int i;
for(i=0;i<n;i++){
v1[i]=v2[i];}
return(v1);
}

double **atuaMatriz(double **v1,double **v2,int m,int n)
{
    int i,j;
    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++){
            v1[i][j]=v2[i][j];
        }
    }
    return(v1);
}

//maximum number of elements in karst conduits
double maximo(DimK *M, int n)
{
int i;
double max;
max=-1;
for(i=0;i<n;i++){
 if(M[i].s>=max){
max=M[i].s; }
}
return max;
}

//relative error
double erro_relativo(double  *y,double *y_old, int n)
{
	int i;
	double sum1, sum2,result;
	sum1=0;
	sum2=0;

	for(i=0;i<n;i++){
sum1=sum1+(y[i]-y_old[i])*(y[i]-y_old[i]);
sum2=sum2+y[i]*y[i];	}
result=sqrt(sum1)/sqrt(sum2);
	return(result);
}

/* Begin main */
int main(int argc, char* argv[])
{

double ti = get_wall_time();
double tolerancia=1e-08;
int iter_max=1000;
int type;
nxg = 60;
nyg = 60;
nzg = 60; 
viscosidade=1e-03;
rho=1000;
double Ppoco=10e06;
double well_index=7.64352339e-13;//6.761e-13;//


/* Set boundary values at each face.
    val_bound[1]: value of boundary (Dirichlet or Neumann) at Left boundary;
    val_bound[2]: value of boundary (Dirichlet or Neumann) at Right boundary;
    val_bound[3]: value of boundary (Dirichlet or Neumann) at Upper boundary;
    val_bound[4]: value of boundary (Dirichlet or Neumann) at Down boundary;
    val_bound[5]: value of boundary (Dirichlet or Neumann) at Top boundary;
    val_bound[6]: value of boundary (Dirichlet or Neumann) at Bottom boundary;
    */
    val_bound[1] = 60e06;
    val_bound[2] = 60e06;
    val_bound[3] = 60e06;
    val_bound[4] = 60e06;
    val_bound[5] = 0.0;
    val_bound[6] = 0.0;
    printf("L=%e, R=%e, U=%e, D=%e, T=%e, B=%e\n", val_bound[1], val_bound[2],val_bound[3],val_bound[4],val_bound[5],val_bound[6]);
   

    N.x = atoi(argv[1]);
    N.y = atoi(argv[2]);
    N.z = atoi(argv[3]);
    N.xy=N.x*N.y;
    N.t=N.x*N.y*N.z;
    N.px=N.x/2; // position of well in x
    N.py=N.y/2; // position of well in y
    N.pz=10; // position of well in z
    I.x = atof(argv[4]);
    I.y = atof(argv[5]);
    I.z = atof(argv[6]);
    h.x = I.x/N.x;
    h.y = I.y/N.y;
    h.z = I.z/N.z; 

    number_karst = atoi(argv[8]);
    type = atoi(argv[9]);
      
        
    printf("number_karst=%d, type=%d\n",number_karst,type);  
    printf("Nx = %d, Ny = %d, Nz = %d, Npz=%d , Npx=%d, Npy=%d \n",N.x,N.y,N.z,N.pz,N.px,N.py);
    printf("Ix = %f, Iy = %f, Iz = %f \n",I.x,I.y,I.z);
    printf("hx = %f, hy = %f, hz = %f \n",h.x,h.y,h.z);

    int i, j, k, s;

    FILE *arq;

DimK *M=(DimK*)malloc((number_karst+1)*sizeof(DimK));
Kc=create_vector(number_karst);
FACEK *BoundKarst=(FACEK*)malloc(sizeof(FACEK));
BoundKarst->l= create_int(number_karst);     
BoundKarst->r= create_int(number_karst);    
   

cont_karst=create_matrix(2,number_karst);
C_contorno=create_matrix(2,number_karst);
double karst_index[number_karst], Ac[number_karst];
n_brechas=create_int(number_karst);


char name[1024];

sprintf(name,"data_%d_%d.dat",number_karst,type);


//LEITURA DE DADOS DOS CONDUTOS
arq=fopen(name,"r");
if(arq==NULL)
{
	printf("ERROR! The file was not opened!\n");
	exit(1);
}
else
{	printf("The file was successfully opened!\n");
}
//POSICAO DO KARST
for(i=0;i<number_karst;i++)
{fscanf(arq,"%i", &M[i].s);}
for(i=0;i<number_karst;i++)
{fscanf(arq,"%i",&M[i].sxi);}
for(i=0;i<number_karst;i++)
{fscanf(arq,"%i",&M[i].syi);}
for(i=0;i<number_karst;i++)
{fscanf(arq,"%i",&M[i].szi);}
for(i=0;i<number_karst;i++)
{fscanf(arq,"%i",&M[i].pos);}

//PERMEABILIDADE DO KARST
for(i=0;i<number_karst;i++)
{fscanf(arq,"%lf",&Kc[i]);}

//TIPO DE CONTORNO 
for(i=0;i<number_karst;i++){
	fscanf(arq,"%i",&BoundKarst->l[i]);}
for(i=0;i<number_karst;i++){
	fscanf(arq,"%i",&BoundKarst->r[i]);}

// VALOR NO CONTORNO
for(j=0;j<2;j++){
	for(i=0;i<number_karst;i++){
fscanf(arq,"%lf",&cont_karst[j][i]);}}

// Constante C no contorno
for(j=0;j<2;j++){
    for(i=0;i<number_karst;i++){
fscanf(arq,"%lf",&C_contorno[j][i]);}}

//karst index and condut's area
for(i=0;i<number_karst;i++){
	fscanf(arq,"%lf",&karst_index[i]);}
for(i=0;i<number_karst;i++){
	fscanf(arq,"%lf",&Ac[i]);}
   
// QUANTIDADE DE BRECHAS
for(i=0;i<number_karst;i++){
	fscanf(arq,"%i",&n_brechas[i]);}

//NUMERO MAXIMO DE BRECHAS
int max_brechas=-1;
for(i=0;i<number_karst;i++){
	if(n_brechas[i]>=max_brechas){
		max_brechas=n_brechas[i];}
}

// BETA
beta=create_matrix(max_brechas,number_karst);
for(j=0;j<max_brechas;j++){
	for(i=0;i<number_karst;i++){
		fscanf(arq,"%lf",&beta[j][i]);
		beta[j][i]=beta[j][i]-1;} //beta=R^4-1
}

// Constante C das brechas
C=create_matrix(max_brechas,number_karst);
for(j=0;j<max_brechas;j++){
    for(i=0;i<number_karst;i++){
        fscanf(arq,"%lf",&C[j][i]);}
}

// LOCALIZACAO DAS BRECHAS (indice j do conduto)
brechas=create_matrix_int(max_brechas,number_karst);
for(j=0;j<max_brechas;j++){
	for(i=0;i<number_karst;i++){
	fscanf(arq,"%d",&brechas[j][i]);}
}

fclose(arq);


    hs=create_vector(number_karst);
    for(i=0;i<number_karst;i++)
    {
	if(M[i].pos==1){
		hs[i]=h.x;
	}
	else if(M[i].pos==2){
		hs[i]=h.y;
	}
	else if(M[i].pos==3){
		hs[i]=h.z;
	}
    
    }


Mmax=maximo(M,number_karst);


double *p, *p_old, *pw;
   pc=create_matrix(Mmax,number_karst);
   pc_old=create_matrix(Mmax,number_karst);
   p=create_vector(N.t);
   p_old=create_vector(N.t);
   

double *IDwi;
pw=create_vector(N.t);
IDwi=create_vector(N.t);

double **KI_c,*KI;    
KI_c=create_matrix(Mmax,number_karst);
KI=create_vector(N.t);


   int aux;
   for(i=N.pz;i<N.z;i++){
   aux=pindex(N.px,N.py,i,N);
   pw[aux]=Ppoco;
   IDwi[aux]=1;
   }
 

    Perm *K = (Perm*)malloc( sizeof(Perm));
    K->X = create_vector(N.t); /* x direction */        
    K->Y = create_vector(N.t); /* y direction */
    K->Z = create_vector(N.t); /* y direction */
    K->S= create_matrix(Mmax,number_karst);      

    /* Set type of permeability
       P->permeability = 1: homogeneous permeability
       P->permeability = 2: read permeability of a file
       P->permeability = 3: You write your own set absolute permeability tensor 
   */
    K->permeability = atoi(argv[7]);
    //fn = argv[10];
    fnwi = argv[11];
    dirname = argv[10];
    

    /* Set type of boundary on each face (Righ, Left, Up, Down, Top, Bottom).
       boundary.x = 1: Neumann boundary;
       boundary.x = 0: Dirichlet boundary;
       boundary.x = 2: Robin Linear
       boundary.x = 3: Robin não linear
       x = r,l,u,d,t,b;
    */

    boundary.l = 0;
    boundary.r = 0;
    boundary.d = 0;
    boundary.u = 0;
    boundary.t = 1;
    boundary.b = 1;
 
ucL_old=create_matrix(max_brechas,number_karst);
ucR_old=create_matrix(max_brechas,number_karst);

ELE* perm = (ELE*)malloc((N.t+1)*sizeof(ELE)); // keff do reservatorio
ELEK *keff = (ELEK*)malloc((Mmax+1)*sizeof(ELEK));// keff do conduite
keff->l = create_matrix(Mmax,number_karst); /* left face */        
keff->r = create_matrix(Mmax,number_karst);  /* right face */
keff->c = create_matrix(Mmax,number_karst);  /* upper face */

uc = (ELEK*)malloc((Mmax+1)*sizeof(ELEK)); // fluxo no conduite
uc->l = create_matrix(Mmax,number_karst); /* left face */        
uc->r = create_matrix(Mmax,number_karst);  /* right face */
uc->c = create_matrix(Mmax,number_karst);  /* upper face */
v = (ELE*)malloc((N.t+1)*sizeof(ELE)); /* flux  no reservatorio*/
indcol1 = (IND*)malloc((Mmax+1)*sizeof(IND));
indcol = (IND*)malloc((N.t+1)*sizeof(IND)); /* auxiliar index to identify the coefficients of pressure matrix */ 

/* Contorno */
    Bound *B = (Bound*)malloc( sizeof(Bound));
    B->l = create_vector(N.t); /* left face */        
    B->r = create_vector(N.t); /* right face */
    B->u = create_vector(N.t); /* upper face */
    B->d = create_vector(N.t); /* down face */
    B->t = create_vector(N.t); /* top face */
    B->b = create_vector(N.t); /* bottom face */

   create_boundary(val_bound, B, N, boundary);


double *WI;
WI=create_vector(N.t);
set_perm_faces(K,I,h,N,Mmax,Kc);


 reservoir(N, K, WI, KI,pc_old,indcol,B,perm, h,boundary,viscosidade,pw,M,p);
 p_old=atuaVector(p_old,p,N.t);

int Nwi=0;
for(i=N.pz;i<N.z;i++){
   aux=pindex(N.px,N.py,i,N);
   WI[aux]=well_index;
    Nwi=Nwi+1; }

 int posx,posy,posz;   
 double **p_aux;
p_aux=create_matrix(Mmax,number_karst);

for(j=0;j<number_karst;j++){
	if(M[j].pos==1){
		posx=M[j].sxi;
		for(i=0;i<M[j].s;i++){
			aux=pindex(posx,M[j].syi,M[j].szi,N);
			KI[aux]=karst_index[j];
			KI_c[i][j]=karst_index[j]/(Ac[j]*hs[j]);
			p_aux[i][j]=p_old[aux];
			K->S[i][j]=Kc[j];
            posx++;
            
		}
	}
	else if(M[j].pos==2){
		posy=M[j].syi;
		for(i=0;i<M[j].s;i++){
			aux=pindex(M[j].sxi,posy,M[j].szi,N);
			KI[aux]=karst_index[j];
			posy++;
			KI_c[i][j]=karst_index[j]/(Ac[j]*hs[j]);
			p_aux[i][j]=p_old[aux];
			K->S[i][j]=Kc[j];
		}
	}
	else if(M[j].pos==3){
		posz=M[j].szi;
		for(i=0;i<M[j].s;i++){
			aux=pindex(M[j].sxi,M[j].syi,posz,N);
			KI[aux]=karst_index[j];
			posz++;
			KI_c[i][j]=karst_index[j]/(Ac[j]*hs[j]);
			p_aux[i][j]=p_old[aux];
			K->S[i][j]=Kc[j];
		}
	}
}


conduit(p_aux,cont_karst,KI_c,viscosidade,M,hs,indcol1,BoundKarst,keff,K,perm,p_old,pc);

double error, *error1;
double *Pc,*Pc_old;
Pc=create_vector(Mmax);
Pc_old=create_vector(Mmax);
error1=create_vector(number_karst+1);
error =1;

int iter=0;

while((iter<iter_max) && (error>tolerancia)){

	pc_old=atuaMatriz(pc_old,pc,Mmax,number_karst);
	reservoir(N, K, WI, KI,pc_old,indcol,B,perm, h,boundary,viscosidade,pw,M,p);
 
    
	error1[number_karst]=erro_relativo(p,p_old,N.t);
    p_old=atuaVector(p_old,p,N.t);


    for(j=0;j<number_karst;j++){
    if(M[j].pos==1){
        posx=M[j].sxi;
        for(i=0;i<M[j].s;i++){
            aux=pindex(posx,M[j].syi,M[j].szi,N);
            posx++;
            p_aux[i][j]=p_old[aux];
        }
    }
    else if(M[j].pos==2){
        posy=M[j].syi;
        for(i=0;i<M[j].s;i++){
            aux=pindex(M[j].sxi,posy,M[j].szi,N);
            posy++;
            p_aux[i][j]=p_old[aux];
        }
    }
    else if(M[j].pos==3){
        posz=M[j].szi;
        for(i=0;i<M[j].s;i++){
            aux=pindex(M[j].sxi,M[j].syi,posz,N);
            posz++;
            p_aux[i][j]=p_old[aux];
            }
    }
    }



    conduit(p_aux,cont_karst,KI_c,viscosidade,M,hs,indcol1,BoundKarst,keff,K,perm,p_old,pc);
    for(i=0;i<number_karst;i++){
        for(j=0;j<M[j].s;j++){
            Pc[j]=pc[j][i];
            Pc_old[j]=pc_old[j][i];
        }
        error1[i]=erro_relativo(Pc,Pc_old,Mmax);
    }
 
    double max=-1;
    for(i=0;i<number_karst+1;i++){
        if(error1[i]>max){
            error=error1[i];
        }
    }
    
	iter++;

}



double tf=get_wall_time();

tf=tf-ti;

printf("iterations number=%d\n",iter);
printf("error=%e\n",error);
printf("cpu time=%f\n",tf);
fluxo1d(pc,keff,uc,hs,M,cont_karst,BoundKarst,K,perm,p_old);
save_pressure_1d(pc,M,hs,h);
save_flux_1d(uc,M,hs,h);
fluxos(p,perm,v,B,h,N);
save_pressure(p,N,h);
save_flux(v,N,h);


char fname[1024];
    FILE *ff;	
   
    /* IP COMPUTATION */
    double Q = 0.0;
    double Vw = Nwi*h.x*h.y*h.z; /* Volume of well: Change when necessary */
    double Ve = h.x*h.y*h.z; // Volume of element 
    double pavg = 0.0, Qborda = 0.0;                                 
    double mi = viscosidade;
    double QR = 0.0, QL = 0.0, QU = 0.0, QD = 0.0, QT = 0.0, QB = 0.0;
    int ii, jj, kk, ss;
    double pmedio=0.0;	

 
for (i=0; i < N.t; i++){
   Q = Q + (WI[i]/mi)*(pw[i] - p[i]);
   pavg = pavg + IDwi[i]*p[i]*Ve/Vw;
    pmedio=pmedio+p[i];
} 
pmedio=pmedio/N.t;


printf("pmedio=%f\n",pmedio); // average pressure of matrix

double Qc3d[number_karst], Qc[number_karst];
double Qcb[number_karst];
double **fon_karst;
fon_karst=create_matrix(Mmax,number_karst);

int dirx,diry,dirz;

FILE *fg;


for(j=0;j<number_karst;j++){

sprintf(fname,"%s/fonte_karst_%i.dat",dirname,j);
fg=fopen(fname,"w");
    Qc3d[j]=0.0;
    //Qcb[j]=0.0;
    Qc[j]=0.0;

dirx=M[j].sxi;
diry=M[j].syi;
dirz=M[j].szi;

for(i=0;i<M[j].s;i++){

                if(M[j].pos==1){
                aux = pindex(dirx,M[j].syi,M[j].szi,N);
                dirx++;  
                }

                else if(M[j].pos==2){
                aux = pindex(M[j].sxi,diry,M[j].szi,N);
                diry++;  
                }

                else if(M[j].pos==3){                
                aux = pindex(M[j].sxi,M[j].syi,dirz,N);
                dirz++;  
                }

                                
                Qc[j]=Qc[j]+(-KI[aux]/mi)*(pc[i][j]-p[aux])*hs[j];
                Qc3d[j]=Qc3d[j]+(KI[aux]/mi)*(pc[i][j]-p[aux]);
                fon_karst[i][j]=(KI[aux]/mi)*(pc[i][j]-p[aux]);                       
                
}              
Qcb[j]=uc->l[0][j]*Ac[j]+uc->r[M[j].s-1][j]*Ac[j];
printf("Qcb[%d]=%e,Qc1d[%d]=%e\n,",j,Qcb[j],j,Qc[j]);
}
           
            double Qcb_total=0.0;
            double Qc_total=0.0;
            double Qc3d_total=0.0;

            for(j=0;j<number_karst;j++){
                Qcb_total=Qcb_total+Qcb[j];
                Qc_total=Qc_total+Qc[j];
                Qc3d_total=Qc3d_total+Qc3d[j];
            }


    /* Face R (right) e L (left) */
	for (ii=0;ii<N.y;ii++){
        for (kk=0;kk<N.z;kk++){

            ss = pindex(N.x-1,ii,kk,N);
			QR = QR + v[ss].r*h.y*h.z;

			ss = pindex(0,ii,kk,N);
			QL = QL + v[ss].l*h.y*h.z;		
       	}
    }

    /* Face U (up) e D (down) */
	for (ii=0;ii<N.x;ii++){
        for (kk=0;kk<N.z;kk++){

            ss = pindex(ii,N.y-1,kk,N);
		    QU = QU + v[ss].u*h.x*h.z;		

			ss = pindex(ii,0,kk,N);
			QD = QD + v[ss].d*h.x*h.z;		
        }
    }

    /* Face T (top) e B (bottom) */
	for (ii=0;ii<N.x;ii++){
        for (kk=0;kk<N.y;kk++){

            ss = pindex(ii,kk,N.z-1,N);
		    QT = QT + v[ss].t*h.x*h.y;		

			ss = pindex(ii,kk,0,N);
			QB = QB + v[ss].b*h.x*h.y;		
        }
    }


  
   
double Dp1=(Ppoco-pmedio);
double IP1;  
IP1=Q/Dp1;
    
printf("IP=%e\n", IP1);

/* Total da vazão nas bordas */
Qborda = QR + QL + QU + QD + QT + QB;

printf("Q=%e\n", Q);
printf("Q+Qc=%e , Qb=%e \n", Q+Qc3d_total,Qborda);
printf("Qcb=%e \n", Qcb_total);
printf("Qc3d=%e \n", Qc3d_total);
printf("Qc1d=%e\n", Qc_total);
        

double *P_pc;
ELE *u_uc;
u_uc=(ELE*)malloc((N.t+1)*sizeof(ELE));
P_pc=create_vector(N.t);

for(j=0;j<N.t;j++){
    u_uc[j].l=v[j].l;
    u_uc[j].r=v[j].r;
    u_uc[j].u=v[j].u;
    u_uc[j].d=v[j].d;
    u_uc[j].b=v[j].b;
    u_uc[j].t=v[j].t;
    P_pc[j]=p[j];
}

for(j=0;j<number_karst;j++){

    if(M[j].pos==1){
        dirx=M[j].sxi;
        for(i=0;i<M[j].s;i++){
        aux = pindex(dirx,M[j].syi,M[j].szi,N);
        P_pc[aux]=pc[i][j];
        u_uc[aux].l=uc->l[i][j];
        u_uc[aux].r=uc->r[i][j];
        dirx++;
        }
    }

    else if(M[j].pos==2){
        diry=M[j].syi;
        for(i=0;i<M[j].s;i++){        
        aux = pindex(M[j].sxi,diry,M[j].szi,N);
        P_pc[aux]=pc[i][j];           
        u_uc[aux].d=uc->l[i][j];
        u_uc[aux].u=uc->r[i][j];
        diry++;
        }
    }

    else if(M[j].pos==3){  
        dirz=M[j].szi;
        for(i=0;i<M[j].s;i++){              
        aux = pindex(M[j].sxi,M[j].syi,dirz,N);
        u_uc[aux].b=uc->l[i][j];
        u_uc[aux].t=uc->r[i][j];
        P_pc[aux]=pc[i][j]; 
        dirz++;
        }
    }


}

vts_StructuredGrid_file(u_uc,P_pc,N,h);


 if(B !=NULL){
        free(B->u);
        free(B->d);
        free(B->l);
        free(B->r);
        free(B->t);
        free(B->b); 
        free(B);
 }

  if(BoundKarst !=NULL){
        free(BoundKarst->l);
        free(BoundKarst->r);      
        free(BoundKarst);
 }

free(keff);
free(cont_karst);
free(C_contorno);
free(n_brechas);
free(hs);
free(pw);
free(IDwi);
free(ucR_old);
free(ucL_old),
free(Pc);
free(Pc_old);
free(error1);
free(fon_karst);
free(u_uc);
free(P_pc);
free(indcol1);
free(uc);
free(pc);
free(pc_old);
free(M);
free(Kc);


if(K !=NULL){
	free(K->X);
	free(K->Y);
    free(K->Z);
    free(K->S);

}

free(perm);
free(indcol);
free(v);
free(p);
free(p_old);
free(KI);
free(KI_c);
free(WI);


    return 0;

} 
