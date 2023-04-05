#ifndef STRUCTS_H
#define STRUCTS_H


// Structure for the variables in each element:

struct ELEstruct
{
    double u;
    double d;
    double l;
    double r;
    double c;
    double t;
    double b;
};
typedef struct ELEstruct ELE;

struct ELEstruct_K
{
    double **l;
    double **r;
    double **c;
};
typedef struct ELEstruct_K ELEK;

struct Boundary
{
    double *u;
    double *d;
    double *l;
    double *r;
    double *t;
    double *b;
};
typedef struct Boundary Bound;

struct Index 
{
    int c;
    int u;
    int d;
    int l;
    int r;
    int t;
    int b;
};
typedef struct Index IND;

struct Faces 
{
    int u;
    int d;
    int l;
    int r;
    int t;
    int b;
};
typedef struct Faces FACE;

struct FacesKarst
{
	int *l;
	int *r;

};
typedef struct FacesKarst FACEK; 

// Structure for the absolute permeabilities:

struct Permeab 
{
    double *X;
    double *Y;
    double *Z;
    double **S;
    int permeability;
};
typedef struct Permeab Perm;

// Structures for grid sizes: 

struct Grid 
{
    double x;
    double y;
    double z;
};
typedef struct Grid grid;

// Structures for dimension sizes: 

struct Dimension 
{
    int x;
    int y;
    int z;
    int xy;
    int t;
    int px;
    int pz;
    int py;
};
typedef struct Dimension Dim;


struct Dimension_karst
{
int s;
int sxi;
int syi;
int szi;
int pos;
};
typedef struct Dimension_karst DimK;

struct Fonte
{
    double r;
    double c;
};
typedef struct Fonte fonte;

#endif
