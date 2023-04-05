# hmfem.3D_AMG

Solve Darcy's equation.

Serial Implementation
C, C++.

The code is mostly written in C.
The AMG solver used to solve the resulting linear system is written in C++.

## Problem

u = u(x), p = p(x), x \in Omega \subset R^d, d = {2,3}

u is the Darcy flux and p is the pressure.

Differential equation in mixed form:

div(u) = 0; u = -P grad(p);

where P = K * lambda

and lambda = 1 always in this code

K is the absolute permeability and lambda is the mobility.

**In the code, lambda is set as 1.0**

Complete the model with boundary conditions.

## Compile
The compilation is in a makefile

To compile:

$ make

### Folders
Object files are in 
/build/obj/

Executable file is in
/build/run/

**make sure that the folders run/ and obj/ exists**

AMG toolbox is in
/lib/libAMG/toolbox/


## Run

Run the script roda.sh

$ ./roda.sh

Inside this script it calls the script exe.sh: 

$ ./exe.sh nx ny nz Ix Iy Iz Perm type wi

Is it necessary to provide the following inputs:

The number of fine cells in each direction (x, y and z) for the subdomain (processor):

**nx**

**ny**

**nz**

The size of the subdomain:

**Ix**

**Iy**

**Iz**

Choose the permeability Type:

**Perm**

1 *homogeneous permeability: set as 1.0e-14 (can be changed inside function)*

The type: a variable for naming the output folder

**type**

The variable:

**wi** 

is the peacemann well value set to zero. **DO NOT CHANGE**

## Output

Pressure and Flux solution, as the output file (output.out) and error file (err.out) are in folder:
/AMG/AMG_Perm1_nxXnyXnz_IxXIyXIz_type/

**The final solution is given in a pressure inside each cell, and the velocities in each face of each
cell**

The variable is a struct named sol where:

sol[s].p is the pressure inside cell

sol[s].r is the velocity on the RIGHT face

sol[s].l is the velocity on the LEFT face

sol[s].u is the velocity on the UP face

sol[s].d is the velocity on the DOWN face

sol[s].t is the velocity on the TOP face

sol[s].b is the velocity on the BOTTOM face

sol[s].x is the velocity on x-direction (= (sol[s].r - sol[s].l)/2)

sol[s].y is the velocity on y-direction (= (sol[s].u - sol[s].d)/2)

sol[s].z is the velocity on z-direction (= (sol[s].t - sol[s].b)/2)

The point value is given by (x,y,z)

**s** is the number of the cell given by: **s = i + j * N + k * N * M**

where i = 0:nx-1, j = 0:ny-1, k = 0:nz-1, N = nx, M = ny;

The pressure file "pressure.dat" is in columns in the following order:

x y z sol[s].p

The flux file "flux.dat" is in columns in the following order:

x y z sol[s].x sol[s].y sol[s].z sol[s].r sol[s].l sol[s].u sol[s].d sol[s].t sol[s].b


