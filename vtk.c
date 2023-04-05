#include "structs.h"
#include "variables.h"
#include "functions.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int vts_StructuredGrid_file(ELE *v, double *x, Dim N, grid h)
{
  	// Print pressure and velocity in one VTK file ( structured grid, .vts)
  	int i,j,k,s;
  	double x_pos, y_pos, z_pos;
  	FILE *fn;
  	char fname[1024];


  	sprintf(fname, "%s/campos.vts",dirname);  
  	fn = fopen(fname, "w");
      
    if (fn == NULL)
      printf("Unable to open file %s.\n",fname);

      
    fprintf(fn,"<?xml version=\"1.0\"?>\n");
    fprintf(fn,"<VTKFile type=\"StructuredGrid\" version=\"0.1\">\n");  
    	fprintf(fn,"<StructuredGrid WholeExtent=\"0 %d 0 %d 0 %d\">\n",N.x-1,N.y-1,N.z-1);
    	    fprintf(fn,"<Piece Extent=\"0 %d 0 %d 0 %d\">\n",N.x-1,N.y-1,N.z-1);
    	      	fprintf(fn,"<PointData  Vectors=\"Velocity\" Scalars=\"pressure\">\n"); 
    	        	fprintf(fn,"<DataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    	        	// print velocity in points
    	        	  	for (k = 0; k< N.z; k++){
    	        	    	for (j = 0; j < N.y; j++){
    	        	      		for (i = 0; i < N.x; i++){
    	        	        		s = pindex(i,j,k,N);                          
    	        	        		fprintf(fn,"%1.8e %1.8e %1.8e\n",(0.5)*(v[s].r - v[s].l), (0.5)*(v[s].u - v[s].d),(0.5)*(v[s].t - v[s].b));
    	        	      		}
    	        	    	}	
    	        	  	}
    	        	fprintf(fn,"</DataArray>\n");
    	        	fprintf(fn,"<DataArray type=\"Float32\" Name=\"pressure\" NumberOfComponents=\"1\" format=\"ascii\">\n");
    	        	// print pressure in points
    	        	  for (k = 0; k< N.z; k++){
    	        	    for (j = 0; j < N.y; j++){
    	        	      for (i = 0; i < N.x; i++){
    	        	        s = pindex(i,j,k,N);                          
    	        	        fprintf(fn,"%1.8e\n",x[s]);
    	        	      }
    	        	    }
    	        	  }
    	        	fprintf(fn,"</DataArray>\n");
    	      	fprintf(fn,"</PointData>\n");
    	      	fprintf(fn,"<CellData>\n");
    	      	fprintf(fn,"</CellData>\n");
    	      	fprintf(fn,"<Points>\n");
    	        	fprintf(fn,"<DataArray type=\"Float32\" NumberOfComponents=\"3\">\n");
    	        // print the values of points x, y, z.
    	          		for (k = 0; k< N.z; k++){
    	            		z_pos = h.z*(k + 0.5);                        
    	              		for (j = 0; j < N.y; j++){
    	                		y_pos = h.y*(j + 0.5);
    	                		for (i = 0; i < N.x; i++){
    	                  			x_pos = h.x*(i + 0.5);
    	                  			fprintf(fn,"%.4f %.4f %.4f \n",x_pos,y_pos,z_pos);
    	                		}
    	              		}
    	            	}
    	        	fprintf(fn,"</DataArray>\n");
    	       	fprintf(fn,"</Points>\n");
    	    fprintf(fn,"</Piece>\n");
    	fprintf(fn,"</StructuredGrid>\n");
    fprintf(fn,"</VTKFile>\n");   

  fclose(fn);

  return 0;
}

int vts_StructuredGrid_flux_file(ELE *v, Dim N, grid h)
{	
	// print velocity only
	int i,j,k,s;
	double x_pos, y_pos,z_pos;
	FILE *fn;
	char fname[1024];

  	sprintf(fname, "%s/flux.vts",dirname);  
  	fn = fopen(fname, "w");
    	
    if (fn == NULL)
    	printf("Unable to open file %s.\n",fname);

    	
    fprintf(fn,"<?xml version=\"1.0\"?>\n");
    fprintf(fn,"<VTKFile type=\"StructuredGrid\" version=\"0.1\">\n");	
    	fprintf(fn,"<StructuredGrid WholeExtent=\"0 %d 0 %d 0 %d\">\n",N.x-1,N.y-1,N.z-1);
    		fprintf(fn,"<Piece Extent=\"0 %d 0 %d 0 %d\">\n",N.x-1,N.y-1,N.z-1);
	 			fprintf(fn,"<PointData Vectors=\"Velocity\">\n");
            		fprintf(fn,"<DataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n");
            		  for (k = 0; k< N.z; k++){
            		    for (j = 0; j < N.y; j++){
            		      for (i = 0; i < N.x; i++){
            		        s = pindex(i,j,k,N);                          
            		        fprintf(fn,"%1.8e %1.8e %1.8e\n",(0.5)*(v[s].r - v[s].l), (0.5)*(v[s].u - v[s].d),(0.5)*(v[s].t - v[s].b));
            		      }
            		    }
            		  }
            		fprintf(fn,"</DataArray>\n");
	 			  fprintf(fn,"</PointData>\n");
          		fprintf(fn,"<CellData>\n");
          		fprintf(fn,"</CellData>\n");
          		fprintf(fn,"<Points>\n");
	 				fprintf(fn,"<DataArray type=\"Float32\" NumberOfComponents=\"3\">\n");
              			for (k = 0; k< N.z; k++){
                			z_pos = h.z*(k + 0.5);  							        
    				  		for (j = 0; j < N.y; j++){
    				  			y_pos = h.y*(j + 0.5);
    				  			for (i = 0; i < N.x; i++){
                     				x_pos = h.x*(i + 0.5);
                        			fprintf(fn,"%.4f %.4f %.4f \n",x_pos,y_pos,z_pos);
    				  			}
    				  		}
    				  	}
	 				fprintf(fn,"</DataArray>\n");
	 			fprintf(fn,"</Points>\n");
	 		fprintf(fn,"</Piece>\n");
    	fprintf(fn,"</StructuredGrid>\n");
    fprintf(fn,"</VTKFile>\n");		

	fclose(fn);

	return 0;
}

int vts_StructuredGrid_pressure_file(double *x, Dim N, grid h)
{
  	// Print pressure only
  	int i,j,k,s;
  	double x_pos, y_pos, z_pos;
  	FILE *fn;
  	char fname[1024];

  	sprintf(fname, "%s/pressure.vts",dirname);  
  	fn = fopen(fname, "w");
      
    if (fn == NULL)
      printf("Unable to open file %s.\n",fname);
      
    fprintf(fn,"<?xml version=\"1.0\"?>\n");
    fprintf(fn,"<VTKFile type=\"StructuredGrid\" version=\"0.1\">\n");  
    	fprintf(fn,"<StructuredGrid WholeExtent=\"0 %d 0 %d 0 %d\">\n",N.x-1,N.y-1,N.z-1);
    		fprintf(fn,"<Piece Extent=\"0 %d 0 %d 0 %d\">\n",N.x-1,N.y-1,N.z-1);
          		fprintf(fn,"<PointData Scalars=\"pressure\">\n");
            		fprintf(fn,"<DataArray type=\"Float32\" Name=\"pressure\" NumberOfComponents=\"1\" format=\"ascii\">\n");
              			for (k = 0; k< N.z; k++){
                			for (j = 0; j < N.y; j++){
                 	 			for (i = 0; i < N.x; i++){
                    				s = pindex(i,j,k,N);                          
                    				fprintf(fn,"%1.8e\n",x[s]);
                  				}
                			}
              			}
            		fprintf(fn,"</DataArray>\n");
          		fprintf(fn,"</PointData>\n");
          		fprintf(fn,"<CellData>\n");
          		fprintf(fn,"</CellData>\n");
          		fprintf(fn,"<Points>\n");
            		fprintf(fn,"<DataArray type=\"Float32\" NumberOfComponents=\"3\">\n");
              			for (k = 0; k< N.z; k++){
                			z_pos = h.z*(k + 0.5);  							        
    				  		for (j = 0; j < N.y; j++){
    				  			y_pos = h.y*(j + 0.5);
    				  			for (i = 0; i < N.x; i++){
                     				x_pos = h.x*(i + 0.5);
                        			fprintf(fn,"%.4f %.4f %.4f \n",x_pos,y_pos,z_pos);
    				  			}
    				  		}
    				  	}
            		fprintf(fn,"</DataArray>\n");
           		fprintf(fn,"</Points>\n");
        	fprintf(fn,"</Piece>\n");
        fprintf(fn,"</StructuredGrid>\n");
    fprintf(fn,"</VTKFile>\n");   

  fclose(fn);

  return 0;
}

