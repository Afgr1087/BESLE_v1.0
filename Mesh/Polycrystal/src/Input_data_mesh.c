#include <stdio.h>
#include <stdlib.h>
#include "Global_variables.h"
/*****************************************************************************/
/*                                                                           */
//  funtion to read data from the voro structure                             */
/*                                                                           */
/*****************************************************************************/


extern struct external_data input, output;

void Input_data_struc()
{
    FILE *fp;
    char* filename = "IO_files/packing4.dat"; 
	char* filename1 = "IO_files/packing10.dat"; 
	char* Mesh_vtk; 

    int i, j;
    
    //struct external_data input, output;
    fp = fopen(filename, "r");

	// reading mesh name    
	fscanf(fp, "%s", &input.Mesh_file);

    // reading ngrains and max_faces    
    fscanf(fp, "%d", &dm);
       
    // reading ngrains and max_faces    
    fscanf(fp, "%d %d", &input.ngrains, &max_faces);
    
    // reading matrix of SUBREGIONS and faces_region
    input.SUBREGIONS = (int *) malloc(input.ngrains *(max_faces+1)*sizeof(int));
    for (i=0; i<input.ngrains; i++){    
        for (j=0; j<max_faces+1; j++){
            fscanf(fp, "%d", &input.SUBREGIONS[i*(max_faces+1)+j]);
        }
    }
    input.faces_region = (int *) malloc(input.ngrains *2*sizeof(int));
    for (i=0; i<input.ngrains; i++){    
        for (j=0; j<2; j++){
            fscanf(fp, "%d", &input.faces_region[i*2+j]);
        }
    }
    
    // reading ntotal_faces and max_vertices_global   
    fscanf(fp, "%d %d", &input.ntotal_faces, &max_vertices_global);

    // reading matrix of FACES_struc and vertices_face
    input.FACES_struc = (int *) malloc(input.ntotal_faces *(max_vertices_global+1)*sizeof(int));
    for (i=0; i<input.ntotal_faces; i++){    
        for (j=0; j<max_vertices_global+1; j++){
            fscanf(fp, "%d", &input.FACES_struc[i*(max_vertices_global+1)+j]);
        }
    }
    input.vertices_face = (int *) malloc(input.ntotal_faces *2*sizeof(int));
    for (i=0; i<input.ntotal_faces; i++){    
        for (j=0; j<2; j++){
            fscanf(fp, "%d", &input.vertices_face[i*2+j]);
        }
    }

    // reading ntotal_vertices 
    fscanf(fp, "%d", &input.ntotal_vertices);

    // reading matrix of VERTICES_struc
    input.VERTICES_struc = (float *) malloc(input.ntotal_vertices*4*sizeof(float));
    for (i=0; i<input.ntotal_vertices; i++){    
        for (j=0; j<4; j++){
            fscanf(fp, "%f", &input.VERTICES_struc[i*4+j]);
        }
    }

    // reading matrix of NORMAL_VECTORS
    input.NORMAL_VECTORS = (float *) malloc(input.ntotal_faces*4*sizeof(float));
    for (i=0; i<input.ntotal_faces; i++){    
        for (j=0; j<4; j++){
            fscanf(fp, "%f", &input.NORMAL_VECTORS[i*4+j]);
        }
    }

    // reading ntotal_interfaces
    fscanf(fp, "%d", &input.ntotal_interfaces);

    // reading matrix of INTERFACES_struc
    input.INTERFACES_struc = (int *) malloc(input.ntotal_interfaces *6*sizeof(int));
    for (i=0; i<input.ntotal_interfaces; i++){    
        for (j=0; j<6; j++){
            fscanf(fp, "%d", &input.INTERFACES_struc[i*6+j]);
        }
    }

    // reading Volumes
    output.Volumes = (float *) malloc(input.ngrains*2*sizeof(float));
    for (i=0; i<input.ngrains; i++){ 
        for (j=0; j<2; j++){
            fscanf(fp, "%f", &output.Volumes[i*2+j]);
        }
    }


	fclose(fp); 
	
	

    //------------------------------------------------------------------------------------------//
    
    printf("'----------------------------------------------------------------------------'\n");
}

