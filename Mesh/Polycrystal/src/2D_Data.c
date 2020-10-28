#include <stdio.h>
#include <stdlib.h>
#include "Global_variables.h"
#include <math.h>
/*****************************************************************************/
/*                                                                           */
/*  funtion to compute the transformation for faces                          */
/*                                                                           */
/*****************************************************************************/
extern struct external_data input, output;
float transform_vector(float *normal,float *vector,float *vector_t);
float transform_vector_inverse(float *normal,float *vector,float *vector_t);

void Planar_points_segments()
{

    int i, j, nvertices, vertice, vertice_1, vertice_2, count, ns, no1, no2;
    int nel, elem_struc, node1, node2, ii, ind, nel_local, k, nvertices_total;
    float x1, y1, x2, y2, total_L, L, Le_bar, Le, delta_x, delta_y;
    float x_i, y_i, z_i;
    float normal[3], vector[3], vector_t[3], Le_bar_face;
      
    printf("\n    3. Transformation into 2D FACES_struc\n\n");
    
    // Matrix of VERTICES_local_struc
    input.VERTICES_local_struc = (float *) malloc(input.ntotal_vertices*4*sizeof(float));
    nvertices_total = 0;
    for (i=0; i<input.ntotal_faces; i++){ // loop over faces
        // normal vector in ith face
        normal[0] = input.NORMAL_VECTORS[i*4+1]; //x   
        normal[1] = input.NORMAL_VECTORS[i*4+2]; //y
        normal[2] = input.NORMAL_VECTORS[i*4+3]; //z
        nvertices = input.vertices_face[i*2+1]; // vertices in ith face
        for (j=0; j<nvertices; j++){ // loop over vertices in ith face
            vertice = input.FACES_struc[i*(max_vertices_global+1)+j+1];
            // vector of coordinates to be transform
            vector[0] = input.VERTICES_struc[vertice*4+1]; //x 
            vector[1] = input.VERTICES_struc[vertice*4+2]; //y
            vector[2] = input.VERTICES_struc[vertice*4+3]; //z      

            transform_vector(normal,vector,vector_t);

            input.VERTICES_local_struc[vertice*4] = vertice;
            input.VERTICES_local_struc[vertice*4+1] = vector_t[0];//x
            input.VERTICES_local_struc[vertice*4+2] = vector_t[1];//y
            input.VERTICES_local_struc[vertice*4+3] = vector_t[2];//z
            nvertices_total = nvertices_total + 1;
        }
    }   

    float vector_L[nvertices_total];

    // Mesh parameters, number of divisions (ns)
    // and matrix of division per segment divisions_segment[]
    input.average_face = (float *) malloc(input.ntotal_faces*2*sizeof(float));
    count=0;
    total_L=0;
    for (i=0; i<input.ntotal_faces; i++){ // loop over faces
        nvertices = input.vertices_face[i*2+1]; // vertices in ith face
        // total_L=0;
        for (j=0; j<nvertices-1; j++){    
            // initial and final points of a segment
            vertice_1 = input.FACES_struc[i*(max_vertices_global+1)+j+1]; 
            vertice_2 = input.FACES_struc[i*(max_vertices_global+1)+j+2];    
            // Coordinates of initial and final points 
            x1 = input.VERTICES_local_struc[vertice_1*4+1];//x of vertice 1
            y1 = input.VERTICES_local_struc[vertice_1*4+2];//y of vertice 1
            x2 = input.VERTICES_local_struc[vertice_2*4+1];//x of vertice 2
            y2 = input.VERTICES_local_struc[vertice_2*4+2];//y of vertice 2
            // length            
            L = sqrt(pow(x1-x2,2)+pow(y1-y2,2));  
            //vector_L[count] = L;  
            //count = count + 1;
            total_L = total_L + L;
        }
        vertice_1 = input.FACES_struc[i*(max_vertices_global+1)+(j-1)+2]; 
        vertice_2 = input.FACES_struc[i*(max_vertices_global+1)+1];
        // Coordinates of initial and final points 
        x1 = input.VERTICES_local_struc[vertice_1*4+1];//x of vertice 1
        y1 = input.VERTICES_local_struc[vertice_1*4+2];//y of vertice 1
        x2 = input.VERTICES_local_struc[vertice_2*4+1];//x of vertice 2
        y2 = input.VERTICES_local_struc[vertice_2*4+2];//y of vertice 2
        // length            
        L = sqrt(pow(x1-x2,2)+pow(y1-y2,2));  
        //vector_L[count] = L;
        //count = count + 1;
  
        total_L = total_L + L;
        //Average length for ith face
        //Le_bar = total_L/nvertices;
        //Amax = 0.5*pow((Le_bar/dm),2); 
        //input.average_face[i*2] = Le_bar;
        //input.average_face[i*2+1] = Amax;
    }
    
    // Segment with the minium length
    //float minor=vector_L[0];
    //for (i=0; i<nvertices_total; i++){
    //   if(vector_L[i]>0){
    //        if(vector_L[i]<minor){
    //            minor = vector_L[i];
    //        }
    //    }
    //}

    Le_bar = total_L/nvertices_total;
	
    // Amax = 0.5*pow((Le_bar/dm),2); 
    //printf("Le_bar = %f\n",Le_bar);

    input.ntotal_points = 0;
    count = -1;
    input.divisions_segment = (int *) malloc(input.ntotal_vertices*2*sizeof(int));
    input.SEGMENTS_struc = (int *) malloc(input.ntotal_vertices*3*sizeof(int));
    
    for (i=0; i<input.ntotal_faces; i++){ // loop over faces
        nvertices = input.vertices_face[i*2+1]; // vertices in ith face
        
        //Average length for ith face
        //Le_bar = input.average_face[i*2];
        for (j=0; j<nvertices-1; j++){    
            // initial and final points of a segment
            vertice_1 = input.FACES_struc[i*(max_vertices_global+1)+j+1]; 
            vertice_2 = input.FACES_struc[i*(max_vertices_global+1)+j+2];    
            // Coordinates of initial and final points 
            x1 = input.VERTICES_local_struc[vertice_1*4+1];//x of vertice 1
            y1 = input.VERTICES_local_struc[vertice_1*4+2];//y of vertice 1
            x2 = input.VERTICES_local_struc[vertice_2*4+1];//x of vertice 2
            y2 = input.VERTICES_local_struc[vertice_2*4+2];//y of vertice 2
            // length of segment          
            Le = sqrt(pow(x1-x2,2)+pow(y1-y2,2));
            // Numer of division on a segment
            ns = dm*round(Le/Le_bar + 1);
            //ns = ns + 3;
            //printf("ns=%d\n",ns);
            // divisions per segment
            count = count + 1;
            input.ntotal_points = input.ntotal_points + ns;
            if (ns==0){
                input.ntotal_points = input.ntotal_points + 1;
            }   
            input.divisions_segment[(count)*2] = count;
            input.divisions_segment[(count)*2+1] = ns;
            // SEGMENTS_struc
            input.SEGMENTS_struc[count*3] = count;
            input.SEGMENTS_struc[count*3+1] = vertice_1;
            input.SEGMENTS_struc[count*3+2] = vertice_2;   
      
        }
        vertice_1 = input.FACES_struc[i*(max_vertices_global+1)+(j-1)+2]; 
        vertice_2 = input.FACES_struc[i*(max_vertices_global+1)+1];
        // Coordinates of initial and final points 
        x1 = input.VERTICES_local_struc[vertice_1*4+1];//x of vertice 1
        y1 = input.VERTICES_local_struc[vertice_1*4+2];//y of vertice 1
        x2 = input.VERTICES_local_struc[vertice_2*4+1];//x of vertice 2
        y2 = input.VERTICES_local_struc[vertice_2*4+2];//y of vertice 2
        // length            
        Le = sqrt(pow(x1-x2,2)+pow(y1-y2,2));     
        // Numer of division on a segment
        ns = dm*round(Le/Le_bar + 1);
        // divisions per segment
        count = count + 1;
        if (ns==0){
            input.ntotal_points = input.ntotal_points + 1;
        }
        input.ntotal_points = input.ntotal_points + ns;
        input.divisions_segment[(count)*2] = count;
        input.divisions_segment[(count)*2+1] = ns;
        // SEGMENTS_struc
        input.SEGMENTS_struc[count*3] = count;
        input.SEGMENTS_struc[count*3+1] = vertice_1;
        input.SEGMENTS_struc[count*3+2] = vertice_2; 
    }
  
    // Matrix of POINTS_struc and ELEM_struc
    elem_struc = -1;
    node2 = 0;
    count = 0;
    input.POINTS_struc = (float *) malloc(input.ntotal_points*4*sizeof(float));
    input.ELEM_struc = (int *) malloc(input.ntotal_points*3*sizeof(int)); 
    input.divisions_face = (int *) malloc(input.ntotal_faces*2*sizeof(int));
    for (ii=0; ii<input.ntotal_faces; ii++){ // loop over faces
        nvertices = input.vertices_face[ii*2+1]; // vertices in ith face
        nel_local = 0;
        for (i=0; i<nvertices; i++){ // loop over segments in the ith face
            // Nodes 1 and 2
            no1 = input.SEGMENTS_struc[count*3+1]; no2 = input.SEGMENTS_struc[count*3+2];
            // Coordinates of node 1
            x1 = input.VERTICES_local_struc[no1*4+1]; y1 = input.VERTICES_local_struc[no1*4+2];
            // Coordinates of node 2
            x2 = input.VERTICES_local_struc[no2*4+1]; y2 = input.VERTICES_local_struc[no2*4+2];
            // space
            delta_x = x2-x1; delta_y = y2-y1;
            // number of elements in a segment
            nel = input.divisions_segment[(count)*2+1];
            nel_local = nel_local + nel;
            // z coordinate
            count = count + 1;
            z_i = input.VERTICES_local_struc[no1*4+3];
            if (nel>0){
                for (j=0; j<nel; j++){ // loop over divisions
                    x_i = x1 + (delta_x/nel)*j;
                    y_i = y1 + (delta_y/nel)*j;            
                    // Matrix of POINTS_struc
                    elem_struc = elem_struc + 1;
                    input.POINTS_struc[elem_struc*4] = elem_struc;
                    input.POINTS_struc[elem_struc*4+1] = x_i;
                    input.POINTS_struc[elem_struc*4+2] = y_i;
                    input.POINTS_struc[elem_struc*4+3] = z_i;
                    // Matrix of ELEM_struc
                    node1 = node2; node2 = node2 + 1;
                    input.ELEM_struc[elem_struc*3] = elem_struc;
                    input.ELEM_struc[elem_struc*3+1] = node1;
                    input.ELEM_struc[elem_struc*3+2] = node2; 
                }
            }
            if (nel==0){
                    nel_local = nel_local + 1;
                    // Matrix of POINTS_struc
                    elem_struc = elem_struc + 1;
                    input.POINTS_struc[elem_struc*4] = elem_struc;
                    input.POINTS_struc[elem_struc*4+1] = x1;
                    input.POINTS_struc[elem_struc*4+2] = y1;
                    input.POINTS_struc[elem_struc*4+3] = z_i;
                    // Matrix of ELEM_struc
                    node1 = node2; node2 = node2 + 1;
                    input.ELEM_struc[elem_struc*3] = elem_struc;
                    input.ELEM_struc[elem_struc*3+1] = node1;
                    input.ELEM_struc[elem_struc*3+2] = node2;
            }
                
        }
        input.divisions_face[ii*2] = ii;
        input.divisions_face[ii*2+1] = nel_local; 
        ind = elem_struc-nel_local+1;
        input.ELEM_struc[elem_struc*3+2] = input.ELEM_struc[ind*3+1];        

    }
    
    count = 0;
    elem_struc = -1;
    for (i=0; i<input.ntotal_faces; i++){ // loop over faces
        nvertices = input.vertices_face[i*2+1]; // vertices in ith face
        nel_local = 0;
        total_L = 0;
        for (j=0; j<nvertices; j++){ 
            nel = input.divisions_segment[(count)*2+1];
            count = count + 1;
            nel_local = nel_local + nel;
            for (k=0; k<nel; k++){ // loop over divisions

                elem_struc = elem_struc + 1;
  
                // initial and final points of a segment        
                node1 = input.ELEM_struc[elem_struc*3+1];
                node2 = input.ELEM_struc[elem_struc*3+2]; 

                // Coordinates of initial and final points 
                x1 = input.POINTS_struc[node1*4+1];//x of vertice 1
                y1 = input.POINTS_struc[node1*4+2];//y of vertice 1
                x2 = input.POINTS_struc[node2*4+1];//x of vertice 2
                y2 = input.POINTS_struc[node2*4+2];//y of vertice 2

                // length            
                L = sqrt(pow(x1-x2,2)+pow(y1-y2,2));    
                total_L = total_L + L;
            }
            
        }
        
        

        // Average length for ith face
        Le_bar = total_L/nel_local;
        Amax = 0.5*pow((Le_bar),2);  
        
        input.average_face[i*2+1] = Amax;
    }
    
    printf("'----------------------------------------------------------------------------'\n");
}


void ELEM_POINTS()
{

    int i, j, k, tam, ntriangles, count, count2, nvertices, segment, npoints;
    int nsegments, count3;
    float z;
    float normal[3], vector[3], vector_t[3];
    
    printf("\n    5. Matrices ELEM and POINTS to export\n\n");
    
    // ============================================================================
    // Matrix ELEM and POINTS =====================================================
    count = 0;
    count2 = 0;
    count3 = 0;
    segment=0;
    output.ELEM = (int *) malloc(output.ntotal_ELEM*5*sizeof(int));
    output.POINTS = (float *) malloc(output.ntotal_POINTS*4*sizeof(float));
    output.SEGMENTS = (int *) malloc(output.ntotal_SEGMENTS*3*sizeof(int));
    tam = input.ntotal_points;
    for (i=0; i<input.ntotal_faces; i++){ // loop over faces
        nvertices = input.vertices_face[i*2+1];
        z = input.VERTICES_local_struc[segment*4+3];
        segment = segment + nvertices;
        ntriangles = output.Store_ELEM[i*tam];
        nsegments = output.store_SEGMENTS[i*tam];
        for (j=0; j<ntriangles; j++){ // loop over triangles
            output.ELEM[count*5] = count;                            // element
            output.ELEM[count*5+1] = i;                              // face
            output.ELEM[count*5+2] = output.Store_ELEM[i*tam+j*3+2]; // node1
            output.ELEM[count*5+3] = output.Store_ELEM[i*tam+j*3+3]; // node2
            output.ELEM[count*5+4] = output.Store_ELEM[i*tam+j*3+4]; // node3
            count = count + 1;
        }

        for (j=0; j<nsegments; j++){ // loop over segments
            output.SEGMENTS[count3*3] = count3;
            output.SEGMENTS[count3*3+1] = output.store_SEGMENTS[i*tam+j*2+2];
            output.SEGMENTS[count3*3+2] = output.store_SEGMENTS[i*tam+j*2+3];
            count3 = count3 + 1;
        }        

        // normal vector in ith face
        normal[0] = input.NORMAL_VECTORS[i*4+1]; //x   
        normal[1] = input.NORMAL_VECTORS[i*4+2]; //y
        normal[2] = input.NORMAL_VECTORS[i*4+3]; //z

        npoints = output.Store_POINTS[i*tam];
        for (j=0; j<npoints; j++){ // loop over points    
            
            // vector of coordinates to be transform
            vector[0] = output.Store_POINTS[i*tam+j*2+1];  //x 
            vector[1] = output.Store_POINTS[i*tam+j*2+2];  //y
            vector[2] = z;                                  //z    

            transform_vector_inverse(normal,vector,vector_t);
            
            output.POINTS[count2*4] = count2;      //point               
            output.POINTS[count2*4+1] = vector_t[0]; // x  
            output.POINTS[count2*4+2] = vector_t[1]; // y 
            output.POINTS[count2*4+3] = vector_t[2]; // z                           
            count2 = count2 + 1;

        }

        
 
    }
    // ============================================================================
    // ============================================================================

    printf("'----------------------------------------------------------------------------'\n");

}
