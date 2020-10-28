/*****************************************************************************/
/*                                                                           */
/*  (tricall.c)                                                              */
/*                                                                           */
/*  Example program that demonstrates how to call Triangle.                  */
/*                                                                           */
/*  Accompanies Triangle Version 1.6                                         */
/*  July 19, 1996                                                            */
/*                                                                           */
/*  This file is placed in the public domain (but the file that it calls     */
/*  is still copyrighted!) by                                                */
/*  Jonathan Richard Shewchuk                                                */
/*  2360 Woolsey #H                                                          */
/*  Berkeley, California  94705-1927                                         */
/*  jrs@cs.berkeley.edu                                                      */
/*                                                                           */
/*****************************************************************************/

/* If SINGLE is defined when triangle.o is compiled, it should also be       */
/*   defined here.  If not, it should not be defined here.                   */

/* #define SINGLE */

#ifdef SINGLE
#define REAL float
#else /* not SINGLE */
#define REAL double
#endif /* not SINGLE */

#include <stdio.h>
#include <stdlib.h>
#include "triangle.h"
#include "Global_variables.h"

extern void Input_data_struc();
extern void Planar_points_segments();
extern void ELEM_POINTS();
extern void Output_file();
extern void Output_file3();
struct external_data input, output;
int max_faces, max_vertices_global, dm;
double Amax;
/*****************************************************************************/
/*                                                                           */
/*  main()   Create and refine a mesh.                                       */
/*                                                                           */
/*****************************************************************************/

int main()
{
    
    int i, j, k, ndivisions, count, count2, tam, npoints, nfaces_reg;
    int nfaces, npoints2, face, face_A, face_B, reg_A, reg_B, ii;
    int ntriangles, num_points, nsegs_in, nsegs_out;
    
    char parameters[50];

    //  funtion to read data from the voro structure   
    Input_data_struc();
    struct triangulateio in, mid, out, vorout;
    
    // Definition of input parameters ------------------------------------------

    Planar_points_segments();
    
    printf("\n    4. Triangulate each face\n\n");
    
    // Triangulation of each face ----------------------------------------------
    
    npoints=0;
    npoints2=0;
    output.ntotal_POINTS = 0;
    output.ntotal_ELEM = 0;
    output.ntotal_SEGMENTS = 0;
    tam = input.ntotal_points;
    output.Store_ELEM = (int *) malloc(input.ntotal_faces*tam*sizeof(int));
    output.Store_POINTS = (float *) malloc(input.ntotal_faces*tam*sizeof(float));
    output.FACES = (int *) malloc(input.ntotal_faces*2*sizeof(int));
    output.store_SEGMENTS = (int *) malloc(input.ntotal_faces*tam*sizeof(int));
    
    count = 0;
    for (i=0; i<input.ngrains; i++){ // Regions
        for (j=0; j<input.faces_region[i*2+1]; j++){ // Faces per region
            face = input.SUBREGIONS[i*(max_faces+1)+j+1]; // ith face
            for (k=0; k<input.ntotal_interfaces; k++){ // interfaces
                face_A = input.INTERFACES_struc[k*6+3]; // Face_A
                face_B = input.INTERFACES_struc[k*6+4]; // Face_B
                // interface
                if (face == face_A){ 
                    reg_A = input.INTERFACES_struc[k*6+1]; // Reg_A
                    reg_B = input.INTERFACES_struc[k*6+2]; // Reg_B
                    if (reg_A < reg_B){
                        
                        ndivisions = input.divisions_face[face_A*2+1]; //division in face_A                    
                        
                        // Definition of input parameter of Triangulate -----------------------
                        in.numberofpoints = ndivisions;
                        in.numberofpointattributes = 0;
                        in.pointlist = (REAL *) malloc(in.numberofpoints * 2 * sizeof(REAL));
                        in.pointmarkerlist = (int *) NULL;
                        in.pointattributelist = (REAL *) NULL;

                        in.numberofsegments = ndivisions;
                        in.segmentlist = (int *) malloc(in.numberofsegments * 2 * sizeof(int));
                        in.segmentmarkerlist = (int *) NULL;

                        in.numberofholes = 0;
                        in.holelist = (REAL *) NULL;

                        in.numberofregions = 0;
                        in.regionlist = (REAL *) NULL;

                        // Definition of output parameters ------------------------------------
                        out.pointlist = (REAL *) NULL; 
                        out.pointmarkerlist = (int *) NULL;
                              
                        out.segmentlist = (int *) NULL; 
                        out.segmentmarkerlist = (int *) NULL;      

                        out.trianglelist = (int *) NULL;
                        // --------------------------------------------------------------------                        
                        
                        count2=0;
                        for (ii=0; ii<ndivisions; ii++){ // loop over divisions

                            // in.segmentlist !!
                            in.segmentlist[ii*2] = count2; // initial node
                            count2 = count2 + 1;
                            in.segmentlist[ii*2+1] = count2; // final node
                          
                            // in.pointlist !!
                            in.pointlist[ii*2] = input.POINTS_struc[count*4+1]; // x
                            in.pointlist[ii*2+1] = input.POINTS_struc[count*4+2]; // y

                            count = count + 1; 
                                            
                        }
                        in.segmentlist[(ii-1)*2+1] = 0; // final node  
                        // --------------------------------------------------------------------    
                        // Triangulate ith face
                        Amax = input.average_face[face_A*2+1];
                        sprintf(parameters, "pqza%.16fQ",Amax);
                        //printf("aqui\n");            
                        triangulate(parameters, &in, &out, (struct triangulateio *) NULL);
                        
                        // --------------------------------------------------------------------    
                        // Saving elements in a global matrix
                        if (face_A==0){
                            output.Store_ELEM[face_A*tam] = out.numberoftriangles;
                            output.Store_ELEM[face_A*tam+1] = out.numberofpoints;
                            for (ii=0; ii<out.numberoftriangles*3; ii++){
                                output.Store_ELEM[face_A*tam+ii+2] = out.trianglelist[ii];
                            }
                        }else{
                            output.Store_ELEM[face_A*tam] = out.numberoftriangles;
                            output.Store_ELEM[face_A*tam+1] = out.numberofpoints;
                            npoints = npoints + output.Store_ELEM[(face_A-1)*tam+1];
                            for (ii=0; ii<out.numberoftriangles*3; ii++){
                                output.Store_ELEM[face_A*tam+ii+2] = out.trianglelist[ii]+npoints;
                            }
                        }

                        // Saving points in a global matrix
                        output.Store_POINTS[face_A*tam] = out.numberofpoints;
                        for (ii=0; ii<out.numberofpoints*2; ii++){
                            output.Store_POINTS[face_A*tam+ii+1] = out.pointlist[ii];
                        }

                        // Matrix of SEGMENTS 
                        if (face_A==0){
                            output.store_SEGMENTS[face_A*tam] = in.numberofsegments; 
                            output.store_SEGMENTS[face_A*tam+1] = out.numberofsegments;   
                            for (ii=0; ii<in.numberofsegments*2; ii++){  
                                output.store_SEGMENTS[face_A*tam+ii+2] = in.segmentlist[ii];    
                            }
                        }else{
                            output.store_SEGMENTS[face_A*tam] = in.numberofsegments;  
                            output.store_SEGMENTS[face_A*tam+1] = out.numberofsegments;  
                            npoints2 = npoints2 + output.store_SEGMENTS[(face_A-1)*tam+1];
                            for (ii=0; ii<in.numberofsegments*2; ii++){  
                                output.store_SEGMENTS[face_A*tam+ii+2] = in.segmentlist[ii]+npoints2;    
                            }
                        }

                        // ============================================================================
                        // Matrix of FACES ============================================================
                        output.FACES[face_A*2] = face_A; // face
                        output.FACES[face_A*2+1] = out.numberoftriangles; // elements in ith face
                        // ============================================================================
                        // ============================================================================
            
                        output.ntotal_POINTS = output.ntotal_POINTS + out.numberofpoints;
                        output.ntotal_ELEM = output.ntotal_ELEM + out.numberoftriangles;
                        output.ntotal_SEGMENTS = output.ntotal_SEGMENTS + in.numberofsegments;
                        
                        // Free memonry ---------------------------------------------------------
                        free(in.pointlist); 
                        free(in.pointattributelist);
                        free(in.pointmarkerlist);
                        free(in.segmentlist);
                        free(in.segmentmarkerlist);
                        free(in.holelist);
                        free(in.regionlist);

                        free(out.pointlist);
                        free(out.pointmarkerlist);
                        free(out.segmentlist);
                        free(out.segmentmarkerlist);
                        free(out.trianglelist);
                        
                        break;

                    }else{
                        
                        ndivisions = input.divisions_face[face_B*2+1]; //division in face_A 
                        count = count + ndivisions;
                        
                        // Saving elements in a global matrix
                        ntriangles = output.Store_ELEM[(face_B)*tam]; // triangles in the face_B
                        num_points = output.Store_ELEM[(face_A-1)*tam+1]; // points in the face_B

                        output.Store_ELEM[face_A*tam] = ntriangles;
                        output.Store_ELEM[face_A*tam+1] = num_points; 
                        for (ii=0; ii<ntriangles*3; ii++){
                            output.Store_ELEM[face_A*tam+ii+2] = output.Store_ELEM[face_B*tam+ii+2];
                        }
                        
                            
                        // Matrix of SEGMENTS
                        nsegs_in =  output.store_SEGMENTS[face_B*tam]; // segments_in at face_B
                        nsegs_out = output.store_SEGMENTS[(face_A-1)*tam+1]; // segments_out at face_B
 
                        output.store_SEGMENTS[face_A*tam] = nsegs_in; 
                        output.store_SEGMENTS[face_A*tam+1] = nsegs_out;   
                        for (ii=0; ii<nsegs_in*2; ii++){  
                            output.store_SEGMENTS[face_A*tam+ii+2] = output.store_SEGMENTS[face_B*tam+ii+2];    
                        }

                        // ============================================================================
                        // Matrix of FACES ============================================================
                        output.FACES[face_A*2] = face_A; // face
                        output.FACES[face_A*2+1] = ntriangles; // elements in ith face
                        // ============================================================================
                        // ============================================================================
                        
                        output.ntotal_ELEM = output.ntotal_ELEM + ntriangles;
                        output.ntotal_SEGMENTS = output.ntotal_SEGMENTS + nsegs_in;
                            
                        break;
                    }
                }
                
                // Non-interface
                if ((k == input.ntotal_interfaces-1)&&(face != face_A)&&(face != face_B)){
                        
                        ndivisions = input.divisions_face[face*2+1]; //division in face                    
                             
                        // Definition of input parameter of Triangulate -----------------------
                        in.numberofpoints = ndivisions;
                        in.numberofpointattributes = 0;
                        in.pointlist = (REAL *) malloc(in.numberofpoints * 2 * sizeof(REAL));
                        in.pointmarkerlist = (int *) NULL;
                        in.pointattributelist = (REAL *) NULL;

                        in.numberofsegments = ndivisions;
                        in.segmentlist = (int *) malloc(in.numberofsegments * 2 * sizeof(int));
                        in.segmentmarkerlist = (int *) NULL;

                        in.numberofholes = 0;
                        in.holelist = (REAL *) NULL;

                        in.numberofregions = 0;
                        in.regionlist = (REAL *) NULL;

                        // Definition of output parameters ------------------------------------
                        out.pointlist = (REAL *) NULL; 
                        out.pointmarkerlist = (int *) NULL;
                              
                        out.segmentlist = (int *) NULL; 
                        out.segmentmarkerlist = (int *) NULL;      

                        out.trianglelist = (int *) NULL;
                        // --------------------------------------------------------------------

                        count2=0;
                        for (ii=0; ii<ndivisions; ii++){ // loop over divisions

                            // in.segmentlist !!
                            in.segmentlist[ii*2] = count2; // initial node
                            count2 = count2 + 1;
                            in.segmentlist[ii*2+1] = count2; // final node
                          
                            // in.pointlist !!
                            in.pointlist[ii*2] = input.POINTS_struc[count*4+1]; // x
                            in.pointlist[ii*2+1] = input.POINTS_struc[count*4+2]; // y

                            count = count + 1; 
                                            
                        }
                        in.segmentlist[(ii-1)*2+1] = 0; // final node  

                        // --------------------------------------------------------------------    
                        // Triangulate ith face
                        Amax = input.average_face[face*2+1];
                        sprintf(parameters, "pqza%.16fQ",Amax);
                         
                        triangulate(parameters, &in, &out, (struct triangulateio *) NULL);
                        // --------------------------------------------------------------------  

                        // Saving elements in a global matrix
                        if (face==0){
                            output.Store_ELEM[face*tam] = out.numberoftriangles;
                            output.Store_ELEM[face*tam+1] = out.numberofpoints;
                            for (ii=0; ii<out.numberoftriangles*3; ii++){
                                output.Store_ELEM[face*tam+ii+2] = out.trianglelist[ii];
                            }
                        }else{
                            output.Store_ELEM[face*tam] = out.numberoftriangles;
                            output.Store_ELEM[face*tam+1] = out.numberofpoints;
                            npoints = npoints + output.Store_ELEM[(face-1)*tam+1];
                            for (ii=0; ii<out.numberoftriangles*3; ii++){
                                output.Store_ELEM[face*tam+ii+2] = out.trianglelist[ii]+npoints;
                            }
                        }

                        // Saving points in a global matrix
                        output.Store_POINTS[face*tam] = out.numberofpoints;
                        for (ii=0; ii<out.numberofpoints*2; ii++){
                            output.Store_POINTS[face*tam+ii+1] = out.pointlist[ii];
                        }

                        // Matrix of SEGMENTS 
                        if (face==0){
                            output.store_SEGMENTS[face*tam] = in.numberofsegments; 
                            output.store_SEGMENTS[face*tam+1] = out.numberofsegments;   
                            for (ii=0; ii<in.numberofsegments*2; ii++){  
                                output.store_SEGMENTS[face*tam+ii+2] = in.segmentlist[ii];    
                            }
                        }else{
                            output.store_SEGMENTS[face*tam] = in.numberofsegments;  
                            output.store_SEGMENTS[face*tam+1] = out.numberofsegments;  
                            npoints2 = npoints2 + output.store_SEGMENTS[(face-1)*tam+1];
                            for (ii=0; ii<in.numberofsegments*2; ii++){  
                                output.store_SEGMENTS[face*tam+ii+2] = in.segmentlist[ii]+npoints2;    
                            }
                        }
                        
                        // ============================================================================
                        // Matrix of FACES ============================================================
                        output.FACES[face*2] = face; // face
                        output.FACES[face*2+1] = out.numberoftriangles; // elements in ith face
                        // ============================================================================
                        // ============================================================================
  
                        // Saving elements in a global matrix
                        
                        output.ntotal_POINTS = output.ntotal_POINTS + out.numberofpoints;
                        output.ntotal_ELEM = output.ntotal_ELEM + out.numberoftriangles;
                        output.ntotal_SEGMENTS = output.ntotal_SEGMENTS + in.numberofsegments;
                        
                        // Free memonry ---------------------------------------------------------
                        free(in.pointlist); 
                        free(in.pointattributelist);
                        free(in.pointmarkerlist);
                        free(in.segmentlist);
                        free(in.segmentmarkerlist);
                        free(in.holelist);
                        free(in.regionlist);

                        free(out.pointlist);
                        free(out.pointmarkerlist);
                        free(out.segmentlist);
                        free(out.segmentmarkerlist);
                        free(out.trianglelist);
                         
                }
            }
        }
    }
    
    // ============================================================================
    // Matrix of El_reg ===========================================================
    count = 0;
    output.El_reg = (int *) malloc(input.ngrains*2*sizeof(int));
    for (i=0; i<input.ngrains; i++){ // loop over region
        nfaces_reg = 0;
        nfaces = input.faces_region[i*2+1];
        for (j=0; j<nfaces; j++){ // loop over faces in ith region
            nfaces_reg = nfaces_reg + output.FACES[count*2+1];     
            count = count + 1;
        }  
        output.El_reg[i*2] = i; 
        output.El_reg[i*2+1] = nfaces_reg;
    }
    // ============================================================================
    // ============================================================================
    printf("'----------------------------------------------------------------------------'\n");
    
    // Matrices ELEM and POINTS

    ELEM_POINTS();
    
    // Write the Input_data.dat file 
    
    float *new_nv, nv_aux[3];
    int ind1,ind2, nel_face;

    new_nv =  (float *) malloc(output.ntotal_ELEM*4*sizeof(int));
    
    ind1 = 0;
    ind2 = 0;
    for (i=0; i<input.ntotal_faces; i++){ 
        nel_face = output.FACES[i*2+1];
        for (j=0; j<3; j++){
            nv_aux[j] = input.NORMAL_VECTORS[i*4+j+1];
        }
                
        ind2 = ind2 + nel_face;
        //printf("ind1=%d, ind2=%d\n",ind1,ind2);
        for (j=ind1; j<ind2; j++){
            for (k=0; k<3; k++){        
                new_nv[j*4+k+1] = nv_aux[k];
            }
            
        }

        ind1 = ind2;
        
    }

    free(input.NORMAL_VECTORS);
    input.NORMAL_VECTORS =  (float *) malloc(output.ntotal_ELEM*4*sizeof(int));

    for (i=0; i<output.ntotal_ELEM; i++){
        for (j=0; j<3; j++){
            input.NORMAL_VECTORS[i*4+j+1] = new_nv[i*4+j+1];           
        }
    }

    free(output.FACES);
    output.FACES =  (int *) malloc(output.ntotal_ELEM*2*sizeof(int));
    for (i=0; i<output.ntotal_ELEM; i++){
        output.FACES[i*2] = i;
        output.FACES[i*2+1] = 1;
        output.ELEM[i*5+1] = i;
    }


    Output_file();
	
	Output_file3();

  return 0;
}
