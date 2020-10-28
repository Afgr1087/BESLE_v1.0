//==========================================================================================//
//------------------------------------------------------------------------------------------//
//                              MODULE FOR DATA ADQUISITION                                 //
//------------------------------------------------------------------------------------------//
//==========================================================================================//

#include <stdio.h>
#include <iostream>
#include <fstream>
#include "Input_data.hh"
#include <vector>
#include <string>
#include <sstream>
#include <stdlib.h>

using namespace std;

void Input_data_structure() // Data from the polycrystalline structure
{

    int reg, nfaces, nvertices;
    int count_faces, max_vertices, max_points;
    int ind_faces_ini, ind_faces_fin, count_vertices=0;
    int ind_1, ind_2;
    char val;
    printf("\n    2. Data adquisition from the polycrystalline structure \n\n");
  
    // open a file in read mode.
    ifstream infile;
    infile.open(file1); 
//------------------------------------------------------------------------------------------//    
    // faces_region 
    faces_region.resize(ngrains*2);
    ntotal_faces = 0;
    ntotal_vertices = 0;
    max_faces=0;
    for (int i=0; i<ngrains; i++){ // Loop over the grains

        // current region ith
        infile >> reg; 
        
        // Number of faces in region ith
        infile >> nfaces; 
       
        faces_region[i*2] = reg;
        faces_region[i*2+1] = nfaces;

        ntotal_faces = ntotal_faces + nfaces; // Counter of faces (global)

        for (int j=0; j<nfaces; j++){ // Loop over the faces in the ith region   
            
            // Number of vertices in region ith
            infile >> nvertices;
            ntotal_vertices = ntotal_vertices + nvertices; // Counter of vertices (global)
        }

        infile.ignore();    
    }
    infile.close();

//------------------------------------------------------------------------------------------//
    // SUBREGIONS = [# region, Face_1, Face_2, ... , Face_n]

    // max_faces -> maximum number of faces in one region
    for (int i=0; i<ngrains; i++){
        if (faces_region[i*2+1]>max_faces){
            max_faces=faces_region[i*2+1];
        }
    }
    
    SUBREGIONS.resize(ngrains*(max_faces+1));
	SUBREGIONS_aux.resize(ngrains*(max_faces+1));
    //std::vector<int> SUBREGIONS_aux(ngrains*(max_faces+1));
    count_faces = 0;
    for (int i=0; i<ngrains; i++){
        nfaces = faces_region[i*2+1];
        reg = faces_region[i*2];
        SUBREGIONS[i*(max_faces+1)]=reg;
        SUBREGIONS_aux[reg*(max_faces+1)]=reg;
        for (int j=0; j<nfaces; j++){
            SUBREGIONS[i*(max_faces+1)+j+1] = count_faces;   
            SUBREGIONS_aux[reg*(max_faces+1)+j+1] = count_faces;           
            count_faces = count_faces + 1;  
        }
    }  

//------------------------------------------------------------------------------------------//
    // vertices_face
    infile.open(file1); 
    vertices_face.resize(ntotal_faces*2);
    count_faces = 0;
    for (int i=0; i<ngrains; i++){ // Loop over the grains
        
        // current region ith
        infile >> reg;  

        // Number of faces in region ith
        infile >> nfaces; 

        for (int j=0; j<nfaces; j++){ // Loop over the faces

            // Number of vertices per face
            infile >> nvertices; 

            vertices_face[count_faces*2+0] = count_faces;
            vertices_face[count_faces*2+1] = nvertices;
                        
            count_faces = count_faces + 1;
        }

    }
    infile.close();

//------------------------------------------------------------------------------------------//
    infile.open(file2);
    // FACES_struc = [#face, vertice_1, vertice_2, ... , vertice_n]
    ind_faces_ini = 0;
    ind_faces_fin = 0;
    count_vertices = 0;
    // max_vertices -> maximum number of vertices (global)
    max_vertices=0;
    for (int j=0; j<ntotal_faces; j++){
        if (vertices_face[j*2+1]>max_vertices){
            max_vertices=vertices_face[j*2+1];
        }
    }
    max_vertices_global = max_vertices;
    
    FACES_struc.resize(ntotal_faces*(max_vertices_global+1));
    VERTICES_struc.resize(ntotal_vertices*4);
    for (int i=0; i<ngrains; i++){ // Loop over the grains

        nfaces = faces_region[i*2+1]; // number of faces of ith region
        ind_faces_ini = ind_faces_fin;
        ind_faces_fin = ind_faces_fin + nfaces;
        // max_vertices -> maximum number of vertices of jth face
        max_vertices=0;
        for (int j=ind_faces_ini; j<ind_faces_fin; j++){
            if (vertices_face[j*2+1]>max_vertices){
                max_vertices=vertices_face[j*2+1];
            }
        }
        
        // FACES_struc_aux - > informations of faces in ith region
        std::vector<int> FACES_struc_aux;
        FACES_struc_aux.resize(nfaces*(max_vertices+1));
        for (int j=0; j<nfaces; j++){ // Loop over the faces of ith region
            nvertices = vertices_face[(ind_faces_ini+j)*2+1]; // number of vertices in jth face
            FACES_struc_aux[j*(max_vertices+1)] = SUBREGIONS[i*(max_faces+1)+j+1];
           
            for (int k=0; k<nvertices; k++){
                infile >> val;
                if(val != '(' || val!= ',' || val != ')'){
                    infile >> FACES_struc_aux[j*(max_vertices+1)+k+1];
                }           
            }
            infile.ignore();
        }
                
        // max_points -> maximum number of points of ith region
        max_points=0;
        for (int j=0; j<nfaces; j++){
            for (int k=0; k<max_vertices; k++){
                if (FACES_struc_aux[j*(max_vertices+1)+k+1]>max_points){
                    max_points=FACES_struc_aux[j*(max_vertices+1)+k+1];
                }
            }
        }
       
        // VERTICES_struc_aux - > informations of coordinates in ith region and jth face
        std::vector<double> VERTICES_struc_aux((max_points+1)*4);
        VERTICES_struc_aux.resize(nfaces*((max_points+1)*4));
        for (int j=0; j<max_points+1; j++){ // Loop over the points of jth face

            VERTICES_struc_aux[j*4] = j;
            for (int k=0; k<3; k++){
                infile >> val;
                if(val != '(' || val!= ',' || val != ')'){
                    infile >> VERTICES_struc_aux[j*4+k+1];
                }           
            }
            infile.ignore();
        }
        
        // ======================================================================================
        // =============================== FACES_struc ==========================================
        
        for (int j=0; j<nfaces; j++){
            nvertices = vertices_face[(ind_faces_ini+j)*2+1]; // number of vertices in jth face
            FACES_struc[(ind_faces_ini+j)*(max_vertices_global+1)] = ind_faces_ini+j;
            for (int k=0; k<nvertices; k++){
                FACES_struc[(ind_faces_ini+j)*(max_vertices_global+1)+k+1] = count_vertices;
                count_vertices = count_vertices + 1;
            }
        }  
    
        // ======================================================================================
        // =============================== VERTICES_struc =======================================

        for (int j=0; j<nfaces; j++){
            nvertices = vertices_face[(ind_faces_ini+j)*2+1]; // number of vertices in jth face
            for (int k=0; k<nvertices; k++){

                ind_1 = FACES_struc_aux[j*(max_vertices+1)+k+1];
                ind_2 = FACES_struc[(ind_faces_ini+j)*(max_vertices_global+1)+k+1];

                VERTICES_struc[ind_2*4] = ind_2;
                VERTICES_struc[ind_2*4+1] = VERTICES_struc_aux[ind_1*4+1];
                VERTICES_struc[ind_2*4+2] = VERTICES_struc_aux[ind_1*4+2];
                VERTICES_struc[ind_2*4+3] = VERTICES_struc_aux[ind_1*4+3];
            }
        }
        // ======================================================================================
        // ======================================================================================
    }
    infile.close();
//------------------------------------------------------------------------------------------//
    infile.open(file3);

    // ==========================================================================================
    // =============================== NORMAL_VECTORS ===========================================
    //std::vector<double> NORMAL_VECTORS(ntotal_vertices*4);
    NORMAL_VECTORS.resize(ntotal_vertices*4);
    for (int i=0; i<ntotal_faces; i++){ // Loop over the total faces
        NORMAL_VECTORS[i*4] = i;
        for (int j=0; j<3; j++){
            infile >> val;
            if(val != '(' || val!= ',' || val != ')'){
                infile >> NORMAL_VECTORS[i*4+j+1];
            }           
        }
        infile.ignore();
    }   

    // ==========================================================================================
    // ==========================================================================================

    infile.close();
//------------------------------------------------------------------------------------------//
    infile.open(file7);
    
    // ==========================================================================================
    // =================================== INTERFACES_struc =====================================
    int num,ind;
    std::vector<int> INTERFACES_aux(ngrains*2*(max_faces+1));
    std::vector<int> interfaces_region(ngrains);
    
    for (int i=0; i<ngrains; i++){ // Loop over the grains
        infile >> reg; // current region ith
        INTERFACES_aux[reg*2*(max_faces+1)] = reg; // current region ith
        nfaces = faces_region[i*2+1]; // faces in the ith region
        ind=0;
        for (int j=0; j<nfaces; j++){
            infile >> num;
            INTERFACES_aux[reg*2*(max_faces+1)+j+1] = num; 
            ind = ind + 1;  
        }
        interfaces_region[reg] = ind;
    }
    infile.close();

    std::vector<int> INTERFACES(2*ntotal_faces*5);
    int count_nint,nint_a,reg_a,reg_b,nint_b;
    count_nint=0;
    for (int i=0; i<ngrains; i++){ // Loop over the grains
        reg_a = i; // Region_A 
        nint_a =  interfaces_region[reg_a]; //interfaces in region a   
        for (int j=0; j<nint_a; j++){ // Loop over interfaces in region a    
            reg_b =  INTERFACES_aux[reg_a*2*(max_faces+1)+j+1];
            if (reg_b>=0){
                INTERFACES[count_nint*5] = count_nint;                  
                INTERFACES[count_nint*5+1] = reg_a; // Region a
                INTERFACES[count_nint*5+2] = reg_b; // Region b
                INTERFACES[count_nint*5+3] = SUBREGIONS_aux[reg_a*(max_faces+1)+j+1]; //face_a
                nint_b =  interfaces_region[reg_b]; //interfaces in region b
                for (int k=0; k<nint_b; k++){ // Loop over interfaces in region b
                    ind = INTERFACES_aux[reg_b*2*(max_faces+1)+k+1];       
                    if (ind==reg_a){
                        INTERFACES[count_nint*5+4] = SUBREGIONS_aux[reg_b*(max_faces+1)+k+1]; //face_b               
                    } 
                }
                count_nint = count_nint +1;
            }
        }
    }

    ntotal_interfaces= count_nint;
    INTERFACES_struc.resize(ntotal_interfaces*6);
    for (int i=0; i<ntotal_interfaces; i++){ // Loop over the grains
        for (int j=0; j<5; j++){
            INTERFACES_struc[i*6+j] = INTERFACES[i*5+j];
        }
    }

//------------------------------------------------------------------------------------------//
    infile.open(file9);

    // ==========================================================================================
    // =============================== Volumes ==================================================
    
    Volumes.resize(ngrains*2);
    for (int i=0; i<ngrains; i++){ // Loop over the total faces
        Volumes[i*2] = i;
        for (int j=0; j<1; j++){
            infile >> Volumes[i*2+j+1];            
        }
        infile.ignore();
    }   

    

    // ==========================================================================================
    // ==========================================================================================

//------------------------------------------------------------------------------------------//

}   
//==========================================================================================//
