//==========================================================================================//
//------------------------------------------------------------------------------------------//
//                        MODULE TO GENERATE THE TRIANGLES ELEMENTS                         //
//------------------------------------------------------------------------------------------//
//==========================================================================================//

#include <vector>
#include "Ouput_data.hh"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;

void Ouput_data_Triangle_Mesh() 
{
 
    // open a file in write mode.
    ofstream outfile;
    outfile.open(file4); 
    
	// Mesh name
	outfile << Mesh_file << endl;

    // Mesh density parameter
    outfile << dm << endl;

    // Write matrix of SUBREGIONS and faces_region
    outfile << ngrains << ' ' << max_faces<< endl;
    
    for (int i=0; i<ngrains; i++){    
        for (int j=0; j<max_faces+1; j++){
            outfile << SUBREGIONS[i*(max_faces+1)+j] << ' ';
        }
        outfile << ' ' << endl;
    }
    for (int i=0; i<ngrains; i++){    
        for (int j=0; j<2; j++){
            outfile << faces_region[i*2+j] << ' ';
        }
        outfile << ' ' << endl;
    }

    // Write matrix of FACES_struc and vertices_faces
    outfile << ntotal_faces << ' ' << max_vertices_global << endl;
    for (int i=0; i<ntotal_faces; i++){    
        for (int j=0; j<max_vertices_global+1; j++){
            outfile << FACES_struc[i*(max_vertices_global+1)+j] << ' ';
        }
        outfile << ' ' << endl;
    }
    for (int i=0; i<ntotal_faces; i++){    
        for (int j=0; j<2; j++){
            outfile << vertices_face[i*2+j] << ' ';
        }
        outfile << ' ' << endl;
    }

    // Write matrix of VERTICES_struc
    outfile << ntotal_vertices << endl;
    for (int i=0; i<ntotal_vertices; i++){    
        for (int j=0; j<4; j++){
            outfile << VERTICES_struc[i*4+j] << ' ';
        }
        outfile << ' ' << endl;
    }

    // Write matrix of NORMAL_VECTORS
    for (int i=0; i<ntotal_faces; i++){    
        for (int j=0; j<4; j++){
            outfile << NORMAL_VECTORS[i*4+j] << ' ';
        }
        outfile << ' ' << endl;
    }

    // Write matrix of INTERFACES_struc
    outfile << ntotal_interfaces << endl;
    for (int i=0; i<ntotal_interfaces; i++){    
        for (int j=0; j<6; j++){
            outfile << INTERFACES_struc[i*6+j] << ' ';
        }
        outfile << ' ' << endl;
    }

       
    // Write matrix of Volumes of each grain
    for (int i=0; i<ngrains; i++){    
        for (int j=0; j<2; j++){
            outfile << Volumes[i*2+j] << ' ';
        }
        outfile << ' ' << endl;
    }

    //------------------------------------------------------------------------------------------//

    //printf("\n'----------------------------------------------------------------------------'\n");

}

