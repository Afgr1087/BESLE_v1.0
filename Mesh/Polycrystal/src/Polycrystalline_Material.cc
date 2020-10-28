

#include "voro++.hh"
#include <time.h>
#include <stdio.h>
#include "Polycrystalline_Material.hh"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

using namespace voro;
using namespace std;

void start_greet()
    {

    printf("\n'****************************************************************************'\n");
    printf("'****************************************************************************'\n\n");
    printf("                3D MESH GENERATOR OF POLYCRYSTALLINE MATERIALS \n\n");
    printf("    3D Polycrystalline structure triangular mesh using the libreries \n");
    printf("    Voro++ and the Triangle mesh generator. \n\n");
    printf("    By Andres F. Galvis \n");
    printf("    School of Mechanical Engineering \n");
    printf("    University of Campinas \n");
    printf("    Date: 19/02/2020 \n\n");
    printf("'****************************************************************************'\n");
    printf("'****************************************************************************'\n\n");

    }

// Set up constants for the container geometry
const double x_min=0;
const double y_min=0;
const double z_min=0;

// This function returns a random double between 0 and 1
double rnd() {return double(rand())/RAND_MAX;}

void Polycrystalline_structure()
    {
        // Set the number of particles that are going to be randomly introduced
        //const double cvol=(x_max-x_min)*(y_max-y_min)*(x_max-x_min);
	    int i;
	    double x,y,z;
        int n_x,n_y,n_z;
        double x_min_cube, y_min_cube, z_min_cube;
        double x_max_cube, y_max_cube, z_max_cube;
        double x_min_cube_aux, y_min_cube_aux;
        double x_max_cube_aux, y_max_cube_aux;
        double space_x, space_y, space_z;
        int k, j, cont;
        
	Mesh_file="Export/Mesh.dat";

        n_x=ngrains_x; n_y=ngrains_y, n_z=ngrains_z;

        ngrains = ngrains_x*ngrains_y*ngrains_z;

	    // Create a container with the geometry given above, and make it
	    // non-periodic in each of the three coordinates. Allocate space for
	    // eight particles within each computational block
	    container con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
			    false,false,false,8);

	    //Randomly add particles into the container
        // srand (time(NULL));
	    //for(i=0;i<ngrains;i++) {
		//    x=x_min+rnd()*(x_max-x_min);
		//    y=y_min+rnd()*(y_max-y_min);
		//    z=z_min+rnd()*(z_max-z_min);
		//    con.put(i,x,y,z);
	    //}
     
        space_x = x_max/ngrains_x;
        space_y = y_max/ngrains_y;
        space_z = z_max/ngrains_z;

        x_min_cube = 0;
        y_min_cube = 0;
        z_min_cube = 0;

        x_max_cube = space_x;
        y_max_cube = space_y;
        z_max_cube = space_z;

        x_min_cube_aux = x_min_cube;
        y_min_cube_aux = y_min_cube;
        
        x_max_cube_aux = x_max_cube;
        y_max_cube_aux = y_max_cube;
                

        // open a file in write mode.
        ofstream outfile;
        outfile.open(file8); 

        cont = 0;

        if (stochastic == 1){
            srand (time(NULL));
        }
        for(k=0;k<ngrains_z;k++){
            if(cont>0){
                z_min_cube = z_min_cube + space_z;
                z_max_cube = z_max_cube + space_z;
            }
            for(j=0;j<ngrains_y;j++){
                if(j>0){
                    y_min_cube_aux = y_min_cube_aux + space_y;
                    y_max_cube_aux = y_max_cube_aux + space_y;
                }
                
                for(i=0;i<ngrains_x;i++){
                    if(i>0){    
                        x_min_cube_aux = x_min_cube_aux + space_x;    
                        x_max_cube_aux = x_max_cube_aux + space_x;    
                    }
                    
                    
                    x=x_min_cube_aux+rnd()*(x_max_cube_aux-x_min_cube_aux);
		            y=y_min_cube_aux+rnd()*(y_max_cube_aux-y_min_cube_aux);
		            z=z_min_cube+rnd()*(z_max_cube-z_min_cube);
                   
                    outfile << cont << ' ' << x <<  ' ' << y << ' ' << z << endl;
                    
                    //con.put(cont,x,y,z);

                    cont = cont + 1;

                }   
                x_min_cube_aux = x_min_cube;
                x_max_cube_aux = x_max_cube;
            }
            y_min_cube_aux = y_min_cube;
            y_max_cube_aux = y_max_cube;
        }

        con.import(file8);
        
        // Do a custom output routine to store the number of vertices, edges,
        // and faces of each Voronoi cell
	con.print_custom("%i", file0);
        con.print_custom("%i %s %a", file1);
        con.print_custom("%t\n%P", file2);
        con.print_custom("%l", file3);
        con.print_custom("%i\n%n", file7);
        con.print_custom("%v", file9);
        

	    // Output the particle positions in gnuplot format
	    con.draw_particles(file5);

	    // Output the Voronoi cells in gnuplot format
	    con.draw_cells_gnuplot(file6);

        
    }

