//====================================================================================//
//====================================================================================//
//                                                                                    //
//               U N I C A M P - UNIVERSITY OF CAMPINAS                               //
//                     F E M - SCHOOL OF MECHANICAL ENGINEERING                       //
//                  D M C - DEPARTMENT OF COMPUTATIONAL MECHANICS                     //
//                                                                                    //
// AUTHOR: Andrés Felipe Galvis Rodriguez                                             //
// SUPERVISOR: Prof. Dr. Paulo Sollero                                                //
//************************************************************************************//
// Program: 3D_Polycrystalline                                                        //
// Description: 3D Polycrystalline structure triangular mesh                          //
//              using Voro++[1] and Triangle mesh[2] generator                        //
//                                                                                    //
// [1] C. H. Rycroft, “Voro++: A three-dimensional voronoi cell library in c++,”      //
//     Chaos, vol. 19, p. 041111, 2009.                                               //
// [2] J. R. Shewchuk, Triangle: Engineering a 2D Quality Mesh Generator and          //
//      Delaunay Triangulator. Springer-Verlag, 1996.                                 //
//====================================================================================//
//====================================================================================//

#include <stdio.h>
#include <fstream>

#include "Global_variables.hh"
using namespace std;
extern void Polycrystalline_structure();
extern void start_greet();
extern void Setup();
extern void Input_data_structure();
extern void Ouput_data_Triangle_Mesh();

//--------------------------------------------------------------------------------------------------------
//========================================================================================================
//========================================== MAIN PROGRAM =================================================
//---------------------------------------------------------------------------------------------------------
int main()
{
	int reg, i, reg_0, val, iteration=0;

    start_greet();

	Setup();

	printf("    1. Generation of polycrystalline structure \n\n");
	while (iteration==0){

		// 1. Generation of polycrystalline structure 
		Polycrystalline_structure();
		
		// check if voro++ gave grains in numerical order
		ifstream infile;
		infile.open(file0);
		std::vector<int> regs;
		regs.resize(ngrains);
		for (i=0; i<ngrains; i++){ // Loop over the grains
		    infile >> reg; // current region ith
			regs[i] = reg;
		}
		infile.close();

		reg_0 = regs[0]; // region
		//printf("reg_0= %d\n",reg_0);
		for (i=1; i<ngrains; i++){
			val = regs[i] - reg_0;
			//printf("val=%d ",regs[i]);
			if (val != 1){
				//printf("break reg=%d, regs=%d\n",reg_0, regs[i]);
				iteration=0;
				break;
			}else{
				if (i==ngrains-1){
					iteration=1;
				}
				reg_0 = regs[i];
			}
		}
		//printf("iteration = %d",iteration);
	
	}

	printf("'----------------------------------------------------------------------------'\n");

    // 2. Data adquisition from the polycrystalline structure
    Input_data_structure();
    
    // 3. Export data to triangle mesh C program
    Ouput_data_Triangle_Mesh();

    return 0;
}
//---------------------------------------------------------------------------------------------------------
//========================================== END PROGRAM ==================================================
