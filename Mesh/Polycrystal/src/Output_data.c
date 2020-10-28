/*****************************************************************************/
/*                                                                           */
/*                   Module to write the output file.                        */
/*                                                                           */
/*****************************************************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "Global_variables.h"

void Output_file()
{
    
    int i, j;
    
    FILE *fp;
	
	fp = fopen("Export/Mesh.dat", "w+");
	
    fprintf(fp,"%10d\n",input.ngrains);

    fprintf(fp,"%10d %10d\n",output.ntotal_POINTS,3);

    for (i=0; i<output.ntotal_POINTS; i++){   
        for (j=0; j<4; j++){
            if (j>0){
                fprintf(fp,"%30.15f",output.POINTS[i*4+j]); 
            }   
        }
        fprintf(fp,"\n");
    }
	
    fprintf(fp,"%10d%10d\n",output.ntotal_ELEM,3);

    for (i=0; i<output.ntotal_ELEM; i++){   
        for (j=0; j<5; j++){
            if (j>1){
            	fprintf(fp,"%10d",output.ELEM[i*5+j]+1);    
        	}
		}
        fprintf(fp,"\n");
    }

    fprintf(fp,"%10d%10d\n",output.ntotal_ELEM,3);

    for (i=0; i<output.ntotal_ELEM; i++){   
        for (j=0; j<4; j++){
            if (j>0){
            	fprintf(fp,"%30.15f",input.NORMAL_VECTORS[i*4+j]);
		}
        }
        fprintf(fp,"\n");
    }

    fprintf(fp,"%10d%10d\n",input.ngrains,1);

    for (i=0; i<input.ngrains; i++){   
        for (j=1; j<2; j++){
            //if (j>0){
                fprintf(fp,"%10d",output.El_reg[i*2+j]);
            //}
        }
        fprintf(fp,"\n");
    }

    //fprintf(fp,"%10d%10d\n",input.ngrains,2);

    //for (i=0; i<input.ngrains; i++){   
    //    fprintf(fp,"%10.15f             ",output.Volumes[i*2+1]);
    //    fprintf(fp,"\n");
    //}
    

    //fprintf(fp,"%10d%10d\n",output.ntotal_SEGMENTS,3);

    //for (i=0; i<output.ntotal_SEGMENTS; i++){   
    //   for (j=0; j<3; j++){
    //        fprintf(fp,"%10d",output.SEGMENTS[i*3+j]+1);    
    //    }
    //    fprintf(fp,"\n");
    //}

    fprintf(fp,"end");

    fclose(fp);    
}

void Output_file3()
{

    int i, j, type_el;
	
    FILE *fp;
	
    fp = fopen("Export/Mesh.vtk", "w+");
	
	fprintf(fp,"# vtk DataFile Version 2.3\n");
	fprintf(fp,"Mesh.vtk\n");
	fprintf(fp,"ASCII\n");
	fprintf(fp,"\n");
	fprintf(fp,"DATASET UNSTRUCTURED_GRID\n");

	fprintf(fp,"POINTS%20d double\n",output.ntotal_POINTS);

   	
	for (i=0; i<output.ntotal_POINTS; i++){   
        for (j=0; j<4; j++){
            if (j>0){
                fprintf(fp,"%27.16f",output.POINTS[i*4+j]); 
            }   
        }
        fprintf(fp,"\n");
    }

	fprintf(fp,"CELLS%20d%20d\n",output.ntotal_ELEM,output.ntotal_ELEM*4);

	type_el = 3;
    for (i=0; i<output.ntotal_ELEM; i++){
		fprintf(fp,"%10d",type_el);   
        for (j=0; j<5; j++){
            if (j>1){
            	fprintf(fp,"%10d",output.ELEM[i*5+j]);    
        	}
		}
        fprintf(fp,"\n");
    }

	fprintf(fp,"CELL_TYPES%20d\n",output.ntotal_ELEM);
	
	type_el = 5;
    for (i=0; i<output.ntotal_ELEM; i++){
		fprintf(fp,"%10d",type_el);   
        
        fprintf(fp,"\n");
    }


    fclose(fp);   

}

