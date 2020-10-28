struct external_data {
    
    //input
    
    int ngrains;                
    int *SUBREGIONS; 
    int *faces_region;                                                                             
    int ntotal_faces; 
    int *FACES_struc; 
    int *vertices_face;
    int ntotal_vertices;   
    float *VERTICES_struc;
    float *NORMAL_VECTORS;
    int ntotal_interfaces;
    int *INTERFACES_struc;    
    float *VERTICES_local_struc; 
    int *divisions_segment;    
    int *SEGMENTS_struc;
    float *POINTS_struc;
    int *ELEM_struc;
    int ntotal_points;
    int *divisions_face;  
    float *average_face; 

    //output
	char Mesh_file;
    int *Store_ELEM;
    float *Store_POINTS; 
    int *FACES;
    int *El_reg;
    int ntotal_ELEM;
    int ntotal_POINTS;    
    int *ELEM;
    float *POINTS;
    int ntotal_SEGMENTS;
    int *store_SEGMENTS;
    int *SEGMENTS;
    float *Volumes;
};
extern struct external_data input, output;
extern int max_faces, max_vertices_global, dm;
extern double Amax;


