#include <vector>

// Number of grains
extern int ngrains;
extern int dm; 


//Matrices in correct format
extern int max_faces, ntotal_vertices, ntotal_faces, max_vertices_global, ntotal_interfaces;
extern std::vector<int> faces_region;
extern std::vector<int> SUBREGIONS;
extern std::vector<int> SUBREGIONS_aux;
extern std::vector<int> vertices_face;
extern std::vector<int> FACES_struc;
extern std::vector<double> VERTICES_struc;
extern std::vector<double> NORMAL_VECTORS;
extern std::vector<int> INTERFACES_struc;
extern std::vector<double> Volumes;
// Write data to the triangle mesh generator
extern const char* file4;
extern const char* Mesh_file;
