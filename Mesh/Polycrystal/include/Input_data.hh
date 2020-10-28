#include <vector>


// Number of grains
extern int ngrains; 

// All geometrical data of the polycrystalline structure
extern const char* file0;
extern const char* file1;
extern const char* file2;
extern const char* file3;
extern const char* file7;
extern const char* file9;
// Matrices in correct format
extern int max_faces, ntotal_vertices, ntotal_faces, max_vertices_global,ntotal_interfaces;
extern std::vector<int> faces_region;
extern std::vector<int> SUBREGIONS;
extern std::vector<int> vertices_face;
extern std::vector<int> SUBREGIONS_aux;
extern std::vector<int> FACES_struc;
extern std::vector<double> VERTICES_struc;
extern std::vector<double> NORMAL_VECTORS;
extern std::vector<int> INTERFACES_struc; 
extern std::vector<double> Volumes;
