
#include <vector>

int ngrains;

// All geometrical data of the polycrystalline structure
const char* file0="IO_files/packing0.dat";
const char* file1="IO_files/packing1.dat";
const char* file2="IO_files/packing2.dat";
const char* file3="IO_files/packing3.dat";
const char* file7="IO_files/packing5.dat";
const char* file8="IO_files/packing6.dat";
const char* file9="IO_files/packing7.dat";
// Write data to the triangle mesh generator
const char* file4="IO_files/packing4.dat";
// Points and vertices to plot the structure in gnuplot
// splot "plot_p.gnu" u 2:3:4, "plot_v.gnu" with lines
const char* file5="IO_files/plot_p.gnu";
const char* file6="IO_files/plot_v.gnu";

// Matrices fomated 
int max_faces, ntotal_vertices, ntotal_faces, max_vertices_global, ntotal_interfaces;
std::vector<int> faces_region;
std::vector<int> SUBREGIONS;
std::vector<int> SUBREGIONS_aux;
std::vector<int> vertices_face;
std::vector<int> vertices_face_aux;
std::vector<int> FACES_struc;
std::vector<double> VERTICES_struc;
std::vector<double> NORMAL_VECTORS; 
std::vector<int> INTERFACES_struc;
std::vector<double> Volumes;
