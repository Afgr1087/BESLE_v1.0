
const char* Mesh_file;

// Number of grains
extern int ngrains; 

// Box dimensions
extern double x_max;
extern double y_max;
extern double z_max;

extern int ngrains_x;
extern int ngrains_y;
extern int ngrains_z;

extern int stochastic;

// All geometrical data of the polycrystalline structure
extern const char* file0;
extern const char* file1;
extern const char* file2;
extern const char* file3;
extern const char* file7;
extern const char* file8;
extern const char* file9;
// Points and vertices to plot the structure in gnuplot
// splot "plot_p.gnu" u 2:3:4, "plot_v.gnu" with lines
extern const char* file5;
extern const char* file6;


