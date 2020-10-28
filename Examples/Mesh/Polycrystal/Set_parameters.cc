#include "Set_parameters.hh"
//====================================== SETUP =========================================
//--------------------------------------------------------------------------------------
void Setup()
    {

		// Set up the number of blocks that the container is divided into
		ngrains_x=3;
		ngrains_y=4;
		ngrains_z=8;

		// Box dimensions
		x_max = 15.0;
		y_max = 15.0;
		z_max = 45.0;

		// Random geometry
		stochastic = 0;

		// Mesh density parameter
		dm = 1;
    }

//--------------------------------------------------------------------------------------
//======================================================================================
