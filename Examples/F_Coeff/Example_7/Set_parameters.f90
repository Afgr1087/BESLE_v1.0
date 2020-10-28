!====================================================================================!
MODULE Set_parameters
!-------------------------------------------------------------------------------!
USE Global_variables
!-------------------------------------------------------------------------------!	
CONTAINS

	SUBROUTINE Setup

	n_materials = 100

	fileplace = "Data/Anisotropic/Orientations/Cu/"
	file_name = "Angles.dat"

	Material = "Anisotropic" 
	Lattice = "cubic" 

	C11=168000.d0; C12=121400.d0; C44=75400.d0	

        z_x_z = 1

	END SUBROUTINE Setup

END MODULE Set_parameters
!====================================================================================!
