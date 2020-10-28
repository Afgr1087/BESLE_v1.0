!====================================================================================!
MODULE Set_parameters
!-------------------------------------------------------------------------------!
USE Global_variables
!-------------------------------------------------------------------------------!	
CONTAINS

	SUBROUTINE Setup

	fileplace = "Data/Anisotropic/"
	file_name = "Material_13.dat"

	Material = "Anisotropic"
	Lattice = "triagonal" 

	C11=87600.d0; C12=6070.d0; 
        C13=13300.d0; C14=17300.d0; 
        C33=106800.d0; C44=57200.d0;
	
        z_x_z = 1
	theta_z=-15.d0; theta_x=135.d0;

	END SUBROUTINE Setup

END MODULE Set_parameters
!====================================================================================!
