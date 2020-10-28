!====================================================================================!
MODULE Set_parameters
!-------------------------------------------------------------------------------!
USE Global_variables
!-------------------------------------------------------------------------------!	
CONTAINS

	SUBROUTINE Setup

	fileplace = "Data/Anisotropic/"
	file_name = "Material_8.dat"

	Material = "Anisotropic" 
	Lattice = "hcp" 

	C11=165000.d0; C12=31100.d0; C13=50000.d0; C33=61800.d0; C44=33600.d0	
	
	theta_x=-125.d0; theta_y= 10.d0; theta_z=332.d0;

	END SUBROUTINE Setup

END MODULE Set_parameters
!====================================================================================!
