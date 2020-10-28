!====================================================================================!
MODULE Set_parameters
!-------------------------------------------------------------------------------!
USE Global_variables
!-------------------------------------------------------------------------------!	
CONTAINS

	SUBROUTINE Setup

	!n_materials = 100

	fileplace = "Data/Anisotropic/"
	file_name = "Material_1.dat"

	Material = "Anisotropic" 
	Lattice = "cubic" 

	!C11=165000.d0; C12=31000.1d0; C13=50000.d0; C33=61000.8d0; C44=33000.6d0	
	C11=230000.d0; C12=135000.d0; C44=117000.d0;	

	!E = 210000.d0; ni = 0.3;

	END SUBROUTINE Setup

END MODULE Set_parameters
!====================================================================================!
