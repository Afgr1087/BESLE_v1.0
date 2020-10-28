!====================================================================================!
MODULE Set_parameters
!-------------------------------------------------------------------------------!
USE Global_variables
!-------------------------------------------------------------------------------!	
CONTAINS

	SUBROUTINE Setup

	fileplace = "Data/Isotropic/"
	file_name = "Material_7.dat"

	Material = "Isotropic" 	

	E = 210000.d0; ni = 0.3;

	END SUBROUTINE Setup

END MODULE Set_parameters
!====================================================================================!
